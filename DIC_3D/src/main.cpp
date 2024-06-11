#include "ncorr.h"
#include "DIC_3D.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>

std::vector<int> range(int start, int stop, int step)
{
  if (step == 0)
  {
    throw std::invalid_argument("step for range must be non-zero");
  }

  std::vector<int> result;
  int i = start;
  while ((step > 0) ? (i < stop) : (i > stop))
  {
    result.push_back(i);
    i += step;
  }

  return result;
}

struct SubsetCorr {
    float u = 0., v = 0.;
    float coef = 1e6;
    bool skipped = false;
};

using namespace ncorr;

bool isValidPixel(int x, int y, int width, int height) {
  return (x >= 0 && x < width && y >= 0 && y < height);
}

// Function to calculate the gradient magnitude at a pixel
float calculateGradient(const std::vector<float>& image, int width, int height, int x, int y) {
  float gx = 0.0, gy = 0.0;

  const int SobelX[3][3] = {{ -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 }};
  const int SobelY[3][3] = {{ -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 }};

  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      int neighborX = x + i;
      int neighborY = y + j;
      if (isValidPixel(neighborX, neighborY, width, height)) {
        gx += SobelX[i + 1][j + 1] * image[neighborY * width + neighborX];
        gy += SobelY[i + 1][j + 1] * image[neighborY * width + neighborX];
      }
    }
  }

  return sqrt(gx * gx + gy * gy);
}

float calculateAverageGradient(const std::vector<float>& image, int width, int height,
                               int startX, int startY, int endX, int endY) {
  float sum = 0.0;
  int numPixels = 0;

  for (int y = startY; y < endY; ++y) {
    for (int x = startX; x < endX; ++x) {
      float gradient = calculateGradient(image, width, height, x, y);
      sum += gradient;
      numPixels++;
    }
  }

  if (numPixels == 0) {
    return 0.0;
  }

  return sum / numPixels;
}

SubsetCorr CalculateSubsetDisp (const DIC3D::Input& dic_in, const std::string& ref_image, 
                                const std::string& def_image, std::array<size_t, 2> s_xy_min, const std::pair<bool, std::vector<float>>& initial_guess_uv) {

    std::vector<Image2D> imgs;
    imgs.push_back(ref_image);
    imgs.push_back(def_image);

    auto ref_img_data = imgs[0].get_gs();
    const size_t img_size_x = dic_in.image_size_xy[0], img_size_y = dic_in.image_size_xy[1];
    std::vector<float>ref_vec(ref_img_data.begin(), ref_img_data.end());

    auto aver_grad = calculateAverageGradient(ref_vec, img_size_x, img_size_y, s_xy_min[0], s_xy_min[1], s_xy_min[0] + dic_in.subset_size, s_xy_min[1] + dic_in.subset_size);

    if (aver_grad < dic_in.gradient_threshold) {
        SubsetCorr res;
        res.skipped = true;
        return res;
    }

    auto ROI = Array2D<double>(img_size_y, img_size_x, 0.0);
    for (size_t j = 0; j < dic_in.subset_size; ++j) {
        for (size_t i = 0; i < dic_in.subset_size; ++i) {
            ROI(s_xy_min[1] + j, s_xy_min[0] + i) = 1.0;
        }
    }

    auto DIC_input = DIC_analysis_input(
                            imgs, 						  // Images
				            ROI2D(ROI > 0.5),		      // ROI
					       dic_in.downsampling_factor,    // scalefactor
					       INTERP::CUBIC_KEYS_PRECOMPUTE, // Interpolation
					       SUBREGION::SQUARE,			  // Subregion shape
					       dic_in.subset_size,            // Subregion radius
					       1,                                         		// # of threads
					       DIC_analysis_config::NO_UPDATE,// DIC configuration for reference image updates
					       false,// Debugging enabled/disabled
                           initial_guess_uv);			  

	// Perform DIC_analysis    
	auto DIC_output = DIC_analysis(DIC_input);
    SubsetCorr res;

    auto arr_u = DIC_output.disps[0].get_u().get_array();
    auto arr_v = DIC_output.disps[0].get_v().get_array();
    float u_aver = 0., v_aver = 0.;
    
    size_t n_u = 0, n_v = 0;
    for (size_t i = 0; i < arr_u.size(); ++i) {
        if (arr_u(i) != 0) {
            u_aver += arr_u(i);
            ++n_u;
        }
        if (arr_v(i) != 0) {
            v_aver += arr_v(i);
            ++n_v;
        }
    }

    if (n_u == 0) {
        return res;
    }

    res.u = u_aver / n_u;
    res.v = v_aver / n_v;
    res.coef = DIC_output.corr_coef;

    return res;
}

std::vector<std::vector<float>> Calculate_2D_disps(const DIC3D::Input& dic_in, const std::string& ref_image, 
                                const std::string& def_image, size_t s_x, size_t s_y)
{

    std::vector<Image2D> imgs;
    imgs.push_back(ref_image);
    imgs.push_back(def_image);

    auto ref_img_data = imgs[0].get_gs();
    const size_t img_size_x = dic_in.image_size_xy[0], img_size_y = dic_in.image_size_xy[1];
    std::vector<float>ref_vec(ref_img_data.begin(), ref_img_data.end());

    auto ROI = Array2D<double>(img_size_y, img_size_x, 0.0);
    for (size_t j = dic_in.roi_xy_min[1]; j < dic_in.roi_xy_max[1]; ++j) {
        for (size_t i = dic_in.roi_xy_min[0]; i < dic_in.roi_xy_max[0]; ++i) {
            ROI(j, i) = 1.0;
        }
    }

    auto DIC_input = DIC_analysis_input(
                        imgs, 						  // Images
				        ROI2D(ROI > 0.5),		      // ROI
					    dic_in.downsampling_factor,    // scalefactor
					    INTERP::CUBIC_KEYS_PRECOMPUTE, // Interpolation
					    SUBREGION::SQUARE,			  // Subregion shape
					    dic_in.subset_size,            // Subregion radius
					    1,                                         		// # of threads
					    DIC_analysis_config::NO_UPDATE,// DIC configuration for reference image updates
					    false,
                        std::make_pair<bool, std::vector<float>>(false, {0.0, 0.0}));//true);				  // Debugging enabled/disabled
    
    auto DIC_output = DIC_analysis(DIC_input);

    auto arr_u = DIC_output.disps[0].get_u().get_array();
    auto arr_v = DIC_output.disps[0].get_v().get_array();

    const size_t img_w = arr_u.width(), img_h = arr_u.height();
    if (arr_v.width() != img_w || arr_v.height() != img_h) {
        throw std::runtime_error("Wrong ncorr output!");
    }

    std::vector<std::vector<float>> averaged_disps(s_x * s_y, std::vector<float>(2, 0.));
    const size_t win_size = dic_in.subset_size / dic_in.downsampling_factor;
    for (size_t i = 0; i < s_y; ++i) {
        for (size_t j = 0; j < s_x; ++j) {
            size_t x_min = (dic_in.roi_xy_min[0] + dic_in.subset_offset * j);
            size_t y_min = (dic_in.roi_xy_min[1] + dic_in.subset_offset * i);
            float u_aver = 0.0, v_aver = 0.0;

            auto aver_grad = calculateAverageGradient(ref_vec, img_size_x, img_size_y, x_min, y_min, x_min + dic_in.subset_size, y_min + dic_in.subset_size);
            if (aver_grad < dic_in.gradient_threshold) {
                continue;
            }
            x_min /= dic_in.downsampling_factor;
            y_min /= dic_in.downsampling_factor;
            for (size_t y = y_min; y < y_min + win_size; ++y) {
                for (size_t x = x_min; x < x_min + win_size; ++x) {
                    u_aver += arr_u(y, x);
                    v_aver += arr_v(y, x);
                }
            }
            u_aver /= win_size * win_size;
            v_aver /= win_size * win_size;

            averaged_disps[i * s_y + s_x] = {u_aver, v_aver};
        }
    }

    return averaged_disps;
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
		throw std::invalid_argument("Must have 1 arg: DIC setting file!");	
	}

	DIC3D::Input dic_in(argv[1]);
	dic_in.debug_print(std::cout);

    const float bad_coef = 1e6;

    const size_t s_x = (dic_in.roi_xy_max[0] - dic_in.roi_xy_min[0] - dic_in.subset_size) / (dic_in.subset_offset) + 1;
    const size_t s_y = (dic_in.roi_xy_max[1] - dic_in.roi_xy_min[1] - dic_in.subset_size) / (dic_in.subset_offset) + 1;
    const size_t subset_number = s_x * s_y;

    auto stack_times = dic_in.times;

    if (dic_in.backward_calculation) {
        std::reverse(stack_times.begin(), stack_times.end());
    }

    size_t s_z = 1;
    std::vector<int> interesting_layers = {0};
    if (!dic_in.is_2D_case) {
        s_z = dic_in.stack_h / dic_in.z_bounce + 1;
        interesting_layers = range(0, dic_in.stack_h + 1, dic_in.z_bounce);
    }

    std::vector<std::thread> threads(stack_times.size() - 1);

    for (size_t idx = 0; idx < stack_times.size() - 1; idx++) {

        if (dic_in.is_2D_case) {
            threads[idx] = std::thread([&](size_t ind) {
                size_t ref_image_stack_number = stack_times[ind];
                size_t def_image_stack_number = stack_times[ind + 1];

                std::vector<std::vector<float>> result_ref(subset_number, std::vector<float>(2, 0.));

                for (size_t i = 0; i < s_y; ++i) {
                    for (size_t j = 0; j < s_x; ++j) {
                        result_ref[i * s_y + j] = {static_cast<float>(dic_in.roi_xy_min[0] + dic_in.subset_offset * j + dic_in.subset_size / 2),
                                                     static_cast<float>(dic_in.roi_xy_min[1] + dic_in.subset_offset * i + dic_in.subset_size / 2)};
                    }
                }
                std::stringstream ref_t_str;
                ref_t_str << std::setw(dic_in.time_leading_zeros) << std::setfill('0') << std::setw(dic_in.time_leading_zeros) << ref_image_stack_number;
                auto ref_image_str = dic_in.images_folder + dic_in.image_name_prefix + ref_t_str.str() + dic_in.image_name_postfix + dic_in.image_extension;
                std::stringstream def_t_str;
                def_t_str << std::setw(dic_in.time_leading_zeros) << std::setfill('0') << std::setw(dic_in.time_leading_zeros) << def_image_stack_number;
                auto def_image_str = dic_in.images_folder + dic_in.image_name_prefix + def_t_str.str() + dic_in.image_name_postfix + dic_in.image_extension;

                printf("Processing 2D correlation between %zu and %zu images...\n", ref_image_stack_number, def_image_stack_number);
                auto start_corr = std::chrono::system_clock::now();

                auto result_disps = Calculate_2D_disps(dic_in, ref_image_str, def_image_str, s_x, s_y);
                auto result_def = result_ref;

                for (size_t p = 0; p < result_disps.size(); ++p) {
                    const auto& disp = result_disps[p];
                    result_def[p][0] += disp[0];
                    result_def[p][1] += disp[1];
                }

                const auto end_corr = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds_2stacks = (end_corr - start_corr);
                printf("Time spent: %fs on stacks %zu and %zu.\n", elapsed_seconds_2stacks.count(), ref_image_stack_number, def_image_stack_number);

                std::ofstream out_file_ref(dic_in.output_folder + "ref_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
                if (!out_file_ref.is_open()) {
                    throw std::runtime_error("Can't open folder " + dic_in.output_folder + " !");
                }
                for (auto& xyz : result_ref) {
                    out_file_ref << xyz[0] << "," << xyz[1] << std::endl;   
                }
                out_file_ref.close();

                std::ofstream out_file_def(dic_in.output_folder + "def_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
                if (!out_file_def.is_open()) {
                    throw std::runtime_error("Can't open folder  " + dic_in.output_folder + " !");
                }
                for (auto& xyz : result_def) {
                    out_file_def << xyz[0] << "," << xyz[1] << std::endl;   
                }
                out_file_def.close();

            }, idx);
        }
        else {
            threads[idx] = std::thread([&](size_t ind) {
                const size_t ref_image_stack_number = stack_times[ind];
                std::stringstream ref_t_str;
                ref_t_str << std::setw(dic_in.time_leading_zeros) << std::setfill('0') << std::setw(dic_in.time_leading_zeros) << ref_image_stack_number;
                const size_t def_image_stack_number = stack_times[ind + 1];
                std::stringstream def_t_str;
                def_t_str << std::setw(dic_in.time_leading_zeros) << std::setfill('0') << std::setw(dic_in.time_leading_zeros) << def_image_stack_number;
                
                auto initial_z_guess = interesting_layers;
                std::vector<std::pair<bool, std::vector<float>>> initial_uv_guess(initial_z_guess.size(), {false, {0.0, 0.0}});

                std::vector<float> result_coefs(s_z * subset_number, size_t(-1));
                std::vector<std::vector<float>> result_ref(s_z * subset_number, std::vector<float>(3, 0.));
                auto result_def = result_ref;

                printf("Processing correlation between %zu and %zu stacks...\n", ref_image_stack_number, def_image_stack_number);
                const auto start_corr = std::chrono::system_clock::now();

                DIC3D::SubsetMap subset_map;
                auto calc_path = DIC3D::GenerateCalculationPath(dic_in.path_type, s_x, s_y, s_z, subset_map);
                
                for (const auto& s_loc : calc_path) {

                    std::array<size_t, 2> s_xy_min;
                    s_xy_min[0] = dic_in.roi_xy_min[0] + dic_in.subset_offset * s_loc.x;
                    s_xy_min[1] = dic_in.roi_xy_min[1] + dic_in.subset_offset * s_loc.y;

                    std::array<size_t, 2> subset_center;
                    subset_center[0] = s_xy_min[0] + dic_in.subset_size / 2;
                    subset_center[1] = s_xy_min[1] + dic_in.subset_size / 2;

                    size_t xyz_table_start = s_loc.y * (s_x * s_z) + s_loc.x * s_z;

                    for (auto k : interesting_layers) {
                        std::stringstream ref_z_str;
                        ref_z_str << std::setw(dic_in.slice_leading_zeros) << std::setfill('0') << std::setw(dic_in.slice_leading_zeros) << k;
                        
                        const auto ref_image_str = dic_in.images_folder + dic_in.image_name_prefix + ref_t_str.str() + dic_in.image_name_postfix + ref_z_str.str() + dic_in.image_extension;
                        const size_t k_bounced = k / dic_in.z_bounce;

                        auto& subset_info = subset_map[s_loc.x][s_loc.y][k_bounced];

                        result_ref[xyz_table_start + k_bounced] = {static_cast<float>(subset_center[0]), static_cast<float>(subset_center[1]), static_cast<float>(k)};
                        if (std::find(dic_in.layers_to_calculate.begin(), dic_in.layers_to_calculate.end(), k) == dic_in.layers_to_calculate.end()) {
                            result_def[xyz_table_start + k_bounced] = result_ref[xyz_table_start + k_bounced];
                            result_coefs[xyz_table_start + k_bounced] = bad_coef;
                            continue;
                        }

                        int init_z = k;
                        std::pair<bool, std::vector<float>> init_uv = std::make_pair<bool, std::vector<float>>(false, {0.0, 0.0});
                        if (subset_info.use_init) {
                            const auto& guess_subset = subset_map[subset_info.init_x][subset_info.init_y][k_bounced];
                            if (guess_subset.complete) {
                                init_z = guess_subset.z;
                                init_uv = std::make_pair<bool, std::vector<float>>(true, {guess_subset.v, guess_subset.u});
                            }
                        }

                        int search_z_min = k - dic_in.z_radius;
                        if (search_z_min < 0)
                            search_z_min = 0;

                        int search_z_max = k + dic_in.z_radius;
                        if (search_z_max > dic_in.stack_h)
                            search_z_max = dic_in.stack_h;

                        float best_coef = bad_coef;
                        float best_u = 0., best_v = 0.;
                        int best_z = init_z;

                        for (int z = init_z; z <= search_z_max; ++z) {
                            std::stringstream def_z_str;
                            def_z_str << std::setw(dic_in.slice_leading_zeros) << std::setfill('0') << std::setw(dic_in.slice_leading_zeros) << z;

                            auto def_image_str = dic_in.images_folder + dic_in.image_name_prefix + def_t_str.str() + dic_in.image_name_postfix + def_z_str.str() + dic_in.image_extension;
                            auto res = CalculateSubsetDisp(dic_in, ref_image_str, def_image_str, s_xy_min, init_uv);

                            if (res.coef <= best_coef) {
                                best_u = res.u;
                                best_v = res.v;
                                best_coef = res.coef;
                                best_z = z;
                            }
                            else if ((best_coef < dic_in.coef_threshold) && (res.coef - best_coef > dic_in.delta_coef_threshold)) {
                                break; 
                            }
                        }

                        if ((best_z > init_z || (init_z - 1) < 0) & (best_coef < dic_in.coef_threshold)) {
                            result_def[xyz_table_start + k_bounced] = {static_cast<float>(subset_center[0] + best_u), static_cast<float>(subset_center[1] + best_v), static_cast<float>(best_z)};
                            result_coefs[xyz_table_start + k_bounced] = best_coef;
                            {
                                subset_info.u = best_u;
                                subset_info.v = best_v;
                                subset_info.z = best_z;
                                subset_info.complete = true;
                            }
                            continue;
                        }

                        for (int z = init_z - 1; z >= search_z_min + 1; --z) {
                            std::stringstream def_z_str;
                            def_z_str << std::setw(dic_in.slice_leading_zeros) << std::setfill('0') << std::setw(dic_in.slice_leading_zeros) << z;

                            auto def_image_str = dic_in.images_folder + dic_in.image_name_prefix + def_t_str.str() + dic_in.image_name_postfix + def_z_str.str() + dic_in.image_extension;
                            auto res = CalculateSubsetDisp(dic_in, ref_image_str, def_image_str, s_xy_min, init_uv);

                            if (res.coef <= best_coef) {
                                best_u = res.u;
                                best_v = res.v;
                                best_coef = res.coef;
                                best_z = z;
                            }
                            else if ((best_coef < dic_in.coef_threshold) && (res.coef - best_coef > dic_in.delta_coef_threshold)) {
                                break;
                            }
                        }

                        result_def[xyz_table_start + k_bounced] = {static_cast<float>(subset_center[0] + best_u), static_cast<float>(subset_center[1] + best_v), static_cast<float>(best_z)};
                        result_coefs[xyz_table_start + k_bounced] = best_coef;
                        {
                            subset_info.u = best_u;
                            subset_info.v = best_v;
                            subset_info.z = best_z;
                            subset_info.complete = (best_coef < dic_in.coef_threshold);
                        }
                    }

                }

                const auto end_corr = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds_2stacks = (end_corr - start_corr);
                printf("Time spent: %fs on stacks %zu and %zu.\n", elapsed_seconds_2stacks.count(), ref_image_stack_number, def_image_stack_number);

                std::ofstream out_file_ref(dic_in.output_folder + "ref_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
                if (!out_file_ref.is_open()) {
                    throw std::runtime_error("Can't open folder " + dic_in.output_folder + " !");
                }
                for (auto& xyz : result_ref) {
                    out_file_ref << xyz[0] << "," << xyz[1] << "," << xyz[2] << std::endl;   
                }
                out_file_ref.close();

                std::ofstream out_file_def(dic_in.output_folder + "def_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
                if (!out_file_def.is_open()) {
                    throw std::runtime_error("Can't open folder  " + dic_in.output_folder + " !");
                }
                for (auto& xyz : result_def) {
                    out_file_def << xyz[0] << "," << xyz[1] << "," << xyz[2] << std::endl;   
                }
                out_file_def.close();

                std::ofstream out_file_coefs(dic_in.output_folder + "coefs_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
                if (!out_file_coefs.is_open()) {
                    throw std::runtime_error("Can't open folder " + dic_in.output_folder + " !");
                }
                for (auto& coef : result_coefs) {
                    out_file_coefs << coef << std::endl;
                }
                out_file_coefs.close();

            }, idx);
        }
    }

    for (auto& thread : threads) {
        thread.join();
    }

  	return 0;
}