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

SubsetCorr CalculateSubsetDisp (const DIC_3D_Input& dic_in, const std::string& ref_image, 
                                const std::string& def_image, std::array<size_t, 2> s_xy_min) {

    std::vector<Image2D> imgs;
    imgs.push_back(ref_image);
    imgs.push_back(def_image);

    auto ref_img_data = imgs[0].get_gs();
    const size_t img_size_x = dic_in.image_size_xy[0], img_size_y = dic_in.image_size_xy[1];
    std::vector<float>ref_vec(ref_img_data.begin(), ref_img_data.end());

    auto aver_grad = calculateAverageGradient(ref_vec, img_size_x, img_size_y, s_xy_min[0], s_xy_min[1], s_xy_min[0] + dic_in.subset_size, s_xy_min[1] + dic_in.subset_size);

    if (aver_grad < 0.01) {
        // std::cout << "Skip: " << "[" << s_xy_min[0] << ", " << s_xy_min[1] << "]" <<std::endl;
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
					       false);//true);				  // Debugging enabled/disabled

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

int main(int argc, char *argv[]) {

    if (argc != 2) {
		throw std::invalid_argument("Must have 1 arg: DIC setting file!");	
	}

	DIC_3D_Input dic_in(argv[1]);
	dic_in.debug_print(std::cout);

    const float coef_treshold = 0.15;
    const float delta_coef_treshold = 0.1 * coef_treshold;
    // float wrong_coef_treshold = 0.20; //0.30
    const float bad_coef = 1e6;

    const size_t s_x = (dic_in.roi_xy_max[0] - dic_in.roi_xy_min[0] - dic_in.subset_size) / (dic_in.subset_offset) + 1;
    const size_t s_y = (dic_in.roi_xy_max[1] - dic_in.roi_xy_min[1] - dic_in.subset_size) / (dic_in.subset_offset) + 1;
    const size_t subset_number = s_x * s_y;

    size_t s_z = 1;
    std::vector<int> interesting_layers = {0};
    if (!dic_in.is_2D_case) {
        s_z = dic_in.stack_h / dic_in.z_bounce + 1;
        interesting_layers = range(0, dic_in.stack_h + 1, dic_in.z_bounce);
    }

    auto stack_times = dic_in.times;

    if (dic_in.backward_calculation) {
        std::reverse(stack_times.begin(), stack_times.end());
    }

    std::vector<std::thread> threads(stack_times.size() - 1);

    for (size_t idx = 0; idx < stack_times.size() - 1; idx++) {

        threads[idx] = std::thread([&](size_t ind) {
            
            auto initial_z_guess = interesting_layers;
            size_t ref_image_stack_number = stack_times[ind];
            size_t def_image_stack_number = stack_times[ind + 1];

            std::vector<float> result_coefs(s_z * subset_number, size_t(-1));
            std::vector<std::vector<float>> result_ref(s_z * subset_number, std::vector<float>(3, 0.));
            auto result_def = result_ref;

            printf("Processing correlation between %d and %d stacks...\n", ref_image_stack_number, def_image_stack_number);
            std::chrono::time_point<std::chrono::system_clock> start_corr = std::chrono::system_clock::now();
            
            for (size_t i = 0; i < s_y; ++i) {
                for (size_t p = 0; p < s_x; ++p) {
                    size_t j = p;
                    if (i % 2 == 0)
                        j = s_x - p - 1;

                    std::array<size_t, 2> s_xy_min;
                    s_xy_min[0] = dic_in.roi_xy_min[0] + dic_in.subset_offset * j;
                    s_xy_min[1] = dic_in.roi_xy_min[1] + dic_in.subset_offset * i;

                    std::array<size_t, 2> subset_center;
                    subset_center[0] = s_xy_min[0] + dic_in.subset_size / 2;
                    subset_center[1] = s_xy_min[1] + dic_in.subset_size / 2;

                    size_t xyz_table_start = i * (s_x * s_z) + j * s_z;

                    for (auto k : interesting_layers) {
                        std::stringstream ref_z_str;
                        if (!dic_in.is_2D_case) {
                            ref_z_str << std::setw(3) << std::setfill('0') << std::setw(3) << k;
                        }
                        
                        auto ref_image_str = dic_in.images_folder + dic_in.image_name_prefix + std::to_string(ref_image_stack_number) + dic_in.image_name_postfix + ref_z_str.str() + dic_in.image_extension;
                        size_t k_bounced = k / dic_in.z_bounce;

                        result_ref[xyz_table_start + k_bounced] = {static_cast<float>(subset_center[0]), static_cast<float>(subset_center[1]), static_cast<float>(k)};
                        if (dic_in.ignore_1st_layer && (k == interesting_layers[0])) {
                            result_def[xyz_table_start + k_bounced] = result_ref[xyz_table_start + k_bounced];
                            result_coefs[xyz_table_start + k_bounced] = bad_coef;
                            continue;
                        }

                        int init_z = initial_z_guess[k_bounced];
                        int search_z_min = k - dic_in.z_radius;
                        if (search_z_min < 0)
                            search_z_min = 0;

                        int search_z_max = k + dic_in.z_radius;
                        if (search_z_max > dic_in.stack_h)
                            search_z_max = dic_in.stack_h; // TODO: -1 is not necessary

                        float best_coef = bad_coef;
                        float best_u = 0., best_v = 0.;
                        int best_z = init_z;

                        for (int z = init_z; z <= search_z_max; ++z) {
                            std::stringstream def_z_str;
                            if (!dic_in.is_2D_case) {
                                def_z_str << std::setw(3) << std::setfill('0') << std::setw(3) << z;
                            }
                            auto def_image_str = dic_in.images_folder + dic_in.image_name_prefix + std::to_string(def_image_stack_number) + dic_in.image_name_postfix + def_z_str.str() + dic_in.image_extension;
                            auto res = CalculateSubsetDisp(dic_in, ref_image_str, def_image_str, s_xy_min);

                            if (res.coef <= best_coef) {
                                best_u = res.u;
                                best_v = res.v;
                                best_coef = res.coef;
                                best_z = z;
                            }
                            else if ((best_coef < coef_treshold) && (res.coef - best_coef > delta_coef_treshold)) {
                                break; 
                            }
                            // else if (best_coef < wrong_coef_treshold) {
                            //     break;
                            // }
                                
                        }

                        if ((best_z > init_z || (init_z - 1) < 0) & (best_coef < coef_treshold)) {
                            result_def[xyz_table_start + k_bounced] = {static_cast<float>(subset_center[0] + best_u), static_cast<float>(subset_center[1] + best_v), static_cast<float>(best_z)};
                            result_coefs[xyz_table_start + k_bounced] = best_coef;
                            initial_z_guess[k_bounced] = best_z;
                            continue;
                        }

                        for (int z = init_z - 1; z >= search_z_min + 1; --z) {
                            std::stringstream def_z_str;
                            if (!dic_in.is_2D_case) {
                                def_z_str << std::setw(3) << std::setfill('0') << std::setw(3) << z;
                            }
                            auto def_image_str = dic_in.images_folder + dic_in.image_name_prefix + std::to_string(def_image_stack_number) + dic_in.image_name_postfix + def_z_str.str() + dic_in.image_extension;
                            auto res = CalculateSubsetDisp(dic_in, ref_image_str, def_image_str, s_xy_min);

                            if (res.coef <= best_coef) {
                                best_u = res.u;
                                best_v = res.v;
                                best_coef = res.coef;
                                best_z = z;
                            }
                            else if ((best_coef < coef_treshold) && (res.coef - best_coef > delta_coef_treshold)) {
                                break;
                            }
                            // else if (best_coef < wrong_coef_treshold) {
                            //     break;
                            // }
                        }

                        result_def[xyz_table_start + k_bounced] = {static_cast<float>(subset_center[0] + best_u), static_cast<float>(subset_center[1] + best_v), static_cast<float>(best_z)};
                        result_coefs[xyz_table_start + k_bounced] = best_coef;
                        initial_z_guess[k_bounced] = best_z;
                    }

                }
            }

            std::chrono::time_point<std::chrono::system_clock> end_corr = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds_2stacks = end_corr - start_corr;
            printf("Time spent: %f on stacks %d and %d.\n", elapsed_seconds_2stacks.count(), ref_image_stack_number, def_image_stack_number);

            std::ofstream out_file_ref(dic_in.output_folder + "ref_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
            if (!out_file_ref.is_open()) {
                throw std::runtime_error("Can't open folder " + dic_in.output_folder + " !");
            }
            for (auto& xyz : result_ref) {
                out_file_ref << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;   
            }
            out_file_ref.close();

            std::ofstream out_file_def(dic_in.output_folder + "def_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
            if (!out_file_def.is_open()) {
                throw std::runtime_error("Can't open folder  " + dic_in.output_folder + " !");
            }
            for (auto& xyz : result_def) {
                out_file_def << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;   
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

    for (size_t i = 0; i < threads.size(); ++i) {
        threads[i].join();
    }

  	return 0;
}