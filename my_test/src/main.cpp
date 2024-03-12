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
};

using namespace ncorr;

SubsetCorr CalculateSubsetDisp (const DIC_3D_Input& dic_in, const std::string& ref_image, 
                                const std::string& def_image, std::array<size_t, 2> s_xy_min) {

    std::vector<Image2D> imgs;
    imgs.push_back(ref_image);
    imgs.push_back(def_image);

    size_t img_size = 2048;
    auto ROI = Array2D<double>(img_size, img_size, 0.0);
    for (size_t j = 0; j < dic_in.subset_size; ++j) {
        for (size_t i = 0; i < dic_in.subset_size; ++i) {
            ROI(s_xy_min[0] + i, s_xy_min[1] + j) = 1.0;
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
    
    size_t n = 0;
    for (size_t i = 0; i < arr_u.size(); ++i) {
        if (arr_u(i) != 0) {
            u_aver += arr_u(i);
            v_aver += arr_v(i);
            ++n;
        }
    }

    if (n == 0) {
        return res;
    }

    res.u = u_aver /= n;
    res.v = v_aver /= n;
    res.coef = DIC_output.corr_coef;

    return res;
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
		throw std::invalid_argument("Must have 1 arg: DIC setting file!");	
	}

	DIC_3D_Input dic_in(argv[1]);
	dic_in.debug_print(std::cout);

    float coef_treshold = 0.15;
    float delta_coef_treshold = 0.01;

    size_t s_x = (dic_in.roi_xy_max[0] - dic_in.roi_xy_min[0] - dic_in.subset_size) / (dic_in.subset_offset) + 1;
    size_t s_y = (dic_in.roi_xy_max[1] - dic_in.roi_xy_min[1] - dic_in.subset_size) / (dic_in.subset_offset) + 1;
    size_t s_z = dic_in.stack_h / dic_in.z_bounce + 1;

    size_t subset_number = s_x * s_y;
    std::vector<int> interesting_layers = range(0, dic_in.stack_h + 1, dic_in.z_bounce);
    std::vector<std::vector<float>> result_ref(s_z * subset_number, std::vector<float>(3, 0.));
    auto result_def = result_ref;
    std::vector<float> result_coefs(s_z * subset_number, size_t(-1));

    for (size_t idx = 0; idx < dic_in.times.size() - 1; idx++) {
        auto initial_z_guess = interesting_layers;
        size_t ref_image_stack_number = dic_in.times[idx];
        size_t def_image_stack_number = dic_in.times[idx + 1];
        
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
                    ref_z_str << std::setw(3) << std::setfill('0') << std::setw(3) << k;
                    auto ref_image_str = dic_in.images_folder + dic_in.image_name_prefix + std::to_string(ref_image_stack_number) + dic_in.image_name_postfix + ref_z_str.str() + dic_in.image_extension;

                    size_t k_bounced = k / dic_in.z_bounce;

                    result_ref[xyz_table_start + k_bounced] = {subset_center[0], subset_center[1], k};

                    int init_z = initial_z_guess[k_bounced];
                    int search_z_min = k - dic_in.z_radius;
                    if (search_z_min < 0)
                        search_z_min = 0;

                    int search_z_max = k + dic_in.z_radius;
                    if (search_z_max >= dic_in.stack_h)
                        search_z_max = dic_in.stack_h - 1;

                    float best_coef = 1e6;
                    float best_u = 0., best_v = 0.;
                    int best_z = init_z;

                    for (int z = init_z; z < search_z_max + 1; ++z) {
                        std::stringstream def_z_str;
                        def_z_str << std::setw(3) << std::setfill('0') << std::setw(3) << z;
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
                            
                    }

                    if ((best_z > init_z || (init_z - 1) < 0) & (best_coef < coef_treshold)) {
                        result_def[xyz_table_start + k_bounced] = {subset_center[0] + best_u, subset_center[1] + best_v, best_z};
                        result_coefs[xyz_table_start + k_bounced] = best_coef;
                        initial_z_guess[k_bounced] = best_z;
                        continue;
                    }

                    for (int z = init_z - 1; z >= search_z_min + 1; --z) {
                        std::stringstream def_z_str;
                        def_z_str << std::setw(3) << std::setfill('0') << std::setw(3) << z;
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
                    }

                    result_def[xyz_table_start + k_bounced] = {subset_center[0] + best_u, subset_center[1] + best_v, best_z};
                    result_coefs[xyz_table_start + k_bounced] = best_coef;
                    initial_z_guess[k_bounced] = best_z;
                }

            }
        }

        std::ofstream out_file_ref(dic_in.output_folder + "ref_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
        if (out_file_ref.is_open()) {
            for (auto& xyz : result_ref) {
                out_file_ref << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;   
            }
        }
        out_file_ref.close();

        std::ofstream out_file_def(dic_in.output_folder + "def_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
        if (out_file_def.is_open()) {
            for (auto& xyz : result_def) {
                out_file_def << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;   
            }
        }
        out_file_def.close();

        std::ofstream out_file_coefs(dic_in.output_folder + "coefs_" + std::to_string(ref_image_stack_number) + "_" + std::to_string(def_image_stack_number) + ".txt");
        if (out_file_coefs.is_open()) {
            for (auto& coef : result_coefs) {
                out_file_coefs << coef << std::endl;   
            }
        }
        out_file_coefs.close();
    }

  	return 0;
}
