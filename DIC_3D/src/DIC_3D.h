#include "ncorr.h"
#include <iostream>
#include <nlohmann/json.hpp>

namespace DIC3D {

    enum class CalculationPathType {
        SnakeTopRight,
        FromRightToLeft,
        FromTopToBottom
    };

    struct CalcUnit {
        size_t x, y;
    };

    struct SubsetInfo {
        size_t init_x, init_y;
        float u, v, z;
        bool use_init;
        bool complete;
    };

    typedef std::vector<std::vector<std::vector<SubsetInfo>>> SubsetMap;
    std::vector<CalcUnit> GenerateCalculationPath(CalculationPathType path_type, size_t s_x, size_t s_y, size_t s_z, SubsetMap& subset_map);

    struct Input {

        Input(const std::string& filename) {
            std::ifstream f(filename, std::ios::in);
            if (!f.is_open()) {
                throw std::runtime_error("Can't open file " + filename + " !");
            }
            nlohmann::json settings_json = nlohmann::json::parse(f);
            settings_json.at("Image size x").get_to(image_size_xy[0]);
            settings_json.at("Image size y").get_to(image_size_xy[1]);
            settings_json.at("Subset size").get_to(subset_size);
            settings_json.at("Subset offset").get_to(subset_offset);
            settings_json.at("Z spacing").get_to(z_bounce);
            settings_json.at("Z search").get_to(z_radius);
            settings_json.at("Images path").get_to(images_folder);
            settings_json.at("DIC Results").get_to(output_folder);
            settings_json.at("Images Time Prefix").get_to(image_name_prefix);
            settings_json.at("Images Time Zeros").get_to(time_leading_zeros);
            settings_json.at("Images Slice Postfix").get_to(image_name_postfix);
            settings_json.at("Images Slice Zeros").get_to(slice_leading_zeros);
            settings_json.at("Times").get_to(times);
            settings_json.at("ROI x min").get_to(roi_xy_min[0]);
            settings_json.at("ROI y min").get_to(roi_xy_min[1]);
            settings_json.at("ROI x max").get_to(roi_xy_max[0]);
            settings_json.at("ROI y max").get_to(roi_xy_max[1]);
            settings_json.at("Downsampling").get_to(downsampling_factor);
            settings_json.at("Stack height").get_to(stack_h);
            settings_json.at("Image extension").get_to(image_extension);
            int is_2d, is_backward;
            settings_json.at("Is 2D").get_to(is_2d);
            settings_json.at("Is Backward").get_to(is_backward);
            settings_json.at("layers to calculate").get_to(layers_to_calculate);
            is_2D_case = bool(is_2d); backward_calculation = bool(is_backward);
            settings_json.at("coef threshold").get_to(coef_threshold);
            settings_json.at("delta coef threshold").get_to(delta_coef_threshold);
            settings_json.at("crack gradient threshold").get_to(gradient_threshold);

            int calc_path_type = 0;
            settings_json.at("Calculation Path").get_to(calc_path_type);
            if (calc_path_type == 0) {
                path_type = CalculationPathType::SnakeTopRight;
            }
            else if (calc_path_type == 1) {
                path_type = CalculationPathType::FromRightToLeft;
            }
            else {
                path_type = CalculationPathType::FromTopToBottom;
            }

            f.close();

            if (is_2D_case) {
                stack_h = 0;
                z_bounce = 1;
                z_radius = 0;
                image_name_postfix = "";
                slice_leading_zeros = 0;
            }
        };

        void debug_print(std::ostream& os) const {
            os << "Images size x: " << image_size_xy[0] << ", y: " << image_size_xy[1] << std::endl;
            os << "subset size: " << subset_size << ", " << "subset offset: " << subset_offset << std::endl;
            os << "z bounce: " << z_bounce << ", " << "z radius: " << z_radius << std::endl;
            os << "images folder: " << images_folder << std::endl;
            os << "outputmfolder: " << output_folder << std::endl;
            os << "image name prefix: " << image_name_prefix << " ,image name postfix: " << image_name_postfix << std::endl;
            os << "time leading zeros: " << time_leading_zeros << " ,slice leading zeros: " << slice_leading_zeros << std::endl;
            os << "times: ";
            for (size_t i = 0; i < times.size(); ++i) {
                os << times[i] << " ";
            }
            os << "\n" << "roi min: " << roi_xy_min[0] << ", " << roi_xy_min[1] << " ,roi max: " << roi_xy_max[0] << ", " << roi_xy_max[1] << std::endl;
            os << "downsampling factor: " << downsampling_factor << std::endl;
            os << "stack height: " << stack_h << std::endl;
            os << "image extension: " << image_extension << std::endl;
            os << "layers to calculate: ";
            for (size_t i = 0; i < layers_to_calculate.size(); ++i) {
                os << layers_to_calculate[i] << " ";
            }
            os << "\n" << "backward calculation: " << backward_calculation << std::endl;
            os << "2D case: " << is_2D_case << std::endl;
            os << "coef threshold: " << coef_threshold << std::endl;
            os << "delta coef threshold: " << delta_coef_threshold << std::endl;
            os << "crack gradient threshold: " << gradient_threshold << std::endl;
        }

        std::array<size_t, 2> image_size_xy;

        size_t subset_size, subset_offset;
        size_t z_bounce, z_radius;

        std::string images_folder;
        std::string output_folder;
        std::string image_name_prefix, image_name_postfix;
        size_t time_leading_zeros, slice_leading_zeros;
        std::vector<size_t> times;
        std::array<size_t, 2> roi_xy_min, roi_xy_max;
        size_t downsampling_factor;
        size_t stack_h;
        std::string image_extension;
        std::vector<size_t> layers_to_calculate;
        bool backward_calculation;
        bool is_2D_case;
        CalculationPathType path_type = CalculationPathType::FromRightToLeft;

        float coef_threshold, delta_coef_threshold, gradient_threshold;
    };

}