#include "ncorr.h"
#include <iostream>
#include <nlohmann/json.hpp>

struct DIC_3D_Input {

    DIC_3D_Input(const std::string& filename) {
        std::ifstream f(filename, std::ios::in);
        if (!f.is_open()) {
            throw std::runtime_error("Can't open file " + filename + " !");
        }
        nlohmann::json settings_json = nlohmann::json::parse(f);
        settings_json.at("Subset size").get_to(subset_size);
        settings_json.at("Subset offset").get_to(subset_offset);
        settings_json.at("Z spacing").get_to(z_bounce);
        settings_json.at("Z search").get_to(z_radius);
        settings_json.at("Images path").get_to(images_folder);
        settings_json.at("DIC Results").get_to(output_folder);
        settings_json.at("Images Time Prefix").get_to(image_name_prefix);
        settings_json.at("Images Slice Postfix").get_to(image_name_postfix);
        settings_json.at("Times").get_to(times);
        settings_json.at("ROI x min").get_to(roi_xy_min[0]);
        settings_json.at("ROI y min").get_to(roi_xy_min[1]);
        settings_json.at("ROI x max").get_to(roi_xy_max[0]);
        settings_json.at("ROI y max").get_to(roi_xy_max[1]);
        settings_json.at("Downsampling").get_to(downsampling_factor);
        settings_json.at("Stack height").get_to(stack_h);
        settings_json.at("Image extension").get_to(image_extension);
        int is_2d, is_backward, ignore1layer;
        settings_json.at("Is 2D").get_to(is_2d);
        settings_json.at("Is Backward").get_to(is_backward);
        settings_json.at("Ignore first layer").get_to(ignore1layer);
        is_2D_case = bool(is_2d); backward_calculation = bool(is_backward); ignore_1st_layer = bool(ignore1layer);

        f.close();

        if (is_2D_case) {
            ignore_1st_layer = false;
            stack_h = 1;
            z_bounce = 1;
            z_radius = 0;
            image_name_postfix = "";
        }
    };

    void debug_print(std::ostream& os) const {
        os << "subset size: " << subset_size << ", " << "subset offset: " << subset_offset << std::endl;
        os << "z bounce: " << z_bounce << ", " << "z radius: " << z_radius << std::endl;
        os << "images folder: " << images_folder << std::endl;
        os << "outputmfolder: " << output_folder << std::endl;
        os << "image name prefix: " << image_name_prefix << " ,image name postfix: " << image_name_postfix << std::endl;
        os << "times: ";
        for (size_t i = 0; i < times.size(); ++i) {
            os << times[i] << " ";
        }
        os << "\n" << "roi min: " << roi_xy_min[0] << ", " << roi_xy_min[1] << " ,roi max: " << roi_xy_max[0] << ", " << roi_xy_max[1] << std::endl;
        os << "downsampling factor: " << downsampling_factor << std::endl;
        os << "stack height: " << stack_h << std::endl;
        os << "image extension: " << image_extension << std::endl;
        os << "ignore 1st_layer: " << ignore_1st_layer << std::endl;
        os << "backward calculation: " << backward_calculation << std::endl;
        os << "2D case: " << is_2D_case << std::endl;
    }

    size_t subset_size, subset_offset;
    size_t z_bounce, z_radius;

    std::string images_folder;
    std::string output_folder;
    std::string image_name_prefix, image_name_postfix;
    std::vector<size_t> times;
    std::array<size_t, 2> roi_xy_min, roi_xy_max;
    size_t downsampling_factor;
    size_t stack_h;
    std::string image_extension;
    bool ignore_1st_layer;
    bool backward_calculation;
    bool is_2D_case;
};