#include "ncorr.h"
#include <iostream>

struct DIC_3D_Input {

    DIC_3D_Input(const std::string& filename) {
        std::ifstream f(filename, std::ios::in);
        if (!f.is_open()) {
            throw std::runtime_error("Can't open file " + filename + " !");
        }
        // f >> is_2D_case;

        f >> subset_size >> subset_offset;
        f >> z_bounce >> z_radius;
        f >> images_folder;
        f >> output_folder;
        f >> image_name_prefix >> image_name_postfix;
        size_t times_count;
        f >> times_count;
        times.resize(times_count);

        for (size_t i = 0; i < times_count; ++i) {
            f >> times[i];
        }

        f >> roi_xy_min[0] >> roi_xy_min[1] >> roi_xy_max[0] >> roi_xy_max[1];
        f >> downsampling_factor;
        f >> stack_h;
        f >> image_extension;
        f >> ignore_1st_layer;
        f >> backward_calculation;
        f >> is_2D_case;

        if (is_2D_case) {
            ignore_1st_layer = false;
            stack_h = 1;
            z_bounce = 1;
            z_radius = 0;
            image_name_postfix = "";
        }

        f.close();
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