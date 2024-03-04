#include "ncorr.h"
#include <iostream>

struct DIC_3D_Input {

    DIC_3D_Input(const std::string& filename) {
        std::ifstream f(filename, std::ios::in);
        if (f.is_open()) {
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
        }

        f.close();
    };

    void debug_print(std::ostream& os) {
        os << "subset size: " << subset_size << ", " << "subset offset: " << subset_offset << std::endl;
        os << "z bounce: " << z_bounce << ", " << "z radius: " << z_radius << std::endl;
        os << "images_folder: " << images_folder << std::endl;
        os << "output_folder: " << output_folder << std::endl;
        os << "image_name_prefix: " << image_name_prefix << " ,image_name_postfix: " << image_name_postfix << std::endl;
        os << "times: ";
        for (size_t i = 0; i < times.size(); ++i) {
            os << times[i] << " ";
        }
        os << "\n" << "roi_min: " << roi_xy_min[0] << ", " << roi_xy_min[1] << " ,roi_max: " << roi_xy_max[0] << ", " << roi_xy_max[1] << std::endl;
        os << "downsampling factor: " << downsampling_factor << std::endl;
    }

    size_t subset_size, subset_offset;
    size_t z_bounce, z_radius;

    std::string images_folder;
    std::string output_folder;
    std::string image_name_prefix, image_name_postfix;
    std::vector<size_t> times;
    std::array<size_t, 2> roi_xy_min, roi_xy_max;
    size_t downsampling_factor;
};