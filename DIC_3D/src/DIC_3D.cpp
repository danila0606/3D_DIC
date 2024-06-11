#include "DIC_3D.h"

namespace DIC3D {

    std::vector<CalcUnit> GenerateCalculationPath(CalculationPathType path_type, size_t s_x, size_t s_y, size_t s_z, SubsetMap& subset_map) {
        std::vector<CalcUnit> path;
        subset_map = SubsetMap(s_x, std::vector<std::vector<SubsetInfo>>(s_y, std::vector<SubsetInfo>(s_z)));

        if (path_type == CalculationPathType::SnakeTopRight) {
            size_t prev_x, prev_y;
            for (size_t y = 0; y < s_y; ++y) {
                for (size_t p = 0; p < s_x; ++p) {
                    size_t x = p;
                    if (y % 2 == 0)
                        x = s_x - p - 1;
                    SubsetInfo info;
                    info.complete = info.use_init = false;
                    if (y || p) {
                        info.use_init = true;
                        info.init_x = prev_x;
                        info.init_y = prev_y;
                    }
                    prev_x = x; prev_y = y;
                    path.emplace_back(CalcUnit{x, y});
                    for (size_t z = 0; z < s_z; ++z)
                        subset_map[x][y][z] = info;
                }
            }
        }
        else if (path_type == CalculationPathType::FromRightToLeft) {
            for (size_t y = 0; y < s_y; ++y) {
                for (int x = s_x - 1; x >= 0; --x) {
                    SubsetInfo info;
                    info.complete = info.use_init = false;
                    if (y || (x != s_x - 1)) {
                        info.use_init = true;
                        if (x == s_x - 1) {
                            info.init_x = s_x - 1;
                            info.init_y = y - 1;
                        }
                        else {
                            info.init_x = x + 1;
                            info.init_y = y;
                        }
                    }
                    path.emplace_back(CalcUnit{x, y});
                    for (size_t z = 0; z < s_z; ++z) {
                        subset_map[x][y][z] = info;
                    }
                }
            }
        }
        else if (path_type == CalculationPathType::FromTopToBottom) {
            for (size_t x = 0; x < s_x; ++x) {
                for (size_t y = 0; y < s_y; ++y) {
                    SubsetInfo info;
                    info.complete = info.use_init = false;
                    if (y || x) {
                        info.use_init = true;
                        if (y == 0) {
                            info.init_x = x - 1;
                            info.init_y = 0;
                        }
                        else {
                            info.init_x = x;
                            info.init_y = y - 1;
                        }
                    }
                    path.emplace_back(CalcUnit{x, y});
                    for (size_t z = 0; z < s_z; ++z)
                        subset_map[x][y][z] = info;
                }
            }
        }

        return path;
    }




}