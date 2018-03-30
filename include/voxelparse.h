#include <string>
#include "CompFab.h"

int generateMtl( std::string filename, uint8_t num_colors);
int generateObj(std::string filename, CompFab::VoxelGrid * voxel_list, uint8_t num_partitions, double scale);
