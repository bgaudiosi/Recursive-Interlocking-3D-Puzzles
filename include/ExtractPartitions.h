#include "CompFab.h"
#include <vector>
#include <tuple>
#include <string>


class Voxel {
    public:
        Voxel( int x_val, int y_val, int z_val);
        int x;
        int y;
        int z;
        int value = 0;
        std::string toString();

};

typedef struct AccessibilityStruct {
    //Square voxels only
    AccessibilityStruct(CompFab::Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ);
    ~AccessibilityStruct();

    inline double & score(unsigned int i, unsigned int j, unsigned int k) {
        return m_scoreArray[k*(m_dimX*m_dimY) + j*m_dimY + i];
    }

    double *m_scoreArray;
    unsigned int m_dimX, m_dimY, m_dimZ, m_size;
    CompFab::Vec3 m_lowerLeft;

} AccessibilityGrid;

std::vector<Voxel> findSeeds( CompFab::VoxelGrid * voxel_list );
unsigned int countNeighbors( CompFab::VoxelGrid * voxel_list, Voxel voxel);
AccessibilityGrid * accessibilityScores( CompFab::VoxelGrid * voxel_list, double alpha, unsigned int recurse);
