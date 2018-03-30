#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <iomanip> // setprecision
#include "../include/CompFab.h"
#include "../include/Mesh.h"
#include "../include/Voxelize.h"
#include "../include/voxelparse.h"
#include "../include/ExtractPartitions.h"

int main(int argc, char **argv)
{
    //fix later
    if(argc < 4)
    {
        std::cout<<"Usage: puzzle InputMeshFilename OutputMeshFilename Dim\n";
        exit(0);
    }
    
    int dim = atoi(argv[3]); //dimension of voxel grid (e.g. 32x32x32)
    std::string filename(argv[2]);

    //CompFab::VoxelGrid * voxel_list = objToVoxelGrid(argv[1], dim);
    CompFab::Vec3 start = CompFab::Vec3(0.0, 0.0, 0.0);
    CompFab::VoxelGrid * voxel_list = new CompFab::VoxelGrid(start, dim, dim, dim, 1.0);
    for (int i = 0; i<dim; i++) {
        for (int j = 0; j<dim; j++) {
            for (int k = 0; k<dim; k++) {
                voxel_list->isInside(i,j,k) = 1;
            }
        }
    }
    std::srand(time(0));
    
    AccessibilityGrid * scores = accessibilityScores(voxel_list, 0.1, 3);
    for (int i = 0; i<dim; i++) {
        for (int j = 0; j<dim; j++) {
            for (int k = 0; k<dim; k++) {
                std::cout << "(" << i << "," << j << "," << k << ") = " << scores->score(i,j,k) << std::endl;
            }
        }
    }
    
    std::vector<Voxel> seeds = findSeeds(voxel_list);
    int seed_choice = (int)(4 * ((double) std::rand() / (RAND_MAX)));
    Voxel seed = seeds[seed_choice];
    std::cout << "seed is " << seed.toString() << std::endl;
    int x;
    int y;
    int z;
    for (int i=0; i < seeds.size(); i++) {
        x = seeds[i].x;
        y = seeds[i].y;
        z = seeds[i].z;
        voxel_list->isInside(x,y,z) = 4;    
    }
    
    CompFab::VoxelGrid * puzzle = voxel_list;
    generateMtl(filename, 4);
    generateObj(filename, puzzle, 4, 1.0);

    return 1;

}
