#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "../include/ExtractPartitions.h"
#include "../include/CompFab.h"

Voxel::Voxel(int x_val, int y_val, int z_val) {
    x = x_val;
    y = y_val;
    z = z_val;

}

//Grid structure for Voxels
AccessibilityStruct::AccessibilityStruct(CompFab::Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ) {
    m_lowerLeft = lowerLeft;
    m_dimX = dimX;
    m_dimY = dimY;
    m_dimZ = dimZ;
    m_size = dimX*dimY*dimZ;

    //Allocate Memory
    m_scoreArray = new double[m_size];

    for(unsigned int i=0; i<m_size; ++i)
    {
        m_scoreArray[i] = 0.0;
    }

}

AccessibilityStruct::~AccessibilityStruct()
{
    delete[] m_scoreArray;
}

std::string Voxel::toString() {
    std::string vox("(" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ")");
    return vox;
}

void printList(std::vector<Voxel> list) {
    for (int i = 0; i < list.size(); i++) {
        std::cout << list[i].toString() << std::endl;
    }
}

std::vector<Voxel> findSeeds( CompFab::VoxelGrid * voxel_list ) {
    std::vector<Voxel> seeds;

    int nx = voxel_list->m_dimX;
    int ny = voxel_list->m_dimY;
    int nz = voxel_list->m_dimZ;
    
    int above = 0;
    int one_side_adjacent;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                if (voxel_list->isInside(i,j,k)) {
                    // First check if any voxels above
                    for (int z = k+1 ; z < nz; z++) {
                        if (voxel_list->isInside(i,j,z)) {
                            above = 1;
                            break;
                        }
                    }
                    if (above) {
                        above = 0;
                        continue;
                    }
                    //Now, check that only two other faces are open and (except bottom)
                    one_side_adjacent = 0;

                    if ( i != 0 ) {
                        if (voxel_list->isInside(i-1,j,k)) {
                           one_side_adjacent++; 
                        } 
                    }
                    if ( i != nx-1 ) {
                        if (voxel_list->isInside(i+1,j,k)) {
                            one_side_adjacent++;
                        }
                    }
                    if (j != 0) {
                        if (voxel_list->isInside(i,j-1,k)) {
                            one_side_adjacent++;  
                        }
                    }
                    if (j != ny-1) {
                        if (voxel_list->isInside(i,j+1,k)) {
                            one_side_adjacent++;
                        }
                    }
                    if (one_side_adjacent == 3) {
                        seeds.push_back(Voxel(i, j, k)); 
                    } 
                }
            }
        }
    }
    return seeds;
}

std::vector<Voxel> getNeighbors(Voxel voxel, CompFab::VoxelGrid * voxel_list) {
    std::vector<Voxel> neighbors;
    int nx = voxel_list->m_dimX;
    int ny = voxel_list->m_dimY;
    int nz = voxel_list->m_dimZ;
    
    if ( voxel.x != 0 ) {
        if (voxel_list->isInside(voxel.x-1,voxel.y,voxel.z)) {
            neighbors.push_back(Voxel(voxel.x-1,voxel.y,voxel.z));
        }
    }   
    if ( voxel.x != nx-1 ) {
        if (voxel_list->isInside(voxel.x+1,voxel.y,voxel.z)) {
            neighbors.push_back(Voxel(voxel.x+1,voxel.y,voxel.z));
        }
    }   
    if (voxel.y != 0) {
        if (voxel_list->isInside(voxel.x,voxel.y-1,voxel.z)) {
            neighbors.push_back(Voxel(voxel.x,voxel.y-1,voxel.z));
        }
    }   
    if (voxel.y != ny-1) {
        if (voxel_list->isInside(voxel.x,voxel.y+1,voxel.z)) {
            neighbors.push_back(Voxel(voxel.x,voxel.y+1,voxel.z));
        }
    }   
    if (voxel.z != 0) {
        if (voxel_list->isInside(voxel.x,voxel.y,voxel.z-1)) {
            neighbors.push_back(Voxel(voxel.x,voxel.y,voxel.z-1));
        }
    }   
    if (voxel.z != nz-1) {
        if (voxel_list->isInside(voxel.x,voxel.y,voxel.z+1)) {
            neighbors.push_back(Voxel(voxel.x,voxel.y,voxel.z+1));
        }
    }
    return neighbors;
}

AccessibilityGrid * accessibilityScores( CompFab::VoxelGrid * voxel_list, double alpha, unsigned int recurse) {
    int nx = voxel_list->m_dimX;
    int ny = voxel_list->m_dimY;
    int nz = voxel_list->m_dimZ;
    CompFab::Vec3 start = CompFab::Vec3(0.0, 0.0, 0.0);
    AccessibilityGrid * scores = new AccessibilityGrid(start, nx, ny, nz);
    std::cout << "recursion level is: " << recurse << std::endl;
    if (recurse == 0) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    scores->score(i, j, k) = getNeighbors(Voxel(i, j, k), voxel_list).size();         
                }
            }
        }
    } else {
        double current_score;
        int x;
        int y;
        int z;
        double multiplier = pow(alpha, recurse);
        std::vector<Voxel> neighbors;
        AccessibilityGrid * old_scores = accessibilityScores( voxel_list, alpha, recurse-1);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    neighbors = getNeighbors(Voxel(i, j, k), voxel_list);
                    current_score = 0;
                    for (int n = 0; n < neighbors.size(); n++) {
                        x = neighbors[n].x;
                        y = neighbors[n].y;
                        z = neighbors[n].z;
                        current_score += old_scores->score(x, y, z);
                    }
                    current_score *= multiplier;
                    current_score += old_scores->score(i,j,k);
                    scores->score(i,j,k) = current_score;
                }
            }
        }
    }
    return scores;
}


