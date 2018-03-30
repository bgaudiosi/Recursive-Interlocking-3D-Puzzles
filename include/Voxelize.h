#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <iostream>
#include <vector>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle);
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir);
CompFab::VoxelGrid * loadMesh(const char *filename, unsigned int dim);
void saveVoxelsToObj(const char * outfile, CompFab::VoxelGrid * voxel_list);
CompFab::VoxelGrid * objToVoxelGrid(const char * filename, int dim);

#endif
