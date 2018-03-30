//Computational Fabrication Assignment #1
// By David Levin 2014
#include <iostream>
#include <vector>
#include <string>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

typedef std::vector<CompFab::Triangle> TriangleList;

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    CompFab::Vec3 e1,e2,h,s,q;
    double a,f,u,v;
    
    e1 = triangle.m_v2 - triangle.m_v1;
    e2 = triangle.m_v3 - triangle.m_v1;
    h = ray.m_direction % e2;
    a = e1*h;
    if (a > -EPSILON && a < EPSILON) {
        return 0;
    }
    
    f = 1.0/a;
    s = ray.m_origin - triangle.m_v1;
    u = f* (s*h);
        
    if (u < 0.0 || u > 1.0) {
        return 0;
    }
    q = s % e1;
    v = f* (ray.m_direction * q);
        
    if (v < 0.0 || u + v > 1.0) {
        return 0;
    }
    double t = f * (e2*q);
    if (t > EPSILON) {
        return 1;
    }

    return 0;

}

//Triangle list (global)
TriangleList g_triangleList;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    CompFab::Ray ray = CompFab::RayStruct(voxelPos, dir);
    for (unsigned int i = 0; i<g_triangleList.size(); i++) {
        if (rayTriangleIntersection(ray, g_triangleList[i])) {
            numHits++;
        }
    }
    return numHits;
}

CompFab::VoxelGrid * loadMesh(const char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    CompFab::VoxelGrid * voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return voxelGrid;
   
}

void saveVoxelsToObj(const char * outfile, CompFab::VoxelGrid * voxelGrid)
{
 
    Mesh box;
    Mesh mout;
    int nx = voxelGrid->m_dimX;
    int ny = voxelGrid->m_dimY;
    int nz = voxelGrid->m_dimZ;
    double spacing = voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

CompFab::VoxelGrid * objToVoxelGrid( const char * filename, int dim) {
    CompFab::VoxelGrid *voxelGrid = loadMesh(filename, dim);
    
    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);

    int intersections;
    for (int i = 0; i<dim; i++) {
        for (int j = 0; j<dim; j++) {
            for (int k = 0; k<dim; k++) {
                voxelPos = CompFab::Vec3Struct(voxelGrid->m_lowerLeft[0] + voxelGrid->m_spacing*i,
                                               voxelGrid->m_lowerLeft[1] + voxelGrid->m_spacing*j,
                                               voxelGrid->m_lowerLeft[2] + voxelGrid->m_spacing*k);

                intersections = numSurfaceIntersections(voxelPos, direction);
                if (intersections % 2 == 1) {
                    voxelGrid->isInside(i,j,k) = 1;
                }
            }
        }
    }

    const char * outfile = "testwrite.obj";
    //Write out voxel data as obj
    saveVoxelsToObj(outfile, voxelGrid);

    return voxelGrid;
}
