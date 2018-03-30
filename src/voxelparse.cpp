#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include "../include/voxelparse.h"

int generateMtl( std::string filename, uint8_t num_colors) {
    std::srand(time(0));
    std::ofstream out_mtl;
    out_mtl.open(filename + ".mtl");
    if (!out_mtl.good()) {
        std::cout << "cannot open mtl" << filename << std::endl;
        return 0 ;
    }

    std::string color_count(std::to_string(num_colors));

    double red_num;
    double green_num;
    double blue_num;

    std::stringstream red;
    std::stringstream green;
    std::stringstream blue;

    std::string mtl("# Blender MTL File: 'None'\n# Material Count: " + color_count + "\n");
    for (int i = 0; i < num_colors; i++) {
        red.str("");
        red.clear();
        green.str("");
        green.clear();
        blue.str("");
        blue.clear();

        red_num = ((double) std::rand() / (RAND_MAX));
        green_num =((double) std::rand() / (RAND_MAX));
        blue_num = ((double) std::rand() / (RAND_MAX));

        red   << std::fixed << std::setprecision(6) << red_num;
        green << std::fixed << std::setprecision(6) << green_num;
        blue  << std::fixed << std::setprecision(6) << blue_num;

        mtl.append("newmtl color." + std::to_string(i) + "\n" +
                   "Ns 0.000000\n" +
                   "Ka 1.000000 1.000000 1.000000\n" +
                   "Kd " + red.str() + " " + green.str() + " " + blue.str() + "\n" +
                   "Ks 0.000000 0.000000 0.000000\n" +
                   "illum 2\n\n"
                   );
    }
    out_mtl << mtl << "#end";
    out_mtl.close();
    return 1;
}

int generateObj(std::string filename, CompFab::VoxelGrid * voxel_list, uint8_t num_partitions, double scale) {
    std::ofstream out(filename + ".obj");
    if(!out.good()){
        std::cout<<"cannot open output file"<<filename<< std::endl;
        return 0;
    }

    std::string obj;
    obj.append("# File generated for my project\n");
    obj.append("mtllib " + filename  + ".mtl\n");
    std::vector<std::string> partitions(num_partitions);

    for (int i = 0; i<partitions.size(); i++) {
        partitions[i].append("o part." + std::to_string(i) + "\n");
    }
    std::vector<uint32_t> partition_sizes(num_partitions, 1);
    std::vector<std::string> partition_vectors(num_partitions);
    std::vector<std::string> partition_faces(num_partitions);

    std::stringstream first_i;
    std::stringstream next_i;
    std::stringstream first_j;
    std::stringstream next_j;
    std::stringstream first_k;
    std::stringstream next_k;

    std::string normals("vn  1.0  0.0  0.0\nvn -1.0  0.0  0.0\nvn  0.0  1.0  0.0\nvn  0.0 -1.0  0.0\nvn  0.0  0.0  1.0\nvn  0.0  0.0 -1.0\n");
    std::string faces;
    std::string colors;

    int nx = voxel_list->m_dimX;
    int ny = voxel_list->m_dimY;
    int nz = voxel_list->m_dimZ;

    uint32_t count = 1;
    unsigned int p;
    for (int current_part = 1; current_part <= num_partitions; current_part++) {
        for (int i = 0; i < nx; i++) {
            for (int j=0; j < ny;j++) {
                for (int k=0; k < nz; k++) {
                    p = voxel_list->isInside(i,j,k);
                    if (p == current_part) {
                        p -= 1;
                        // Add vectors
                        first_i.str("");
                        first_i.clear();
                        first_j.str("");
                        first_j.clear();
                        first_k.str("");
                        first_k.clear();
                        next_i.str("");
                        next_i.clear();
                        next_j.str("");
                        next_j.clear();
                        next_k.str("");
                        next_k.clear();

                        first_i << std::fixed << std::setprecision(6) << (double)(i * scale);
                        next_i << std::fixed << std::setprecision(6) << (double)((i + 1) * scale);
                        first_j << std::fixed << std::setprecision(6) << (double)(j * scale);
                        next_j << std::fixed << std::setprecision(6) << (double)((j + 1) * scale);
                        first_k << std::fixed << std::setprecision(6) << (double)(k * scale);
                        next_k << std::fixed << std::setprecision(6) << (double)((k + 1) * scale);

                        partition_vectors[p].append("v " + first_i.str() + " " + first_j.str() + " " + first_k.str() + "\n");
                        partition_vectors[p].append("v " + first_i.str() + " " + next_j.str() + " " + first_k.str() + "\n");
                        partition_vectors[p].append("v " + first_i.str() + " " + first_j.str() + " " + next_k.str() + "\n");
                        partition_vectors[p].append("v " + first_i.str() + " " + next_j.str() + " " + next_k.str() + "\n");
                        partition_vectors[p].append("v " + next_i.str() + " " + first_j.str() + " " + first_k.str() + "\n");
                        partition_vectors[p].append("v " + next_i.str() + " " + next_j.str() + " " + first_k.str() + "\n");
                        partition_vectors[p].append("v " + next_i.str() + " " + first_j.str() + " " + next_k.str() + "\n");
                        partition_vectors[p].append("v " + next_i.str() + " " + next_j.str() + " " + next_k.str() + "\n");

                        // Add faces
                        partition_faces[p].append("f " + std::to_string(count + 4) + "//" + "1 " +
                                                         std::to_string(count + 6) + "//" + "1 " +
                                                         std::to_string(count + 7) + "//" + "1" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 4) + "//" + "1 " +
                                                         std::to_string(count + 7) + "//" + "1 " +
                                                         std::to_string(count + 5) + "//" + "1" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 0) + "//" + "2 " +
                                                         std::to_string(count + 3) + "//" + "2 " +
                                                         std::to_string(count + 2) + "//" + "2" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 0) + "//" + "2 " +
                                                         std::to_string(count + 1) + "//" + "2 " +
                                                         std::to_string(count + 3) + "//" + "2" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 2) + "//" + "3 " +
                                                         std::to_string(count + 7) + "//" + "3 " +
                                                         std::to_string(count + 6) + "//" + "3" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 2) + "//" + "3 " +
                                                         std::to_string(count + 3) + "//" + "3 " +
                                                         std::to_string(count + 7) + "//" + "3" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 0) + "//" + "4 " +
                                                         std::to_string(count + 4) + "//" + "4 " +
                                                         std::to_string(count + 5) + "//" + "4" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 0) + "//" + "4 " +
                                                         std::to_string(count + 5) + "//" + "4 " +
                                                         std::to_string(count + 1) + "//" + "4" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 1) + "//" + "5 " +
                                                         std::to_string(count + 5) + "//" + "5 " +
                                                         std::to_string(count + 7) + "//" + "5" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 1) + "//" + "5 " +
                                                         std::to_string(count + 7) + "//" + "5 " +
                                                         std::to_string(count + 3) + "//" + "5" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 0) + "//" + "6 " +
                                                         std::to_string(count + 6) + "//" + "6 " +
                                                         std::to_string(count + 4) + "//" + "6" + "\n");
                        partition_faces[p].append("f " + std::to_string(count + 0) + "//" + "6 "  +
                                                         std::to_string(count + 2) + "//" + "6 "  +
                                                         std::to_string(count + 6) + "//" + "6" + "\n");
                        count += 8;
                    }
                }
            }
        }
    }

    for (int i = 0; i<partitions.size(); i++) {
        partitions[i].append(partition_vectors[i]);
        partitions[i].append(normals);
        partitions[i].append("usemtl color." + std::to_string(i) + "\n");
        partitions[i].append("s 1\n");
        partitions[i].append(partition_faces[i]);
        obj.append(partitions[i]);
    }
    out << obj << "#end";
    out.close();
    return 1;
}

