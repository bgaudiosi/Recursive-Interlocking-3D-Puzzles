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
    
    int num_voxels = 0;
    for (int i = 0; i<dim; i++) {
        for (int j = 0; j<dim; j++) {
            for (int k = 0; k<dim; k++) {
                if (voxel_list->isInside(i,j,k) == 1) {
                    num_voxels++;
                }
            }
        }
    }
    int num_pieces = 10;
    //int m = 25;
    int m = num_voxels/ num_pieces;

    AccessibilityGrid * scores = accessibilityScores(voxel_list, 0.1, 3);
    
    std::vector<Voxel> seeds = findSeeds(voxel_list);
    
    int seed_choice = (int)(4 * ((double) std::rand() / (RAND_MAX)));
    Voxel seed = seeds[seed_choice];
    std::cout << "seed is " << seed.toString() << std::endl;
    
    //highlight seeds
    
    //for (int i=0; i < seeds.size(); i++) {
    //    voxel_list->isInside(seeds[i].x, seeds[i].y, seeds[i].z) = 4;    
    //}
    
    Voxel bad_normal(0, 0, 1);
    Voxel normal = findNormal( voxel_list, seed, bad_normal);
    
    std::cout << "normal is " << normal.toString() << std::endl;   
    
    std::vector<Voxel> anchors = findAnchors(voxel_list, seed, normal, bad_normal);
    
    std::vector<VoxelPair> interlock = bfs(voxel_list, scores, seed, bad_normal, 50, 10);
    
    int blocker;
    std::vector<Voxel> key = filterKey(voxel_list, scores, seed, interlock, normal, bad_normal, &blocker);
    std::cout << blocker << std::endl;
    
    anchors.push_back(finalAnchor(voxel_list, seed, interlock[blocker], normal));
    for (int i = 0; i < key.size(); i ++) {
        std::cout << "key before: " << key[i].toString() << std::endl;
    }
    key = expandPiece( voxel_list, scores, key, anchors, m, Voxel(0, 0, 1));
    for (int i = 0; i < key.size(); i ++) {
        std::cout << "key after: " << key[i].toString() << std::endl;
        voxel_list->isInside(key[i].x, key[i].y, key[i].z) = 2;
    }
    bool okay = verifyPiece(voxel_list, key);
    if (!okay) {
        std::cout << "That's bad" << std::endl;
    } else {
        std::cout << "That's good" << std::endl;
    }
    
    std::vector<Voxel> candidates;
    std::vector<Voxel> prevPiece = key;
    Voxel prevNormal = Voxel(0, 0, 1);
    std::vector<Voxel> nextPiece;
    Voxel nextNormal;
    std::vector<Voxel> anchorList;
    std::vector<Voxel> testPiece;
    int expand;
    
    std::vector<Voxel> normal_list;
    normal_list.push_back(prevNormal);
    for (int p = 3; p <= num_pieces; p++) {
        scores = accessibilityScores(voxel_list, 0.1, 3);
        anchorList.clear();
        candidates = findCandidateSeeds(voxel_list, scores, prevPiece, prevNormal);
        nextPiece = createInitialPiece(voxel_list, scores, prevPiece, candidates);
        nextNormal = findNormalDirection( voxel_list, nextPiece[0], prevPiece);
        nextPiece = ensureInterlocking(voxel_list, scores, prevPiece, nextPiece, p-1, prevNormal, &anchorList);
        nextPiece = ensurePieceConnectivity(voxel_list, nextPiece, nextNormal);
        for (int i = 0; i < nextPiece.size(); i++) {
            std::cout << "piece " << std::to_string(p-1) << " " << nextPiece[i].toString() << std::endl;
            voxel_list->isInside(nextPiece[i].x, nextPiece[i].y, nextPiece[i].z) = p;
        }
        okay = false;
        expand = m;
        while (!okay) {
            testPiece = expandPiece( voxel_list, scores, nextPiece, anchorList, expand, nextNormal);
            okay =  verifyPiece(voxel_list, testPiece);
            expand++;
            if (okay) {
                nextPiece = testPiece;
            }
        }
        for (int i = 0; i < nextPiece.size(); i++) {
            std::cout << "upon expansion, piece " << std::to_string(p-1) << " is now " << nextPiece[i].toString() << std::endl;
            voxel_list->isInside(nextPiece[i].x, nextPiece[i].y, nextPiece[i].z) = p;
        }
        
        normal_list.push_back(nextNormal);
        prevNormal = nextNormal;
        prevPiece = nextPiece;
    }
    
    std::cout << "Solution is: " << std::endl;
    printList(normal_list);
    generateMtl(filename, 10);
    generateObj(filename, voxel_list, 10, 1.0);

    return 1;

}
