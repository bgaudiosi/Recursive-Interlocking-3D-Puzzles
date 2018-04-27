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
        std::cout<<"Usage: puzzle InputMeshFilename OutputMeshFilename Dim NumPieces\n";
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
    int num_pieces = atoi(argv[4]);
    int m = num_voxels/ num_pieces;

    AccessibilityGrid * scores = accessibilityScores(voxel_list, 0.1, 3);
    std::vector<Voxel> seeds = findSeeds(voxel_list);
    
    int seed_choice;
    Voxel seed;
    Voxel normal;
    Voxel prevNormal(0,0,1);
    std::vector<Voxel> anchors;
    std::vector<VoxelPair> interlock;
    std::vector<Voxel> key;
    std::vector<Voxel> testKey;
    int blocker;
    bool okay;
    int expand = m;
    while (key.size() < (int)(0.75*m)) {
        seed_choice = (int)(seeds.size() * ((double) std::rand() / (RAND_MAX)));
        seed = seeds[seed_choice];
        std::cout << "Seed is " << seed.toString() << std::endl;
    
        normal = findNormal( voxel_list, seed, prevNormal);
        anchors = findAnchors(voxel_list, seed, normal, prevNormal);
        interlock = bfs(voxel_list, scores, seed, prevNormal, 50, 10);
    
        key = filterKey(voxel_list, scores, seed, interlock, normal, prevNormal, &blocker);
    
        anchors.push_back(finalAnchor(voxel_list, seed, interlock[blocker], normal));
    
        okay = false;

        while (!okay) {
            testKey = expandPiece( voxel_list, scores, key, anchors, m, prevNormal);
            expand++;
            okay = verifyPiece(voxel_list, testKey);
            if (okay) {
                key = testKey;
            }
        }
    }

    for (int i = 0; i < key.size(); i ++) {
        voxel_list->isInside(key[i].x, key[i].y, key[i].z) = 2;
    }
    
    
    
    std::vector<Voxel> candidates;
    std::vector<Voxel> prevPiece = key;
    std::vector<Voxel> nextPiece;
    Voxel nextNormal;
    std::vector<Voxel> anchorList;
    std::vector<Voxel> testPiece;
    int voxels_left = num_voxels;
    std::vector<Voxel> normal_list;
    normal_list.push_back(prevNormal);
    int index;
    std::cout << "entering the loop" << std::endl;
    int iterations;
    bool recurse = true;
    if (recurse) {
        iterations = num_pieces;
    } else {
        iterations = num_pieces - 2 / 3;
        m *= 3;
    }
    for (int p = 3; p <= iterations; p++) {
        std::cout << "On piece " << std::to_string(p-1) << std::endl;
        scores = accessibilityScores(voxel_list, 0.1, 3);
        candidates = findCandidateSeeds(voxel_list, scores, prevPiece, prevNormal);
        nextPiece.clear();
        while (nextPiece.size() < (int)(0.75*m) ) {
            nextPiece.clear();
            anchorList.clear();
            if (candidates.size() == 0) {
                std::cout << "There are no seeds! Guess it's time to stop" << std::endl;
                break;
            }
            nextPiece = createInitialPiece(voxel_list, scores, prevPiece, candidates, &index);
            if (nextPiece.size() == 0) {
                std::cout << "Failed creating initial piece" << std::endl;
                break;
            }
            nextNormal = findNormalDirection( voxel_list, nextPiece[0], prevPiece);
            nextPiece = ensureInterlocking(voxel_list, scores, prevPiece, nextPiece, p-1, prevNormal, &anchorList);
            /*
            std::cout << "after ensuring interlocking, piece is:" << std::endl;
            printList(nextPiece);
            */
            nextPiece = ensurePieceConnectivity(voxel_list, nextPiece, nextNormal);
            for (int i = 0; i < nextPiece.size(); i++) {
                std::cout << "piece " << std::to_string(p-1) << " " << nextPiece[i].toString() << std::endl;
                voxel_list->isInside(nextPiece[i].x, nextPiece[i].y, nextPiece[i].z) = p;
            }
            okay = false;
            expand = m;
            while (!okay) {
                testPiece = expandPiece( voxel_list, scores, nextPiece, anchorList, m, nextNormal);
                okay =  verifyPiece(voxel_list, testPiece);
                //okay = true;
                expand++;
                if (okay) {
                    nextPiece = testPiece;
                }
            }
            // If failure, erase and remove this candidate
            if (nextPiece.size() < (int)(0.75*m)) {
                candidates.erase (candidates.begin()+index);
                for (int i = 0; i < nextPiece.size(); i++) {
                    std::cout << "piece " << std::to_string(p-1) << " " << nextPiece[i].toString() << std::endl;
                    voxel_list->isInside(nextPiece[i].x, nextPiece[i].y, nextPiece[i].z) = 1;
                }
            }
        }

        for (int i = 0; i < nextPiece.size(); i++) {
            //std::cout << "upon expansion, piece " << std::to_string(p-1) << " is now " << nextPiece[i].toString() << std::endl;
            voxel_list->isInside(nextPiece[i].x, nextPiece[i].y, nextPiece[i].z) = p;
        }
        std::cout << "Piece " << std::to_string(p-1) << " successfully made!" << std::endl;
        normal_list.push_back(nextNormal);
        prevNormal = nextNormal;
        prevPiece = nextPiece;
    }
    
    std::cout << "Solution is: " << std::endl;
    printList(normal_list);
    generateMtl(filename, 10);
    generateObj(filename, voxel_list, 10, 5.0);
    
    return 1;

}
