#include "CompFab.h"
#include <vector>
#include <tuple>
#include <string>


class Voxel {
    public:
        Voxel();
        Voxel( int x_val, int y_val, int z_val);
        int x;
        int y;
        int z;
        int value = 0;
        std::string toString();
        inline void operator+=(const Voxel& a) {
            x += a.x;
            y += a.y;
            z += a.z;
        }
};

inline bool operator==(const Voxel& a, const Voxel& b) { return ((a.x == b.x) && (a.y == b.y) && (a.z == b.z)); }
inline Voxel operator+(const Voxel& a, const Voxel& b) { return Voxel(a.x + b.x, a.y + b.y, a.z + b.z); }
inline Voxel operator-(const Voxel& a, const Voxel& b) { return Voxel(a.x - b.x, a.y - b.y, a.z - b.z); }
inline bool operator!=(const Voxel& a, const Voxel& b) { return !(a == b); }

class VoxelPair {
    public:
        VoxelPair( Voxel blocker_vox, double blocker_score, Voxel blockee_vox, double blockee_score);
        Voxel blocker;
        double blocker_accessibility;

        Voxel blockee;
        double blockee_accessibility;
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

class VoxelSort {
    public:
        VoxelSort( Voxel voxel, double score );
        Voxel voxel;
        double score;
};

void printList(std::vector<Voxel> list);
std::vector<Voxel> findSeeds( CompFab::VoxelGrid * voxel_list );
unsigned int countNeighbors( CompFab::VoxelGrid * voxel_list, Voxel voxel);
AccessibilityGrid * accessibilityScores( CompFab::VoxelGrid * voxel_list, double alpha, unsigned int recurse);
Voxel findNormal( CompFab::VoxelGrid * voxel_list, Voxel voxel, Voxel bad_normal);
std::vector<VoxelPair>  bfs(CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, Voxel seed, Voxel normal, int nb_one, int nb_two);
std::vector<Voxel> shortestPath(CompFab::VoxelGrid * voxel_list, Voxel seed, VoxelPair goal, std::vector<Voxel> anchors);
std::vector<Voxel> filterKey(CompFab::VoxelGrid * voxel_list, 
                            AccessibilityGrid * scores, 
                            Voxel seed, 
                            std::vector<VoxelPair> candidates, 
                            Voxel normal_one, 
                            Voxel normal_two, 
                            int * index);
std::vector<Voxel> findAnchors(CompFab::VoxelGrid * voxel_list, Voxel seed, Voxel normal_one, Voxel normal_two);
Voxel finalAnchor(CompFab::VoxelGrid * voxel_list, Voxel seed, VoxelPair blocks, Voxel normal);
std::vector<Voxel> expandPiece( CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, std::vector<Voxel> key, std::vector<Voxel> anchors, int num_voxels, Voxel normal);
bool verifyPiece( CompFab::VoxelGrid * voxel_list, std::vector<Voxel> piece);
Voxel findNormalDirection( CompFab::VoxelGrid * voxel_list, Voxel voxel, std::vector<Voxel> piece);
std::vector<Voxel> findCandidateSeeds(CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, std::vector<Voxel> piece, Voxel perpendicular);
std::vector<Voxel> seedSorter(CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, std::vector<Voxel> seeds, std::vector<Voxel> piece);
std::vector<Voxel> createInitialPiece(CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, std::vector<Voxel> prevPiece, std::vector<Voxel> candidates);
std::vector<Voxel> ensureInterlocking(CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, std::vector<Voxel> prevPiece, std::vector<Voxel> currentPiece, int prevPieceId, Voxel prevNormal, std::vector<Voxel> * theAnchors);
std::vector<Voxel> bfsTwo(CompFab::VoxelGrid * voxel_list, AccessibilityGrid * scores, Voxel seed, Voxel toBlock, Voxel normal, int nb_one, int nb_two, Voxel * anchor, std::vector<Voxel> anchorList);
std::vector<Voxel> ensurePieceConnectivity(CompFab::VoxelGrid * voxel_list, std::vector<Voxel> piece, Voxel normal);
