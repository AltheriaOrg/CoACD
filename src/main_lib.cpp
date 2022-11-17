#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <thread>

#include "main_lib.h"
#include "process.h"
#include "DebugCPP.h"

using namespace std;

struct ConvexHull
{
    double* points;
    unsigned int* triangles;
    unsigned int nPoints;
    unsigned int nTriangles;
};

// load mesh in dynamic memory
EXTERN void* LoadMesh(const float* const p_points,
    const uint32_t countPoints,
    const uint32_t* const p_triangles,
    const uint32_t countTriangles)
{
    Model* m = new Model();
    m->Load(p_points,countPoints,p_triangles,countTriangles);
    return m;
}

EXTERN void DestroyMesh(Model* m)
{
    delete m;
}

EXTERN void DestroyConvexHulls(std::vector<struct ConvexHull>* cvxs)
{
    delete cvxs;
}

EXTERN uint32_t GetNConvexHulls(std::vector<struct ConvexHull>* cvxs) {
    return cvxs->size();
}

EXTERN void GetConvexHull(
    std::vector<struct ConvexHull>* convex_hulls,
    const uint32_t index,
    struct ConvexHull* chM) {

    *chM = (*convex_hulls)[index];
}

// compute with no output
EXTERN void* GenerateConvexDecomposition(Model* mesh, double threshold)
{
    const auto processor_count = std::thread::hardware_concurrency();
    omp_set_num_threads(processor_count - 1);

    Params params;
    params.seed = (unsigned)time(NULL);
    params.threshold = min(max(threshold, 0.01), 1.0);
    
    ofstream of;
    Model temp;
    vector<Model>* convex_hulls;
    vector<struct ConvexHull>* convex_hulls_structs = new vector<struct ConvexHull>();
    ManifoldPreprocess(params, *mesh, temp);

    temp.PCA();
    temp.Normalize();

    convex_hulls = ComputeLib(of, temp, params);

    of.close();

    // copy values to convex hull structs in order to
    // have similar structure to VHACD unity library
    for (int i = 0; i < convex_hulls->size(); ++i) {
        Model* current_ch = &((*convex_hulls)[i]);
        
        // revert model transforms
        current_ch->Recover();
        current_ch->RevertPCA();

        struct ConvexHull new_ch;
        new_ch.nPoints = current_ch->points.size();
        new_ch.nTriangles = current_ch->triangles.size();
        new_ch.points = (double*)malloc(new_ch.nPoints * sizeof(double) * 3);
        for (int j = 0; j < new_ch.nPoints; ++j) {
            new_ch.points[j * 3] = current_ch->points[j][0];
            new_ch.points[j * 3+1] = current_ch->points[j][1];
            new_ch.points[j * 3+2] = current_ch->points[j][2];
        }

        new_ch.triangles = (unsigned int*)malloc(new_ch.nTriangles * sizeof(unsigned int) * 3);
        for (int j = 0; j < new_ch.nTriangles; ++j) {
            new_ch.triangles[j * 3] = (unsigned int)current_ch->triangles[j][0];
            new_ch.triangles[j * 3 + 1] = (unsigned int)current_ch->triangles[j][1];
            new_ch.triangles[j * 3 + 2] = (unsigned int)current_ch->triangles[j][2];
        }
        convex_hulls_structs->push_back(new_ch);
    }

    // return convex hulls -> vector<struct ConvexHull>
    return convex_hulls_structs;
}
