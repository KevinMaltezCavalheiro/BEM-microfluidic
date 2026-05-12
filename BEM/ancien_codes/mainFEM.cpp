#include "functionsFEM.hpp"

#include <gmsh.h>
#include <iostream>
#include <omp.h>

using namespace gmsh;
using namespace std;
using namespace Eigen;

int mainFEM(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // If compiled with OpenMP support, initialize
    // also sets the number of threads to "General.NumThreads".
    int nthreads = omp_get_max_threads();
    initialize();
    omp_set_num_threads(nthreads);

    open(argv[1]);
    model::mesh::generate(2);

    initParallel();

    map<int, double> electrostaticPressure;

    solverFEM(argc, argv, electrostaticPressure,0);

    fltk::run();
    
    finalize();
    return 0;
}
