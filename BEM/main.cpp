#include "functionsBEM.hpp"
#include "functionsFEM.hpp"

#include <gmsh.h>
#include <iostream>
#include <omp.h>

int main(int argc, char **argv)
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

    Eigen::initParallel();

    map<int, double> electrostaticPressure;

    int nbViews = 0;
    double start = omp_get_wtime();

    solverBEM(argc, argv, electrostaticPressure, nbViews);
    double bem_end = omp_get_wtime();
    cout << "time for BEM: " << bem_end - start << "\n";
    /*map<int, double>::iterator it;

    for (it = electrostaticPressure.begin(); it != electrostaticPressure.end(); it++)
    {
        cout << it->first    // string (key)
                << ':'
                << it->second   // string's value 
                << endl;
    }*/

    solverFEM(argc, argv, electrostaticPressure, nbViews);

    double fem_end = omp_get_wtime();
    cout << "time for FEM: " << fem_end - bem_end << "\n";

    fltk::run();
    
    finalize();
    return 0;
}
