#include "functionsBEM.hpp"
#include "functionsFEM.hpp"

#include <gmsh.h>
#include <iostream>
#include <omp.h>

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // If compiled with OpenMP support, gmsh::initialize 
    // also sets the number of threads to "General.NumThreads".
    int nthreads = omp_get_max_threads();
    gmsh::initialize();
    omp_set_num_threads(nthreads);

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    Eigen::initParallel();

    bool postProcessing = 0; // pass it as an argument to the solver, only compute post processing at last iteration

    std::map<int, double> electrostaticPressure;

    std::map<int,std::pair<double,double>> boundaryDisplacementMap;

    int nbViews = 0;
    //double start = omp_get_wtime();

    int iteration;

    int viewTagU = gmsh::view::add("u"); // à modifier plus tard
    int viewTagF = gmsh::view::add("f");

    for(iteration = 1; iteration < 25; iteration ++)
    {
        solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
        solverFEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);
    }

    postProcessing = true;
    //solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
    //solverFEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
    

    //double fem_end = omp_get_wtime();
    //std::cout << "time for FEM: " << fem_end - bem_end << "\n";

    gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}