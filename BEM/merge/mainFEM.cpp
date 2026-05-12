#include "functionsFEM.hpp"

#include <gmsh.h>
#include <iostream>
#include <omp.h>

int main(int argc, char **argv)
{

    bool nonLinearSolver = true; // use non-linear solver or not

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

    std::map<int, double> electrostaticPressure;

    if(nonLinearSolver)
    {
        std::map<int,std::pair<double,double>> boundaryDisplacementMap;

        int viewTagU = gmsh::view::add("u"); // à modifier plus tard
        int viewTagF = gmsh::view::add("f");

        solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, 0, true, 1, viewTagU, viewTagF); //iteration randomly set to 1
    }
    else
    {
        solverFEM(argc, argv, electrostaticPressure, 0);
    }

    gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}