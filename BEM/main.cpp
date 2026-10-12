#include "functionsBEM.hpp"

#include <gmsh.h>
#include <iostream>
#include <omp.h>

//pour changer le solver linéaire/non linéaire -> nonLinearSolver n'est plus mis à zéro et pour utiliser solver BEM avec ou sans FEM aller flag ou deflag solverFEM

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // If compiled with OpenMP support, gmsh::initialize 
    // also sets the number of threads to "General.NumThreads".
    double start = omp_get_wtime();
    int nthreads = omp_get_max_threads();
    std::cout << "number of threads: " << nthreads << std::endl;
    gmsh::initialize();
    omp_set_num_threads(nthreads);

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    Eigen::initParallel();

    std::cout << "-----------------------------------------\n";
    std::cout << "LINEAR ONE-WAY COUPLED SOLVER\n";

    std::map<int, double> electrostaticPressure;

    std::map<int,std::pair<double,double>> boundaryDisplacementMap;

    int nbViews = 0;

    solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap, true, 1, false);

    double total_time = omp_get_wtime() - start;
    std::cout << "total execution time: " << total_time << " [s]. \n";

    gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}

////todo: calculer la force DEP sur la sphère, FEM-BEM couplage pour une goutte ? (avec interface = ressort + untanglemesh), BEM 3D, BEMStokes, electrodesCAD, permittivité complexe
