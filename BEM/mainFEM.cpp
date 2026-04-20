#include "functionsFEM.hpp"

#include <gmsh.h>
#include <iostream>
#include <omp.h>

// this program implements a FEM solver, either linear or non-linear (by taking large displacements into account)
// to choose the linear solver, please add: SetNumber("Non_linear_solver",0); in the .geo file.
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    double start = omp_get_wtime();

    // If compiled with OpenMP support, gmsh::initialize 
    // also sets the number of threads to "General.NumThreads".
    int nthreads = omp_get_max_threads();
    gmsh::initialize();
    omp_set_num_threads(nthreads);

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    Eigen::initParallel();

    // determine if non linear solver must be set
    bool nonLinearSolver = getNonLinearParameter(); 

    // maps the elementTags to a corresponding electrostaticPressure (only useful for coupled solver)
    std::map<int, double> electrostaticPressure;

    int nbViews = 0; // number of views passed to the BEM solver (only useful for coupled solver)

    if(nonLinearSolver)
    {
        std::map<int,std::pair<double,double>> boundaryDisplacementMap; //(only useful for coupled solver)

        int viewTagU = gmsh::view::add("u"); //viewtag of the displacement view
        int viewTagF = gmsh::view::add("f"); //viewtag of the nodal forces view

        solverFEMnonLinear(electrostaticPressure, boundaryDisplacementMap, nbViews, true, 1, viewTagU, viewTagF, false);
    }
    else
    {
        solverFEM(electrostaticPressure, nbViews);
    }

    double total_time = omp_get_wtime() - start;
    std::cout << "-----------------------------------------\n";
    std::cout << "total execution time: " << total_time << " [s]. \n";

    gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}