#include "functionsBEM.hpp"
#include "functionsFEM.hpp"

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

    bool untangleMesh = true; // untangle geometry (deformed geometry, only available for nonLinear)

    // If compiled with OpenMP support, gmsh::initialize 
    // also sets the number of threads to "General.NumThreads".
    double start = omp_get_wtime();
    int nthreads = omp_get_max_threads();
    gmsh::initialize();
    omp_set_num_threads(nthreads);

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    Eigen::initParallel();

    // determine if non linear solver must be set
    bool nonLinearSolver = getNonLinearParameter();
    nonLinearSolver = 0;

    //bool nonLinearSolver(0);

    if(nonLinearSolver)
    {

        std::cout << "-----------------------------------------\n";
        std::cout << "NON-LINEAR TWO-WAY COUPLED ITERATIVE SOLVER\n";

        bool postProcessing = 0; // pass it as an argument to the solver, only compute post processing at last iteration

        std::map<int, double> electrostaticPressure;

        std::map<int,std::pair<double,double>> boundaryDisplacementMap;

        int nbViews = 0;
        double start = omp_get_wtime();

        int iteration;

        int viewTagU = gmsh::view::add("u"); // passed as an argument to the FEM solver
        int viewTagF = gmsh::view::add("f");

        solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap, postProcessing, 1, false);
        solverFEMnonLinear(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1, viewTagU, viewTagF, false);
        
        std::vector<int> nodeTags(boundaryDisplacementMap.size());
        std::vector<double> displacements(boundaryDisplacementMap.size());
        std::vector<double> previousDisplacements(boundaryDisplacementMap.size());
        std::vector<double> relativeDifference(boundaryDisplacementMap.size());

        double relativeDisplacement;
        const double criticalNorm = 1e-7; 
        const int maxNbIteration = 50;

        // see report section 5.2 "stopping criterion" for the theoretical explanation
        int i = 0;
        for(auto p : boundaryDisplacementMap)
        {
            nodeTags[i] = p.first;
            displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
            i++;
        }

        for(iteration = 2; iteration <= maxNbIteration; iteration++)
        {        
            solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap, postProcessing, iteration, false);
            solverFEMnonLinear(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF, false);

            i = 0;
            relativeDisplacement = 0;
            for(auto p : boundaryDisplacementMap)
            {
                previousDisplacements[i] = displacements[i];
                displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
                
                if(previousDisplacements[i] != 0.0)
                    relativeDifference[i] = abs( (displacements[i] - previousDisplacements[i]) / previousDisplacements[i] );
                else
                    relativeDifference[i] = 0;
                
                relativeDisplacement += relativeDifference[i]*relativeDifference[i];
                i++;
            }
            relativeDisplacement = sqrt(relativeDisplacement/i);

            std::cout << "Iteration: " << iteration << "\t-> Averaged relative displacement = " << relativeDisplacement << "\n";
            if(relativeDisplacement < criticalNorm)
                break;
        }


        nbViews = 0;
        solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap, postProcessing, iteration, false);
        postProcessing = true;
        solverFEMnonLinear(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF, untangleMesh);

        if(untangleMesh)
        {
            std::map<int,std::pair<double,double>> boundaryDisplacementMap2;
            solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap2, postProcessing, iteration, untangleMesh);
        }
        else
        {
            solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap, postProcessing, iteration, false);
        }
    }
    else
    {

        std::cout << "-----------------------------------------\n";
        std::cout << "LINEAR ONE-WAY COUPLED SOLVER\n";

        std::map<int, double> electrostaticPressure;

        std::map<int,std::pair<double,double>> boundaryDisplacementMap;

        int nbViews = 0;

        solverBEM(electrostaticPressure, nbViews, boundaryDisplacementMap, true, 1, false);
        //solverFEM(electrostaticPressure, nbViews);

    }

    //double fem_end = omp_get_wtime();
    //std::cout << "time for FEM: " << fem_end - bem_end << "\n";

    double total_time = omp_get_wtime() - start;
    std::cout << "total execution time: " << total_time << " [s]. \n";

    gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}

////todo: calculer la force DEP sur la sphère, FEM-BEM couplage pour une goutte ? (avec interface = ressort + untanglemesh), BEM 3D, BEMStokes, electrodesCAD
