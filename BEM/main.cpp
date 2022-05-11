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

    bool nonLinearSolver = false; // use non-linear solver or not -> SET AS AN ARGUMENT OF MAIN, much easier
    bool untangleMesh = false; // untangle geometry (deformed geometry, only available for nonLinear)

    // If compiled with OpenMP support, gmsh::initialize 
    // also sets the number of threads to "General.NumThreads".
    double start = omp_get_wtime();
    int nthreads = omp_get_max_threads();
    gmsh::initialize();
    omp_set_num_threads(nthreads);

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    Eigen::initParallel();

    if(nonLinearSolver)
    {

        bool postProcessing = 0; // pass it as an argument to the solver, only compute post processing at last iteration

        std::map<int, double> electrostaticPressure;

        std::map<int,std::pair<double,double>> boundaryDisplacementMap;

        int nbViews = 0;
        //double start = omp_get_wtime();

        int iteration;

        int viewTagU = gmsh::view::add("u"); // à modifier plus tard
        int viewTagF = gmsh::view::add("f");
    


        solverBEM(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1, false);
        solverFEMnonLinear(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1, viewTagU, viewTagF, false);
        
        std::vector<int> nodeTags(boundaryDisplacementMap.size());
        std::vector<double> displacements(boundaryDisplacementMap.size());
        std::vector<double> previousDisplacements(boundaryDisplacementMap.size());
        std::vector<double> relativeDifference(boundaryDisplacementMap.size());

        double relativeDisplacement;
        const double criticalNorm = 0.0000001;
        const int maxNbIteration = 50; // à modifier

        int i = 0;
        for(auto p : boundaryDisplacementMap)
        {
            nodeTags[i] = p.first;
            displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
            i++;
        }

        for(iteration = 2; iteration <= maxNbIteration; iteration++)
        {        
            solverBEM(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, false);
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
            //std::cout << iteration << ", " << relativeDisplacement << ";\n";
            std::cout << "Iteration: " << iteration << "\t-> Averaged relative displacement = " << relativeDisplacement << "\n";
            if(relativeDisplacement < criticalNorm)
                break;
        }


        nbViews = 0;
        solverBEM(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, false);
        postProcessing = true;
        solverFEMnonLinear(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF, untangleMesh);

        if(untangleMesh)
        {
            std::map<int,std::pair<double,double>> boundaryDisplacementMap2;
            solverBEM(electrostaticPressure, boundaryDisplacementMap2, nbViews, postProcessing, iteration, untangleMesh);
        }
        else
        {
            solverBEM(electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, false);
        }
    }
    else
    {

        std::map<int, double> electrostaticPressure;

        std::map<int,std::pair<double,double>> boundaryDisplacementMap;

        int nbViews = 0;

        solverBEM(electrostaticPressure, boundaryDisplacementMap, nbViews, true, 1, false);
        solverFEM(electrostaticPressure, nbViews);

    }

    //double fem_end = omp_get_wtime();
    //std::cout << "time for FEM: " << fem_end - bem_end << "\n";

    double total_time = omp_get_wtime() - start;
    std::cout << "total execution time: " << total_time << " [s]. \n";

    gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}






// #include "functionsBEM.hpp"
// #include "functionsFEM.hpp"

// #include <gmsh.h>
// #include <iostream>
// #include <omp.h>

// int main(int argc, char **argv)
// {
//     if (argc < 2)
//     {
//         std::cout << "Usage: " << argv[0] << " <geo_file>\n";
//         return 0;
//     }

//     bool nonLinearSolver = true; // use non-linear solver or not

//     // If compiled with OpenMP support, gmsh::initialize 
//     // also sets the number of threads to "General.NumThreads".
//     int nthreads = omp_get_max_threads();
//     gmsh::initialize();
//     omp_set_num_threads(nthreads);

//     gmsh::open(argv[1]);
//     gmsh::model::mesh::generate(2);

//     Eigen::initParallel();

//     if(nonLinearSolver)
//     {

//         bool postProcessing = 0; // pass it as an argument to the solver, only compute post processing at last iteration

//         std::map<int, double> electrostaticPressure;
//         std::map<int,std::pair<double,double>> boundaryDisplacementMap;

//         int nbViews = 0;
//         //double start = omp_get_wtime();

//         int iteration;

//         int viewTagU = gmsh::view::add("u"); // à modifier plus tard
//         int viewTagF = gmsh::view::add("f");

//         solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1);
//         solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1, viewTagU, viewTagF);
        
//         std::vector<int> nodeTags(boundaryDisplacementMap.size());
//         std::vector<double> displacements(boundaryDisplacementMap.size());
//         std::vector<double> previousDisplacements(boundaryDisplacementMap.size());
//         std::vector<double> relativeDifference(boundaryDisplacementMap.size());

//         double norm;
//         const double criticalNorm = 0.0001;
//         const int maxNbIteration = 40;

//         int i = 0;
//         for(auto p : boundaryDisplacementMap)
//         {
//             nodeTags[i] = p.first;
//             displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
//             i++;
//         }

//         for(iteration = 2; iteration <= maxNbIteration; iteration++)
//         {        
//             solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
//             solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);

//             i = 0;
//             norm = 0;
//             for(auto p : boundaryDisplacementMap)
//             {
//                 previousDisplacements[i] = displacements[i];
//                 displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
                
//                 if(previousDisplacements[i] != 0.0)
//                     relativeDifference[i] = abs( (displacements[i] - previousDisplacements[i]) / previousDisplacements[i] );
//                 else
//                     relativeDifference[i] = 0;
                
//                 norm += relativeDifference[i]*relativeDifference[i];
//                 i++;
//             }
//             norm = sqrt(norm/i);
//             std::cout << "Norm (iteration '" << iteration << "'): " << norm << "\n";
//             if(norm < criticalNorm)
//                 break;
//         }

//         postProcessing = true;
//         solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
//         solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);
//     }
//     else
//     {

//         std::map<int, double> electrostaticPressure;

//         std::map<int,std::pair<double,double>> boundaryDisplacementMap;

//         int nbViews = 0;

//         solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, true, 1);
//         solverFEM(argc, argv, electrostaticPressure, nbViews);

//     }

//     //double fem_end = omp_get_wtime();
//     //std::cout << "time for FEM: " << fem_end - bem_end << "\n";

//     gmsh::fltk::run();
    
//     gmsh::finalize();
//     return 0;
// }
