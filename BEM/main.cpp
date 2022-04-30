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

    bool nonLinearSolver = false; // use non-linear solver or not

    // If compiled with OpenMP support, gmsh::initialize
    // also sets the number of threads to "General.NumThreads".
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

        // for(iteration = 1; iteration < 23; iteration ++) // implement break condition
        // {
        //     solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
        //     solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);
        // }



        solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1);
        solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1, viewTagU, viewTagF);

        std::vector<int> nodeTags(boundaryDisplacementMap.size());
        std::vector<double> displacements(boundaryDisplacementMap.size());
        std::vector<double> previousDisplacements(boundaryDisplacementMap.size());
        std::vector<double> relativeDifference(boundaryDisplacementMap.size());

        double relativeDisplacement;
        const double criticalNorm = 0.0000001;
        const int maxNbIteration = 100;

        int i = 0;
        for(auto p : boundaryDisplacementMap)
        {
            nodeTags[i] = p.first;
            displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
            i++;
        }

        for(iteration = 2; iteration <= maxNbIteration; iteration++)
        {
            solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
            solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);

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
            // std::cout << "Iteration: " << iteration << "\t-> Averaged relative displacement = " << relativeDisplacement << "\n";
            if(relativeDisplacement < criticalNorm)
                break;
        }



        postProcessing = true;
        solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
        solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);

        /*--------get physical groups--------------*/
        std::map<std::string, std::pair<int, int>> groups;
        gmsh::vectorpair dimTags;
        gmsh::model::getPhysicalGroups(dimTags);
        for (std::size_t i = 0; i < dimTags.size(); ++i)
        {
            int dim = dimTags[i].first;
            int tag = dimTags[i].second;
            std::string name;
            gmsh::model::getPhysicalName(dim, tag, name);
            groups[name] = {dim, tag};
        }

        std::string groupname = "BEM_domain";
        int dim = groups[groupname].first;
        int tag = groups[groupname].second;

        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

        gmsh::vectorpair dimTagsBis(tags.size());
        for(size_t i = 0; i < tags.size(); i++)
        {
            dimTagsBis[i] = std::pair<int,int>(dim, tags[i]);
        }

        for(auto &p:boundaryDisplacementMap)
        {
            std::vector<double> coord, parametricCoord;
            int dim, tag;
            gmsh::model::mesh::getNode(p.first, coord, parametricCoord, dim, tag);
            std::vector<double> total_disp = {coord[0] + p.second.first, coord[1] + p.second.second, coord[2]};
            gmsh::model::mesh::setNode(p.first, total_disp, {});
        }
        gmsh::model::mesh::optimize("UntangleMeshGeometry", false, 1, dimTagsBis);
        std::map<int,std::pair<double,double>> boundaryDisplacementMap2;
        solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap2, nbViews, postProcessing, iteration);
    }
    else
    {
        std::map<int, double> electrostaticPressure;

        std::map<int,std::pair<double,double>> boundaryDisplacementMap;

        int nbViews = 0;

        solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, true, 1);
        solverFEM(argc, argv, electrostaticPressure, nbViews);

    }

    //double fem_end = omp_get_wtime();
    //std::cout << "time for FEM: " << fem_end - bem_end << "\n";

    gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}






// #include "functionsBEM.hpp"
// #include "functionsFEM.hpp"
//
// #include <gmsh.h>
// #include <iostream>
// #include <omp.h>
//
// int main(int argc, char **argv)
// {
//     if (argc < 2)
//     {
//         std::cout << "Usage: " << argv[0] << " <geo_file>\n";
//         return 0;
//     }
//
//     bool nonLinearSolver = true; // use non-linear solver or not
//
//     // If compiled with OpenMP support, gmsh::initialize
//     // also sets the number of threads to "General.NumThreads".
//     int nthreads = omp_get_max_threads();
//     gmsh::initialize();
//     omp_set_num_threads(nthreads);
//
//     gmsh::open(argv[1]);
//     gmsh::model::mesh::generate(2);
//
//     Eigen::initParallel();
//
//     if(nonLinearSolver)
//     {
//
//         bool postProcessing = 0; // pass it as an argument to the solver, only compute post processing at last iteration
//
//         std::map<int, double> electrostaticPressure;
//         std::map<int,std::pair<double,double>> boundaryDisplacementMap;
//
//         int nbViews = 0;
//         //double start = omp_get_wtime();
//
//         int iteration;
//
//         int viewTagU = gmsh::view::add("u"); // à modifier plus tard
//         int viewTagF = gmsh::view::add("f");
//
//         solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1);
//         solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, 1, viewTagU, viewTagF);
//
//         std::vector<int> nodeTags(boundaryDisplacementMap.size());
//         std::vector<double> displacements(boundaryDisplacementMap.size());
//         std::vector<double> previousDisplacements(boundaryDisplacementMap.size());
//         std::vector<double> relativeDifference(boundaryDisplacementMap.size());
//
//         double norm;
//         const double criticalNorm = 0.0001;
//         const int maxNbIteration = 40;
//
//         int i = 0;
//         for(auto p : boundaryDisplacementMap)
//         {
//             nodeTags[i] = p.first;
//             displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
//             i++;
//         }
//
//         for(iteration = 2; iteration <= maxNbIteration; iteration++)
//         {
//             solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
//             solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);
//
//             i = 0;
//             norm = 0;
//             for(auto p : boundaryDisplacementMap)
//             {
//                 previousDisplacements[i] = displacements[i];
//                 displacements[i] = sqrt(p.second.first*p.second.first + p.second.second*p.second.second);
//
//                 if(previousDisplacements[i] != 0.0)
//                     relativeDifference[i] = abs( (displacements[i] - previousDisplacements[i]) / previousDisplacements[i] );
//                 else
//                     relativeDifference[i] = 0;
//
//                 norm += relativeDifference[i]*relativeDifference[i];
//                 i++;
//             }
//             norm = sqrt(norm/i);
//             std::cout << "Norm (iteration '" << iteration << "'): " << norm << "\n";
//             if(norm < criticalNorm)
//                 break;
//         }
//
//         postProcessing = true;
//         solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration);
//         solverFEMnonLinear(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, postProcessing, iteration, viewTagU, viewTagF);
//     }
//     else
//     {
//
//         std::map<int, double> electrostaticPressure;
//
//         std::map<int,std::pair<double,double>> boundaryDisplacementMap;
//
//         int nbViews = 0;
//
//         solverBEM(argc, argv, electrostaticPressure, boundaryDisplacementMap, nbViews, true, 1);
//         solverFEM(argc, argv, electrostaticPressure, nbViews);
//
//     }
//
//     //double fem_end = omp_get_wtime();
//     //std::cout << "time for FEM: " << fem_end - bem_end << "\n";
//
//     gmsh::fltk::run();
//
//     gmsh::finalize();
//     return 0;
// }
