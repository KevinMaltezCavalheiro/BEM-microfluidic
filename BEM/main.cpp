#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense> 
#include <cmath>
#include <iomanip>
# include <stdio.h>
#include <chrono>
#include "functions.hpp"
#include <omp.h>

using namespace std;
using namespace gmsh;
using namespace Eigen;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    initialize();

    open(argv[1]);
    model::mesh::generate(2);

//    #pragma omp parallel
//    #pragma omp critical
//    cout << "Greetings from thread "<< omp_get_thread_num() << endl;

    vector<int> elementTypes;
    vector<std::vector<size_t>> elementTags;
    vector<std::vector<size_t>> nodeTags;
    model::mesh::getElements(elementTypes, elementTags, nodeTags, 1);

    int nb1dElements = 0;
    for(int i = 0; i < elementTypes.size(); i++)
    {
        nb1dElements += elementTags[i].size();
    }

    vector<size_t> finalElementTags(nb1dElements);
    vector<double> electrostaticPressure(nb1dElements);

    solverBEMFun(argc, argv, finalElementTags, electrostaticPressure);

//    for(int i = 0; i < finalElementTags.size(); i++)
//    {
//        cout << "Element tag '" << finalElementTags[i] << "'\t-> Electrostatic pressure = " << electrostaticPressure[i] << endl;
//    }

    finalize();

    return 0;
}
