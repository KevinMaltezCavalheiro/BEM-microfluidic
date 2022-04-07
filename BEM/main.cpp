#include "functions.hpp"

#include <gmsh.h>
#include <iostream>

using namespace std;
using namespace gmsh;
using namespace Eigen;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    gmsh::initialize();

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    // will map the 1d element tag onto the corresponding electroPressure resulting from BEM code
    std::map<int, double> electrostaticPressure;

    int nbViews = 0;

    //solverBEM(argc, argv, finalElementTags, electrostaticPressure, nbViews);
    solverBEM(argc, argv, electrostaticPressure, nbViews);

    std::map<int, double>::iterator it;

    for (it = electrostaticPressure.begin(); it != electrostaticPressure.end(); it++)
    {
        std::cout << it->first    // string (key)
                << ':'
                << it->second   // string's value
                << std::endl;
    }

    gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
