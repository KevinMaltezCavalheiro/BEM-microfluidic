#include "functionsBEM.hpp"

#include <gmsh.h>
#include <iostream>

using namespace std;
using namespace gmsh;
using namespace Eigen;

int mainBEM(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    initialize();

    open(argv[1]);
    model::mesh::generate(2);

    // will map the 1d element tag onto the corresponding electroPressure resulting from BEM code
    map<int, double> electrostaticPressure;

    int nbViews = 0;

    //solverBEM(argc, argv, finalElementTags, electrostaticPressure, nbViews);
    solverBEM(argc, argv, electrostaticPressure, nbViews);

    map<int, double>::iterator it;

    for (it = electrostaticPressure.begin(); it != electrostaticPressure.end(); it++)
    {
        cout << it->first    // string (key)
                << ':'
                << it->second   // string's value
                << endl;
    }

    fltk::run();

    finalize();
    return 0;
}
