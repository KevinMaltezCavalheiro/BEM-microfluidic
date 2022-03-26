#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense> 
#include <cmath>
#include <iomanip>
#include "code_kevin.hpp"


using namespace std;
using namespace gmsh;
using namespace Eigen;

int main(int argc, char **argv)
{
    //cout << setprecision(2) << fixed;

    //code_kevin(argc,argv);

    code_py(argc,argv);

    return 0;
}
