#include <gmsh.h>       // gmsh main header
#include <iostream>     // for cout
#include <map>
#include <sstream>  // for stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices (QUESTION: OK to use it ?)
#include <Eigen/SparseCholesky> // solving sparse linear systems (QUESTION: OK to use it ?)
#include <cmath>
//#include "code_louis_function.hpp"
//#include "code_louis.hpp"
#include "code_kevin.hpp"
#include "code_louis.hpp"


using namespace std;
using namespace gmsh;
using namespace Eigen;

int main(int argc, char **argv)
{
    //code_louis(argc,argv);

    code_kevin(argc,argv);

    return 0;
}
