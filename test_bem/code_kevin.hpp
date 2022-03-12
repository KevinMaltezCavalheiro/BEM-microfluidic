#include <gmsh.h>       // gmsh main header
#include <iostream>     // for cout
#include <map>
#include <sstream>  // for stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices (QUESTION: OK to use it ?)
#include <Eigen/SparseCholesky> // solving sparse linear systems (QUESTION: OK to use it ?)
#include <cmath>
#include <iomanip>

using namespace std;
using namespace gmsh;
using namespace Eigen;

int code_kevin(int argc, char **argv);

void show_mat(MatrixXd &A);

void show_vec(VectorXd &A);



