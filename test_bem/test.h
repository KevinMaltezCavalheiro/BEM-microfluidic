#include <gmsh.h>
#include <iostream>
#include <sstream>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Core> 
#include <Eigen/Sparse> 
#include <set>
#include <cmath>
#include <cstdlib>
#include <algorithm>
//#include <QImage>
#include <cstdio>
#include <vector>

typedef Eigen::Triplet<double> T;

using namespace Eigen;
using namespace std;
using namespace gmsh;

void cheeseHole(double x, double y, double z, double r, double lc, vector<int> &shells, vector<int> &volumes);

void x1(int argc, char **argv);

void x2(int argc, char **argv);

void x3(int argc, char **argv);

void x4(int argc, char **argv);

void x5(int argc, char **argv);

void x6(int argc, char **argv);

void x7(int argc, char **argv);

void t1(int argc, char **argv);

void t2(int argc, char **argv);

void t3(int argc, char **argv);

void t4(int argc, char **argv);

void t5(int argc, char **argv);

void t6(int argc, char **argv);

void t7(int argc, char **argv);

void t8(int argc, char **argv);

void t9(int argc, char **argv);

void t10(int argc, char **argv);

void t11(int argc, char **argv);

void t12(int argc, char **argv);

void t13(int argc, char **argv);

void t14(int argc, char **argv);

void t15(int argc, char **argv);

void t16(int argc, char **argv);

void t17(int argc, char **argv);

void t18(int argc, char **argv);

void t19(int argc, char **argv);

void t20(int argc, char **argv);

void t21(int argc, char **argv);

void functions_test(int argc, char **argv);

void test_group(int argc, char **argv);

void sphere_avec_trous(int argc, char **argv);

void structure_avec_torsion(int argc, char **argv);

void code1_loadgeo(int argc, char **argv);

void code2_getgroups(int argc, char **argv);

void code3_buildphysmap(int argc, char **argv);

void code4_hierarchy(int argc, char **argv);

void code5_elprops(int argc, char **argv);

void code6_allnodes(int argc, char **argv);

void code7_parameters(int argc, char **argv);

void code8_jacobians(int argc, char **argv);

void code9_views(int argc, char **argv);

void mycode(int argc, char **argv);

void Tutorial_sparse_example_details(int argc, char **argv);

void Tutorial_sparse_example(int argc, char **argv);

void buildProblem(vector<T>& coefficients, VectorXd& b, int n);

//void saveAsBitmap(const VectorXd& x, int n, const char* filename);

void insertCoefficient(int id, int i, int j, double w, vector<T>& coeffs, VectorXd& b, const VectorXd& boundary);