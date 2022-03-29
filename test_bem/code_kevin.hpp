#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <stdio.h>

using namespace std;
using namespace gmsh;
using namespace Eigen;

void show_mat(MatrixXd &A);

void show_vec(VectorXd &A);

double compute_angle(double &deltaX1, double &deltaX2, double &deltaY1, double &deltaY2);

int code_py(int argc, char **argv);

// Each element has all its information gathered as in this structure.
struct elementStruct
{
    int tag;            // tag of the element
    string bcType; // "dirichlet" or "neumann" (at the end, should replace this by an int or boolean)
    double bcValue [2]; // values of [u, u_n]
    int nodeTags [2];   // tags of the [first, second] nodes of the element
};

void get_bc_from_onelab(vector<elementStruct> &elementVector, int &dirichletCount, bool &dispHierarchy, vector<double> &allCoord);

void compute_g(double &jMidX, double &jDeltaX, double &jMidY, double &jDeltaY, double &iMidX, double &iMidY, double &jLength, double &g);

void compute_internal_phi(vector<double> &allCoord,vector<size_t> &NodeTags2d, vector<elementStruct> &elementVector, vector<vector<double>> &data_boundary_nodes, vector<vector<double>> &data2d, map<int,double> &nodeMapDomain);

void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d, vector<int> &elementTypes2d, vector<vector<size_t>> &nodeTags2d, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data);


