#include <gmsh.h>       // gmsh main header
#include <iostream>     // for cout
#include <map>
#include <sstream>  // for stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices (QUESTION: OK to use it ?)
#include <Eigen/SparseCholesky> // solving sparse linear systems (QUESTION: OK to use it ?)
#include <cmath>

using namespace std;
using namespace gmsh;
using namespace Eigen;

int code_kevin(int argc, char **argv);

int initialisation(map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys, string &groupname,vector<int> &tags,string &integration_rule, vectorpair &dimTags);

int compute_G_matrix(map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys, string &groupname,vector<int> &tags,string &integration_rule, vectorpair &dimTags);

void sorting_nodes(vector<size_t> &nodeTags,vector<double> &boundary_nodes_coordinates_x,vector<double> &boundary_nodes_coordinates_y,vector<double> &boundary_sort_nodes_tags, vector<double> &boundary_sort_length_elements);



