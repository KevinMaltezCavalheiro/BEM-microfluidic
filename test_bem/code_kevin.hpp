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

int code_simple(int argc, char **argv);

int code_kevin(int argc, char **argv);

void initialise_H_matrix(Matrix<double, 3, 3> &H_matrix, double nu, double E);

void initialisation(double &E,double &nu,double &rho,double &bx,double &by,map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys);

int compute_K_matrix(vector<size_t> &all_nodeTags, map<string, pair<int, int>> &groups, string &integration_rule, double bx, double rho, double by, Matrix<double, 3, 3> &H_matrix, SparseMatrix<double> &full_K_matrix,vector<double> &full_f_vector, string groupname);

int boundary_conditions(vector<string> &keys, string &groupname, double &tx, double &ty, map<string, pair<int, int>> &groups, vector<int> &tags, string &integration_rule, vector<double> &full_f_vector);

int dirichlet(vector<int> &new_indices, vector<double> &dir_BC, vector<string> &keys, double &ux, double &uy, string &groupname, map<string, pair<int, int>> &groups, vector<int> &tags, int &number_removed_lines);

void filling(vector<size_t> &all_nodeTags, int &number_removed_lines, vector<int> &new_indices, SparseMatrix<double> &full_K_matrix, vector<double> &full_f_vector, vector<double> &dir_BC,SparseMatrix<double> &reduced_K_matrix,vector<double> &reduced_f_vector);

void reduce_system_solver(vector<size_t> &all_nodeTags, int &number_removed_lines, vector<double> &reduced_f_vector, SparseMatrix<double> &reduced_K_matrix, VectorXd &final_reduced_d, VectorXd &final_reduced_f);

void nodal_values(vector<double> &full_d_vector, vector<int> &new_indices, vector<double> &dir_BC, VectorXd &final_reduced_d);

void plot_and_save(vector<size_t> &all_nodeTags, vector<double> &full_d_vector, vector<string> &names);

void post_processing(string &integration_rule, vector<double> &full_d_vector, Matrix<double, 3, 3> &H_matrix, double &nu, double &E, vector<string> &names);

