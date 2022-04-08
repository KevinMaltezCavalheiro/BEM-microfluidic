#include <iostream>
#include <vector>   // Why should I include this ?
#include <map>
#include <Eigen/Dense>
#include <list>
#include <Eigen/Sparse> // Eigen library for sparse matrices
#ifdef _MSC_VER
#include <gmsh.h_cwrap>       // gmsh main header
#else
#include <gmsh.h>       // gmsh main header
#endif

using namespace gmsh;
using namespace std;
using namespace Eigen;

Matrix<double, 3, 3> Compute_H_matrix(const double E, const double nu);
void Assembly_K_triplet(const int i, const int el, const vector<vector<size_t>> & nodeTags, const int numNodes, map<int, int> & node_map, list<Triplet <double>> & k_tripletList, Matrix<double, Dynamic, Dynamic> & K_matrix, vector<double> & f_vector);
void Assembly_f_thread(const int i, const int el, const vector<vector<size_t>> & nodeTags, const int numNodes, map<int, int> & node_map, list<pair<int,double>> & thread_f_vector, vector<double> & f_vector);
void Compute_B_matrix(const int j, const int numNodes, const vector<double> & basisFunctionsGradient, const vector<double> & determinants, const Matrix<double, 2, 2> & jacobinvtrans, const Matrix<double, 3, 3> & H_matrix, Matrix<double, Dynamic, Dynamic> & B_matrix);
void Get_physical_Prop_Domain(const vector<string> & keys, double & E, double & nu, double & rho, double & bx, double & by);
void Loop_K_F(const Matrix<double, 3, 3> H_matrix, const double rho, const double bx, const double by, const string & integration_rule, int & dim, vector<int> & tags, map<int, int>  & node_map, vector<list<Triplet <double>>> & my_k_list, vector<list<pair<int,double>>> & thread_f_vector);
void Compute_Neumann_BC(const int numNodes, const int i, const vector<vector<size_t>> & elementTags, const double tx, const double ty, map<int, int> & node_map, vector<double> & f_vector, const vector<double> & weights, const vector<double> & basisFunctions, const vector<double> & determinants, const vector<vector<size_t>> & nodeTags, vector<double> & full_f_vector);
void Neumann_BC(const vector<string> & keys, const string & integration_rule, map<int, int> & node_map, map<string, pair<int, int>> & groups, vector<double> & full_f_vector);
void Fill_New_indices_Dir_BC(const vector<string> & keys, map<int, int> & node_map, map<string, pair<int, int>> & groups, vector<int> & new_indices, vector<double> & dir_BC);
void Fill_K_F_Reduced(const vector<double> & full_f_vector, const SparseMatrix<double> & full_K_matrix, const vector<int> & new_indices, const vector<double> & dir_BC, SparseMatrix<double> & reduced_K_matrix, vector<double> & reduced_f_vector, list<Triplet <double>> & reduced_k_tripletList);
void Reaction_Forces(double & Rx, double & Ry, map<string, pair<int, int>> & groups, map<int, int> & node_map, VectorXd & nodal_forces, const vector<string> keys);
vector<string> Plot_Save_ux_uy(const vector<size_t> & FEM_nodeTags, const vector<double> & full_d_vector, const int nbViews);
void Compute_Stress_Strain_Elemental_Nodal(const vector<int> & elementTypes, const vector<vector<size_t>> &elementTags, const vector<vector<size_t>> &elnodeTags, const vector<double> &full_d_vector, vector<size_t> & elems2D, const Matrix<double, 3, 3> &H_matrix, const double nu, const double E, map<int, int> & node_map, vector<vector<double>> & data3, vector<vector<double>> & data4, vector<vector<double>> & data5, vector<vector<double>> & data6, vector<vector<double>> & data7, vector<vector<double>> & data8, vector<vector<double>> & data9, vector<vector<double>> & data10, int & k);
void applyElecPressure(const string & integration_rule, map<int, int> &  node_map, map<string, pair<int, int>> & groups, map<int, double> &electrostaticPressure, vector<double> & full_f_vector);

int solverFEM(int argc, char **argv, map<int, double> &electrostaticPressure, const int nbViews);
