#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace gmsh;
using namespace Eigen;

// Each element has all its information gathered as in this structure.
struct elementStruct{
    int tag;            // tag of the element
    string bcType; // "dirichlet" or "neumann" (at the end, should replace this by an int or boolean)
    double bcValue [2]; // values of [u, u_n]
    size_t nodeTags [2];   // tags of the [first, second] nodes of the element
    size_t midNodeTags;
    double length;
    double x1, y1, x2, y2;
    double midX, midY;
    double deltaX, deltaY;
};


int solverBEM(int argc, char **argv, std::map<int, double> &electrostaticPressure, int & nbViews);
void dispNodeCoordFun(const vector<size_t> &allNodeTags, const vector<double> &allCoord);
void dispMatricesFun(const MatrixXd &A, const MatrixXd &B, const MatrixXd &c, const MatrixXd &b, const MatrixXd &x);
void dispValuesFun(const vector<elementStruct> &elementVector);
void getBC(map<size_t, int> allNodeIndex, const vector<double> &allCoord, vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy, vector<size_t> &NodeTags2d_domain, vector<double> &NodeCoord2d_domain,            vector<int> &elementTypes2d_domain_BEM, vector<vector<size_t>> &elementTags2d_domain_BEM, vector<vector<size_t>> &elnodeTags2d_domain_BEM, int &entities_tags_for_physical_group_domain_BEM);
void computeGeometricParam(elementStruct &element, map<std::size_t, int> allNodeIndex, const vector<double> &allCoord);
double computeH(const double &x, const double &y, const elementStruct &element);
double computeG(const double &x, const double &y, const elementStruct &element, const string &integrationType);
void computeData(const vector<elementStruct> &elementVector, const vector<int> &elementTypes, const vector<vector<size_t>> &elementTags, vector<size_t> &tags, vector<vector<double>> &data, const string integrationType);
void displayData(const string viewName, const string modelName, const string fileName, const vector<size_t> &tags, const vector<vector<double>> &data);
void compute_internal_phi(vector<double> &allCoord,vector<size_t> &NodeTags2d_domain, vector<elementStruct> &elementVector, map<int,double> &nodeMapDomain, const string integrationType);
void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d_domain_BEM, vector<int> &elementTypes2d_domain_BEM, vector<vector<size_t>> &nodeTags2d_domain_BEM, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data, int &entities_tags_for_physical_group_domain_BEM);
void disp_element_node_data(vector<int> &nodes_commun_tag, vector<int> &elements_commun_tag);
void display_time(chrono::time_point<chrono::system_clock> &start, chrono::time_point<chrono::system_clock> &end, string &function_name, bool &show_time);



