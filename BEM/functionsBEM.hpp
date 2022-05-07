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

// Each element has all its information gathered as in this structure.
struct elementStruct{
    int tag;            // tag of the element
    std::string bcType; // "dirichlet" or "neumann" (at the end, should replace this by an int or boolean)
    double bcValue [2]; // values of [u, u_n]
    size_t nodeTags [2];   // tags of the [first, second] nodes of the element
    size_t midNodeTags;
    double length;
    double x1, y1, x2, y2;
    double midX, midY;
    double deltaX, deltaY;
};


void solverBEM(std::map<int, double> &electrostaticPressure, std::map<int,std::pair<double,double>> & boundaryDisplacementMap, int & nbViews, bool postProcessing, const int iteration);
void singleDomainSolverBEM(std::map<int, double> &electrostaticPressure, std::map<int,std::pair<double,double>> & boundaryDisplacementMap, int & nbViews, bool postProcessing, const int iteration, const int domainTag);
void dispNodeCoordFun(const std::vector<size_t> &allNodeTags, const std::vector<double> &allCoord);
void dispMatricesFun(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const Eigen::MatrixXd &c, const Eigen::MatrixXd &b, const Eigen::MatrixXd &x);
void dispValuesFun(const std::vector<elementStruct> &elementVector);
void getBC(std::map<size_t, int> &allNodeIndex, const std::vector<double> &allCoord, std::vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy, std::vector<size_t> &NodeTags2d_domain, std::vector<double> &NodeCoord2d_domain, std::vector<int> &elementTypes2d_domain_BEM, std::vector<std::vector<size_t>> &elementTags2d_domain_BEM, std::vector<std::vector<size_t>> &elnodeTags2d_domain_BEM, int &entities_tags_for_physical_group_domain_BEM, std::vector<int> &entities_tags_2d_domain_BEM, std::map<int, std::pair<double, double>> &boundaryDisplacementMap, std::map<std::size_t, std::pair<bool,bool>> &nodeOnBoundary, const int domainTag);
void computeGeometricParam(elementStruct &element, std::map<std::size_t, int> &allNodeIndex, const std::vector<double> &allCoord, std::pair<double, double> &nodalDisp1, std::pair<double, double> &nodalDisp2);
double computeH(const double &x, const double &y, const elementStruct &element);
double computeG(const double &x, const double &y, const elementStruct &element, const std::string &integrationType);
void computeData(const std::vector<elementStruct> &elementVector, std::vector<size_t> &tags, std::vector<std::vector<double>> &data, const std::string integrationType, std::vector<int> &entities_tags_2d_domain_BEM, bool &solution_type);
void displayData(const std::string viewName, const std::string modelName, const std::string fileName, const std::vector<size_t> &tags, const std::vector<std::vector<double>> &data);
void compute_internal_phi(std::vector<double> &allCoord, std::vector<size_t> &NodeTags2d_domain, std::map<size_t, int> &allNodeIndex, std::vector<elementStruct> &elementVector, std::map<int,double> &nodeMapDomain, const std::string integrationType, std::map<std::size_t, std::pair<bool,bool>> &nodeOnBoundary, bool &solution_type);
void element_node_data(std::vector<std::vector<std::vector<double>>> &data, std::vector<std::vector<size_t>> &elementTags2d_domain_BEM, std::vector<int> &elementTypes2d_domain_BEM, std::vector<std::vector<size_t>> &nodeTags2d_domain_BEM, std::vector<elementStruct> &elementVector, std::map<int, double> &nodeMapDomain, bool &show_element_node_data, int &entities_tags_for_physical_group_domain_BEM, std::map<std::size_t, std::pair<bool,bool>> &nodeOnBoundary);
void disp_element_node_data(std::vector<int> &nodes_commun_tag, std::vector<int> &elements_commun_tag);
void display_time(std::chrono::time_point<std::chrono::system_clock> &start, std::chrono::time_point<std::chrono::system_clock> &end, std::string &function_name, bool &show_time);
void convergenceFun(std::vector<int> &entities_tags_2d_domain_BEM, std::vector<size_t> &tags, std::vector<std::vector<double>> &phi);
double computeGradHx(const double &x, const double &y, const elementStruct &element, const std::string &integrationType);
double computeGradHy(const double &x, const double &y, const elementStruct &element, const std::string &integrationType);
double computeGradGx(const double &x, const double &y, const elementStruct &element, const std::string &integrationType);
double computeGradGy(const double &x, const double &y, const elementStruct &element, const std::string &integrationType);
double computeGradHx_analytique(const double &x, const double &y, const elementStruct &element);
double computeGradHy_analytique(const double &x, const double &y, const elementStruct &element);
double computeGradGx_analytique(const double &x, const double &y, const elementStruct &element);
double computeGradGy_analytique(const double &x, const double &y, const elementStruct &element);
double computeG_analytique(const double &x, const double &y, const elementStruct &element);
void electricField(int &nb2dElements, std::vector<int> &entities_tags_2d_domain_BEM, const std::vector<elementStruct> &elementVector, std::vector<size_t> &tags, const std::string &integrationType, std::vector<std::vector<double>> &data_Electric_field, bool &solution_type);
