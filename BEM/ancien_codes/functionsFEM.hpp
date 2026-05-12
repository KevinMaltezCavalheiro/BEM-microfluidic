#include <iostream>
#include <vector>   // Why should I include this ?
#include <map>
#include <Eigen/Dense>
#include <list>
#include <Eigen/Sparse> // Eigen library for sparse matrices

struct elementData{
    // will store all the kinematics data
    std::vector<size_t> nodes; // will contain the tags of the nodes
    std::vector<double> x; // will contain the x global coordinate of the nodes
    std::vector<double> y; // will contain the y global coordinate of the nodes

    double xc; // x global coordinate of center
    double yc; // y global coordinate of center

    std::vector<double> xi; // x local coordinate of nodes
    std::vector<double> yi; // y local coordinate of nodes

    std::vector<double> u; // x global displacement of the nodes
    std::vector<double> v; // y global displacement of the nodes

    double uc; // x global displacement of center
    double vc; // y global displacement of center

    double theta; // angle of rotation of element

    std::vector<double> ui; // x local displacement of nodes
    std::vector<double> vi; // y local displacement of nodes

    Eigen::Matrix<double, Eigen::Dynamic, 1> pl; //local displacement under vector form

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Kl; // local K matrix

    Eigen::Matrix<double, Eigen::Dynamic, 1> fl; //local current (iterative) forces under vector form

    Eigen::Matrix<double, Eigen::Dynamic, 1> A;
    Eigen::Matrix<double, 1, Eigen::Dynamic> G;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> E;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B;
};

Eigen::Matrix<double, 3, 3> computeHmatrix(const double E, const double nu);
void assembleFthread(const int i, const int el, const std::vector<std::vector<std::size_t>> & nodeTags, const int numNodes, std::map<int, int> & nodeIndexMap, std::list<std::pair<int,double>> & thread_f_vector, std::vector<double> & f_vector);
void computeBmatrix(const int j, const int numNodes, const std::vector<double> & basisFunctionsGradient, const std::vector<double> & determinants, const Eigen::Matrix<double, 2, 2> & jacobinvtrans, const Eigen::Matrix<double, 3, 3> & H_matrix, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & B_matrix);
void getPhysicalProperties(const std::vector<std::string> & keys, double & E, double & nu, double & rho, double & bx, double & by);
void fillElementalKandF(const Eigen::Matrix<double, 3, 3> H_matrix, const double rho, const double bx, const double by, const std::string & integration_rule, int & dim, std::vector<int> & tags, std::map<int, int>  & nodeIndexMap, std::vector<std::list<Eigen::Triplet <double>>> & my_k_list, std::vector<std::list<std::pair<int,double>>> & thread_f_vector, std::map<int, elementData> & FEM_kinematicsMap);
void computeNeumannBC(const int numNodes, const int i, const std::vector<std::vector<std::size_t>> & elementTags, const double tx, const double ty, std::map<int, int> & nodeIndexMap, std::vector<double> & f_vector, const std::vector<double> & weights, const std::vector<double> & basisFunctions, const std::vector<double> & determinants, const std::vector<std::vector<std::size_t>> & nodeTags, std::vector<double> & full_f_vector);
void includeNeumannBC(const std::vector<std::string> & keys, const std::string & integration_rule, std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups, std::vector<double> & full_f_vector);
void fillDOFindicesDirBC(const std::vector<std::string> & keys, std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups, std::vector<int> & new_DOFindices, std::vector<double> & dir_BC);
void computeReactionForces(std::map<std::string, std::pair<int, int>> & groups, std::map<int, int> & nodeIndexMap, std::vector<double> & nodal_forces, const std::vector<std::string> & keys);
void plotUxUy(const std::vector<std::size_t> & FEM_nodeTags, const std::vector<double> & full_d_vector, const int nbViews, std::vector<std::string> & names);
void computeStressStrainElNodal(const std::vector<int> & elementTypes, const std::vector<std::vector<std::size_t>> &elementTags, const std::vector<std::vector<std::size_t>> &elnodeTags, std::vector<std::size_t> & elems2D, const Eigen::Matrix<double, 3, 3> &H_matrix, const double nu, const double E, std::map<int, elementData> & FEM_kinematicsMap, std::vector<std::vector<double>> & data3, std::vector<std::vector<double>> & data4, std::vector<std::vector<double>> & data5, std::vector<std::vector<double>> & data6, std::vector<std::vector<double>> & data7, std::vector<std::vector<double>> & data8, std::vector<std::vector<double>> & data9, std::vector<std::vector<double>> & data10, int & k);
void applyElecPressure(const std::string & integration_rule, std::map<int,std::pair<double,double>> & nodeCoordsMap, std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups, std::map<int, double> &electrostaticPressure, std::vector<double> & full_f_vector);
int initKinematics(std::vector<std::size_t> & FEM_elementTags, std::map<int, elementData> & FEM_kinematicsMap, std::map<int,std::pair<double,double>> & FEM_nodeCoordsMap, int & dim, std::vector<int> & tags);
void updateKinematics(std::map<int, elementData> & FEM_kinematicsMap, std::vector<double> & full_d_vector, std::map<int, int> & nodeIndexMap, std::vector<std::size_t> & FEM_elementTags);
void updateKinMatrices(elementData* element, bool computeFl);
void Assembly_global_f(std::vector<double> & global_f_vector, std::map<int, elementData> & FEM_kinematicsMap, std::map<int, int> & nodeIndexMap, std::vector<int> & new_DOFindices, std::vector<std::size_t> & FEM_elementTags);
double computeResidualNorm(std::vector<double> & global_f_vector, std::vector<double> & f_app_vector, Eigen::Matrix<double, Eigen::Dynamic, 1> & residual);
void Assembly_global_Kg(Eigen::SparseMatrix<double> & global_tangent_matrix, std::map<int, elementData> & FEM_kinematicsMap, std::map<int, int> &  nodeIndexMap, std::vector<int> & new_DOFindices, std::vector<std::size_t> & FEM_elementTags);
void updateNodalDispVector(std::vector<double> & full_d_vector, std::vector<double> & global_p_vector, std::vector<int> & new_DOFindices);
void retrieveFinalNodalForces(std::vector<double> & final_nodal_forces, std::map<int, elementData> & FEM_kinematicsMap, std::map<int, int> & nodeIndexMap, std::vector<std::size_t> & FEM_elementTags);
void retrieveBoundaryDisplacements(std::map<int,std::pair<double,double>> & boundaryDisplacementMap, const std::vector<double> & full_d_vector, std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups);

int solverFEM(int argc, char **argv, std::map<int, double> &electrostaticPressure, std::map<int,std::pair<double,double>> & boundaryDisplacementMap, const int nbViews, bool postProcessing, const int iteration, const int viewTagU, const int viewTagF);