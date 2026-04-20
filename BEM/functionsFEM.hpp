#include <iostream>
#include <vector> 
#include <map>
#include <Eigen/Dense>
#include <list>
#include <Eigen/Sparse>

// data structure used for the non-linear solver "solverFEMnonLinear"
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

    Eigen::Matrix<double, Eigen::Dynamic, 1> A; // see report for exact definition
    Eigen::Matrix<double, 1, Eigen::Dynamic> G; // see report for exact definition
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P; // see report for exact definition
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> E; // see report for exact definition
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C; // see report for exact definition
};

/*------------ "SOLVER" functions -----------------*/

// solves the non-linear (large displacements, no large deformations) elastic problem.
// if coupled to the BEM solver, "electrostaticPressure" is the pressure computed from the electric field
// on the BEM_FEM boundary (mapping element tag -> pressure). "boundaryDisplacementMap" is filled for all
// nodes on the BEM_FEM_boundary as (node tag -> (u,v) displacement)
// "nbViews" is the current number of views already allocated in the gmsh window.
// "postProcessing" boolean activates the post processing operations (must only be computed once at the end of the
// coupled solver). "iteration" denotes the current iteration of the two-way coupled solver (starting at 1).
// "viewTagU" and "viewTagF" correspond to the view tags of the nodal displacements and the nodal forces in the gmsh window.
// when "untangleMesh" is activated, the nodal coordinates of the nodes in the FEM_domain are updated according to
// their final nodal displacements.
void solverFEMnonLinear(std::map<int, double> &electrostaticPressure, 
                    std::map<int,std::pair<double,double>> & boundaryDisplacementMap, int &nbViews, bool postProcessing, 
                    const int iteration, const int viewTagU, const int viewTagF, bool untangleMesh);

// solves the linear elastic problem in the FEM domain.
// "electrostaticPressure" contains the pressure resulting from the electric field computation on the 
// BEM_FEM_boundary. "nbViews" is the current number of views already allocated in the gmsh window.
void solverFEM(std::map<int, double> &electrostaticPressure, int &nbViews);

/*------------ More specific functions called inside the "SOLVER" functions -----------------*/

// returns the boolean non_linear_solver parameter based on what is set in the .geo file,
// add SetNumber("Non_linear_solver", 0); in .geo file to use linear solver.
bool getNonLinearParameter();

// fills the "FEM_nodeTags" vector with the nodes in the FEM domain, fills the "FEM_nodeCoordsMap" with the coordinates (x,y)
// corresponding to each nodeTag, looking at entity tags "tags" of dimension "dim" corresponding to the FEM domain.
void RetrieveFEMNodeTagsAndCoords(std::vector<std::size_t> & FEM_nodeTags, 
                            std::map<int,std::pair<double,double>> & FEM_nodeCoordsMap, 
                            const int & dim, const std::vector<int> & tags);

// fills the "FEM_elementTags" vector with all element tags of the elements of dimension "dim" (in practice: 2)
// in the entities associated to the entity tags "tags" (corresponding to the FEM domain).
// associates an initialized elementData structure to each element tag through the "FEM_kinematicsMap" map,
// using the "FEM_nodeCoordsMap" map gathering the coordinates of each node in the FEM domain.
int initKinematics(std::vector<std::size_t> & FEM_elementTags, std::map<int, elementData> & FEM_kinematicsMap, 
                std::map<int,std::pair<double,double>> & FEM_nodeCoordsMap, 
                const int & dim, const std::vector<int> & tags);

// updates the kinematics matrices of one element "element" based on its updated nodal positions and rotation angle.
// if "computeFl" is set as true, also computes the local nodal forces vector "fl = Kl*pl", only OK if 
// the local stiffness matrix Kl has already been initialized, using "fillElementalKandF".
void updateKinMatrices(elementData* element, const bool computeFl);

// gets physical properties of FEM_domain: Young's modulus "E", Poisson's ratio "nu",
// density "rho", volumic force along x "bx", volumic force along y "by".
void getPhysicalProperties(double & E, double & nu, double & rho, double & bx, double & by);

// computes the Hooke's matrix in plane stress "H_matrix" based on the Young's modulus "E" and Poisson's ratio "nu".
Eigen::Matrix<double, 3, 3> computeHmatrix(const double E, const double nu);

// looping over every element in the FEM_domain (dimension "dim" and entity tags "tags"), filling "full_f_vector",
// with the local volumic forces computed based on the volumic force ("bx","by") and density "rho" with the help of
// "nodeIndexMap" mapping each nodeTag to its global index.
// filling the local stiffness matrices stored in the elementData structure contained in "FEM_kinematicsMap".
void fillElementalKandF(const Eigen::Matrix<double, 3, 3> H_matrix, const double & rho, const double & bx, 
                const double & by, const int & dim, const std::vector<int> & tags, std::map<int, int> & nodeIndexMap,
                std::vector<double> & full_f_vector, std::map<int, elementData> & FEM_kinematicsMap);

// computing the strain displacement "B_matrix" at the "gaussIndex"-th gaussian point of the element
// consisting of "numNodes" nodes, based on the "basisFunctionsGradient", "determinants" and the inverse of the jacobian
// matrix "jacobinvtrans".
void computeBmatrix(const int &gaussIndex, const int &numNodes, const std::vector<double> & basisFunctionsGradient,
                      const std::vector<double> & determinants, 
                      const Eigen::Matrix<double, 2, 2> & jacobinvtrans,
                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & B_matrix);

// assembles the local "f_vector" into the full f vector, "elIndex" is the element index with respect to "nodeTags"
// having "numNodes" nodes, "nodeTags" gathering all the nodeTags of the current element type in the current entity.
// "nodeIndexMap" associates a global index to each nodetag.
void assembleFlocal(const int &elIndex, const std::vector<std::size_t> & nodeTags, const int &numNodes, 
                std::map<int, int> & nodeIndexMap, std::vector<double> & full_f_vector, const std::vector<double> & f_vector);

// fills the "full_f_vector" with the surfacic contribution of the neumann boundary conditions based on the 
// "nodeIndexMap" associating a global index to each node tag, the local integration of the surfacic forces is 
// performed using Gauss integration based on "integration rule".
// "BCkeys" gathers the keys associated to boundary conditions in the .geo file (using SetNumber(...)).
// "groups" associates a dimension and a tag to each physical group.
void includeNeumannBC(const std::vector<std::string> & BCkeys, const std::string & integration_rule,
                std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups,
                std::vector<double> & full_f_vector);

// for the coupled solver, fills the "full_f_vector" with the electrostatic pressure contained in the 
// "electrostaticPressure" map for different element tags, using Gauss integration "integration_rule", the "nodeIndexMap" 
// associating a global nodal index to each nodeTag, "nodeCoordsMap" gathering the coordinates of each node and "groups"
// containing the different physical groups of the domain.
void applyElecPressure(const std::string & integration_rule, std::map<int,std::pair<double,double>> & nodeCoordsMap, 
                    std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups, 
                    std::map<int, double> &electrostaticPressure, std::vector<double> & full_f_vector);

// looping over the boundary conditions imposed in the .geo file and gathered in "BCkeys" related to the
// physical groups "groups", the goal is to fill the two vectors "new_DOFindices" and "dir_BC":
// "new_DOFindices" vector will contain the new indices of the K matrix after the "number_removed_lines" rows/columns
// corresponding to dirichlet boundary conditions have been removed. A removed index will correspond to "new_DOFindices"= -1.
// dir_BC[index] will contain the dirichlet boundary condition prescribed to the node whose global index is retrieved using
// "nodeIndexMap".
void fillDOFindicesDirBC(const std::vector<std::string> & BCkeys, std::map<int, int> & nodeIndexMap,
                        std::map<std::string, std::pair<int, int>> & groups, 
                        std::vector<int> & new_DOFindices, std::vector<double> & dir_BC, int & number_removed_lines);

// updating the kinematics of the "FEM_elementTags" contained in the "FEM_kinematicsMap",
// based on the "full_d_vector" gathering the nodal displacements associated to the global
// node indices contained in "nodeIndexMap".
void updateKinematics(std::map<int, elementData> & FEM_kinematicsMap, const std::vector<double> & full_d_vector, 
                std::map<int, int> & nodeIndexMap, const std::vector<std::size_t> & FEM_elementTags);

// retrieving the global nodal forces based on local nodal forces of the 
// "FEM_elementTags" contained in the "FEM_kinematicsMap", 
// each node Tag corresponds to a global index given by "nodeIndexMap".
void retrieveTotalNodalForces(std::vector<double> & final_nodal_forces, std::map<int, elementData> & FEM_kinematicsMap, 
                        std::map<int, int> & nodeIndexMap, const std::vector<std::size_t> & FEM_elementTags);

// assembles the "global_f_vector" based on the local nodal forces vector of the 
// "FEM_elementTags" contained in the "FEM_kinematicsMap", 
// without taking into account the DOFs fixed by a dirichlet boundary condition,
// the nodes are renumbered using "new_DOFindices".
// each node Tag corresponds to a global index given by "nodeIndexMap".
void assembleGlobalF(std::vector<double> & global_f_vector, std::map<int, elementData> & FEM_kinematicsMap, 
                std::map<int, int> & nodeIndexMap, const std::vector<int> & new_DOFindices, 
                const std::vector<std::size_t> & FEM_elementTags);

// assembling the global tangent stiffness matrix based on the local elemental stiffness matrices contained
// in the "FEM_kinematicsMap" associated to each "FEM_elementTags".
// "new_DOFindices" maps the global initial node index (from "nodeIndexMap") to the index
// obtained without considering the fixed DOFs.
void assembleGlobalKg(Eigen::SparseMatrix<double> & global_tangent_matrix, std::map<int, elementData> & FEM_kinematicsMap, 
                    std::map<int, int> & nodeIndexMap, const std::vector<int> & new_DOFindices, 
                    const std::vector<std::size_t> & FEM_elementTags);

// updating the "full_d_vector" based on the "global_p_vector", using the global and reduced index mapping:
// "new_DOFindices".
void updateNodalDispVector(std::vector<double> & full_d_vector, const std::vector<double> & global_p_vector, 
                        const std::vector<int> & new_DOFindices);

// retrieving the nodal displacements on the BEM_FEM_boundary in order to communicate them to the
// BEM solver in the two-way coupled solver. Using the "groups" to find the physical group associated
// to the boundary, all nodes of the boundary are associated through the "boundaryDisplacementMap" their
// nodal displacement (u,v) retrieved from "full_d_vector" using the indexation of "nodeIndexMap".
void retrieveBoundaryDisplacements(std::map<int,std::pair<double,double>> & boundaryDisplacementMap, 
                            const std::vector<double> & full_d_vector, std::map<int, int> & nodeIndexMap, 
                            std::map<std::string, std::pair<int, int>> & groups);

// computing the reaction forces on all the physical groups "groups" for which a dirichlet boundary
// condition has been prescribed (using the "BCkeys" to retrieve the groups from the .geo file),
// based on the "final_nodal_forces" and the "nodeIndexMap" to retrieve the global index associated
// to one given node tag.
void computeReactionForces(std::map<std::string, std::pair<int, int>> & groups,
                     std::map<int, int> & nodeIndexMap, const std::vector<double> & final_nodal_forces,
                     const std::vector<std::string> & BCkeys);

// plotting the nodal displacement values associated to "FEM_nodeTags" in the "names[0]" model,
// based on the "full_d_vector".
void plotUxUy(const std::vector<std::size_t> & FEM_nodeTags, const std::vector<double> & full_d_vector, 
            const std::vector<std::string> & names);

// computing the elemental-nodal strains and stresses in the FEM domain (dimension "dim" and entity tags "tags"),
// based on Hooke's matrix in plane stress "H_matrix", Young's modulus "E", Poisson's ratio "nu", looping over the
// "nbelems" elements contained in "FEM_kinematicsMap".
// returning the "max_VM_stress" and the corresponding nodetag "max_VM_nodeTag".
void computeStressStrainElNodal(const int & dim, const std::vector<int> & tags, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           const double nu, const double E, std::map<int, elementData> & FEM_kinematicsMap,
                           const int & nbelems, double & max_VM_stress, int & max_VM_nodeTag);

// computing the work done by external forces, as the nodal forces are energy consistent, it
// is sufficient to compute the "externalWork" as the dot product of the nodal displacements
// "full_d_vector" and the nodal forces "final_nodal_forces".
void computeWorkDoneByExternalForces(double &externalWork, const std::vector<double> & full_d_vector, 
                                const std::vector<double> & final_nodal_forces);

// computing the total (internal) strain energy "strain_energy" in the whole FEM domain (dimension "dim" and 
// entity tags "tags"), using the elementData stored in "FEM_kinematicsMap" gathering the informations for all
// FEM elements. The Hooke's matrix in plane stress "H_matrix" is used to compute the stresses from the strains.
// Computed as a sum of integrals over the finite elements, using Gauss integration "integration_rule".
void computeStrainEnergy(double & strainEnergy, const int & dim, const std::vector<int> & tags,
                           std::map<int, elementData> & FEM_kinematicsMap, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           const std::string & integration_rule);

//displays the time of the desired code section "function_name", if "display_time" is set as true.
void displayTime(const double &start, const double &end, const std::string &function_name, const bool &display_time);

/*-------------- FUNCTIONS SPECIFIC TO LINEAR FEM ------------------*/

// assembling the full stiffness matrix "full_K_matrix" using first a triplet vector, assembling the volumic part 
// of the nodal forces "full_f_vector", in the FEM_domain (dimension "dim" and entity tag "tags") using the global
// node indexation defined by "nodeIndexMap". "H_matrix" is the plane stress Hooke's matrix, "rho" the density of the
// material, "bx" the volumic force along x, "by" the volumic force along y.
void assembleKVolumicF(const Eigen::Matrix<double, 3, 3> & H_matrix, const double & rho, const double & bx, const double & by,
                const int & dim, const std::vector<int> & tags, std::map<int, int> & nodeIndexMap, 
                Eigen::SparseMatrix<double> & full_K_matrix, std::vector<double> & full_f_vector);

// assembles the elemental "K_matrix" into a vector of triplets "k_tripletList" using 
// the node tags provided by "nodeTags" vector (index of the element is "elIndex").
// The element has "numNodes" nodes. "nodeIndexMap" provides the global nodal index corresponding to a given node tag.
void assembleKtriplet(const int &elIndex, const std::vector<std::size_t> & nodeTags, const int numNodes, 
                    std::map<int, int> & nodeIndexMap, std::vector<Eigen::Triplet <double>> & k_tripletList,
                    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & K_matrix);

// the "full_f_vector" and "full_K_matrix" are reduced in "reduced_f_vector" and "reduced_K_matrix" respectively,
// taking into account the dirichlet boundary conditions contained in "dir_BC".
// "new_DOFindices" gathers the new indices of the DOFs when taking the boundary conditions into account.
void fillReducedKf(const std::vector<double> & full_f_vector, const Eigen::SparseMatrix<double> & full_K_matrix,
                const std::vector<int> & new_DOFindices, const std::vector<double> & dir_BC,
                Eigen::SparseMatrix<double> & reduced_K_matrix, std::vector<double> & reduced_f_vector);

// computing the elemental-nodal strains and stresses in the FEM domain (dimension "dim" and entity tags "tags"),
// based on Hooke's matrix in plane stress "H_matrix", Young's modulus "E", Poisson's ratio "nu", and the full nodal
// displacement vector "full_d_vector" and the global nodal index map "nodeIndexMap".
// "nbelems" is the number of 2D elements contained in the FEM domain.
// returning the "max_VM_stress" and the corresponding nodetag "max_VM_nodeTag".
void computeStressStrainLinearWay(const int & dim, const std::vector<int> & tags, const std::vector<double> &full_d_vector,
                           const Eigen::Matrix<double, 3, 3> &H_matrix, const double &nu, const double &E,
                           std::map<int, int> & nodeIndexMap, const int &nbelems, double &max_VM_stress, int &max_VM_nodeTag);

// computing the total (internal) strain energy "strain_energy" in the whole FEM domain (dimension "dim" and 
// entity tags "tags"), using the nodal displacements in "full_d_vector". The nodes are indexed using the "nodeIndexMap".
// The Hooke's matrix in plane stress "H_matrix" is used to compute the stresses from the strains.
// Computed as a sum of integrals over the finite elements, using Gauss integration "integration_rule".
void computeStrainEnergyLinearWay(double & strainEnergy, const int & dim, const std::vector<int> & tags,
                           const std::vector<double> &full_d_vector, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           std::map<int, int> & nodeIndexMap, const std::string & integration_rule);