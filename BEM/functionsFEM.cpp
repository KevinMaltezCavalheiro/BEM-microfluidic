#ifdef _MSC_VER
#include <gmsh.h_cwrap>
#else
#include <gmsh.h>
#endif

#include "functionsFEM.hpp"

#include <iostream> 
#include <map>
#include <sstream>
#include <Eigen/Dense> 
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky> // solving sparse linear systems
#include <Eigen/Core>
#include <cmath>
#include <algorithm> 
#ifdef _OPENMP
#include <omp.h>
#endif
#define PI M_PI

// returns the boolean non_linear_solver parameter based on what is set in the .geo file,
// add SetNumber("Non_linear_solver", 0); in .geo file to use linear solver.
bool getNonLinearParameter()
{
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "Non_linear_solver");
    if(keys.size())
    {
        std::vector<double> value;
        gmsh::onelab::getNumber(keys[0], value);
        if(value[0])
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        std::cout << "you did not choose the type of solver (linear or not) in the .geo file, 'SetNumber('Non_linear_solver',x);'. Automatically set to non linear type \n";
        return true;
    }
}

// fills the "FEM_nodeTags" vector with the nodes in the FEM domain, fills the "FEM_nodeCoordsMap" with the coordinates (x,y)
// corresponding to each nodeTag, looking at entity tags "tags" of dimension "dim" corresponding to the FEM domain.
void RetrieveFEMNodeTagsAndCoords(std::vector<std::size_t> & FEM_nodeTags, 
                            std::map<int,std::pair<double,double>> & FEM_nodeCoordsMap, 
                            const int & dim, const std::vector<int> & tags)
{
    //loop over the entities of the FEM domain to retrieve all nodeTags by taking care not to take some nodes twice
    //using sort and set_union methods + retrieving the node coordinates and storing them in a map
    std::vector<double> FEM_nodecoord;
    std::vector<double> FEM_nodeparametricCoord;
    std::vector<std::vector<std::size_t>> tmp_nodeTags(tags.size()); // will store the nodeTags associated to one given entity of the domain
    gmsh::model::mesh::getNodes(tmp_nodeTags[0], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[0], true); 
    for(std::size_t j = 0; j < tmp_nodeTags[0].size(); j++)
    {
        FEM_nodeCoordsMap[tmp_nodeTags[0][j]] = std::pair<double, double>(FEM_nodecoord[3*j], FEM_nodecoord[3*j+1]);
    }
    std::sort(tmp_nodeTags[0].begin(), tmp_nodeTags[0].end()); // need to sort the vectors in order to merge them afterwards
    for (std::size_t i = 1; i < tags.size(); i++)
    { // in the case the domain contains multiple entities
        gmsh::model::mesh::getNodes(tmp_nodeTags[i], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[i], true);
        for(std::size_t j = 0; j < tmp_nodeTags[i].size(); j++)
        {
            FEM_nodeCoordsMap[tmp_nodeTags[i][j]] = std::pair<double, double>(FEM_nodecoord[3*j], FEM_nodecoord[3*j+1]);
        }
        std::sort(tmp_nodeTags[i].begin(), tmp_nodeTags[i].end());
        std::set_union(tmp_nodeTags[i-1].begin(), tmp_nodeTags[i-1].end(),
                       tmp_nodeTags[i].begin(), tmp_nodeTags[i].end(), std::back_inserter(FEM_nodeTags)); // works only for 2 domains -> see little trick just below
        
    }
    if(tags.size() == 1){ // if the domain contains one single entity
        FEM_nodeTags = tmp_nodeTags[0];
    }
    else{ // little trick for more than 2 domains
        std::sort(FEM_nodeTags.begin(), FEM_nodeTags.end());
        auto last = std::unique(FEM_nodeTags.begin(), FEM_nodeTags.end());
        FEM_nodeTags.erase(last, FEM_nodeTags.end());
    }
}

// fills the "FEM_elementTags" vector with all element tags of the elements of dimension "dim" (in practice: 2)
// in the entities associated to the entity tags "tags" (corresponding to the FEM domain).
// associates an initialized elementData structure to each element tag through the "FEM_kinematicsMap" map,
// using the "FEM_nodeCoordsMap" map gathering the coordinates of each node in the FEM domain.
int initKinematics(std::vector<std::size_t> & FEM_elementTags, std::map<int, elementData> & FEM_kinematicsMap, 
                std::map<int,std::pair<double,double>> & FEM_nodeCoordsMap, 
                const int & dim, const std::vector<int> & tags)
{
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;

    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    
    int nbelems=0; // will store total number of elements in FEM domain
    for(std::size_t i=0; i< tags.size(); ++i)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(std::size_t j=0; j< elementTags.size(); ++j)
        {
            nbelems+=elementTags[j].size();
            // filling the FEM_elementTags vector
            FEM_elementTags.insert(FEM_elementTags.end(), elementTags[j].begin(), elementTags[j].end());

            // getting element properties, in particular: number of nodes
            gmsh::model::mesh::getElementProperties(elementTypes[j],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
            for(std::size_t k = 0; k < elementTags[j].size(); k++)
            {
                elementData element;
                std::vector<size_t> tmp_nodes(numNodes);
                std::vector<double> tmp_x(numNodes);
                std::vector<double> tmp_y(numNodes);
                std::vector<double> tmp_u(numNodes);
                std::vector<double> tmp_v(numNodes);
                std::vector<double> tmp_ui(numNodes);
                std::vector<double> tmp_vi(numNodes);
                double tmp_xc = 0; 
                double tmp_yc = 0;
                for(int l = 0; l < numNodes; l++)
                {
                    tmp_nodes[l] = elnodeTags[j][numNodes*k + l];
                    tmp_x[l] = FEM_nodeCoordsMap[tmp_nodes[l]].first;
                    tmp_y[l] = FEM_nodeCoordsMap[tmp_nodes[l]].second;
                    tmp_u[l] = 0; // initializing the displacement at zero
                    tmp_v[l] = 0;
                    tmp_ui[l] = 0;
                    tmp_vi[l] = 0;
                    tmp_xc += tmp_x[l];
                    tmp_yc += tmp_y[l];
                }
                tmp_xc = tmp_xc / numNodes;
                tmp_yc = tmp_yc / numNodes;

                std::vector<double> tmp_xi(numNodes);
                std::vector<double> tmp_yi(numNodes);
                for(int l = 0; l < numNodes; l++)
                {
                    tmp_xi[l] = tmp_x[l] - tmp_xc;
                    tmp_yi[l] = tmp_y[l] - tmp_yc;
                }
                element.nodes = tmp_nodes;
                element.x = tmp_x;
                element.y = tmp_y;
                element.xc = tmp_xc;
                element.yc = tmp_yc;
                element.xi = tmp_xi;
                element.yi = tmp_yi;
                element.u = tmp_u;
                element.v = tmp_v;
                element.uc = 0;
                element.vc = 0;
                element.theta = 0;
                element.ui = tmp_ui;
                element.vi = tmp_vi;
                element.pl = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(2*numNodes,1);
                element.A = Eigen::Matrix<double, Eigen::Dynamic, 1>(2*numNodes,1);
                element.G = Eigen::Matrix<double, 1, Eigen::Dynamic>(1,2*numNodes);
                element.E = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(2*numNodes, 2*numNodes);

                FEM_kinematicsMap[elementTags[j][k]] = element;
                updateKinMatrices(& FEM_kinematicsMap[elementTags[j][k]], false);
            }
        }
    }
    return nbelems;
}

// updates the kinematics matrices of one element "element" based on its updated nodal positions and rotation angle.
// if "computeFl" is set as true, also computes the local nodal forces vector "fl = Kl*pl", only OK if 
// the local stiffness matrix Kl has already been initialized, using "fillElementalKandF".
void updateKinMatrices(elementData* element, const bool computeFl)
{
    // for the different formulas, please refer to the report
    size_t numNodes = element->nodes.size();
    // update A matrix
    for(size_t i = 0; i < numNodes; i++)
    {
        element->A(2*i, 0) = -(element->vi[i]+element->yi[i]);
        element->A(2*i + 1, 0) = element->ui[i]+element->xi[i];
    }

    // update G matrix
    double denominator = 0;
    for(size_t i = 0; i < numNodes; i++)
    {
        element->G(0, 2*i) = -element->yi[i];
        element->G(0, 2*i + 1) = element->xi[i];
        denominator += element->xi[i]*(element->ui[i]+element->xi[i]) + element->yi[i]*(element->vi[i]+element->yi[i]);
    }
    element->G *= 1 / denominator;

    // update P matrix 
    element->P = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(2*numNodes,2*numNodes) - element->A*element->G;

    // update E matrix
    double tmp_cos = cos(element->theta);
    double tmp_sin = sin(element->theta);
    for(size_t i = 0; i < numNodes; i++)
    {
        element->E(2*i, 2*i) = tmp_cos;
        element->E(2*i, 2*i + 1) = -tmp_sin;
        element->E(2*i + 1, 2*i) = tmp_sin;
        element->E(2*i + 1, 2*i + 1) = tmp_cos;
    }

    //update C matrix
    element->C = element->P*element->E.transpose();

    //update fl vector, only if Kl has been initialized
    if(computeFl)
    {
        element->fl = element->Kl*element->pl;
    }
}

// gets physical properties of FEM_domain: Young's modulus "E", Poisson's ratio "nu",
// density "rho", volumic force along x "bx", volumic force along y "by".
void getPhysicalProperties(double & E, double & nu, double & rho, double & bx, double & by)
{
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Volumic Forces|Materials).+");
    for (auto &key : keys)
    {
        // get corresponding value
        std::vector<double> value;
        gmsh::onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("Young") == 0){
                E = value[0];
            }
            if(words[2].compare("Poisson") == 0){
                nu = value[0];
            }
            if(words[2].compare("rho") == 0){
                rho = value[0];
            }
            if(words[2].compare("bx") == 0){
                bx = value[0];
            }
            if(words[2].compare("by") == 0){
                by = value[0];
            }
        }
    }
}

// computes the Hooke's matrix in plane stress "H_matrix" based on the Young's modulus "E" and Poisson's ratio "nu".
Eigen::Matrix<double, 3, 3> computeHmatrix(const double E, const double nu)
{
    Eigen::Matrix<double, 3, 3> H_matrix;
    H_matrix(0,0) = 1.;
    H_matrix(1,0) = nu;
    H_matrix(2,0) = 0.;
    H_matrix(0,1) = nu;
    H_matrix(1,1) = 1.;
    H_matrix(2,1) = 0.;
    H_matrix(0,2) = 0.;
    H_matrix(1,2) = 0.;
    H_matrix(2,2) = (1-nu)/2;
    H_matrix = E/(1-nu*nu)*H_matrix;
    return H_matrix;
}

// looping over every element in the FEM_domain (dimension "dim" and entity tags "tags"), filling "full_f_vector",
// with the local volumic forces computed based on the volumic force ("bx","by") and density "rho" with the help of
// "nodeIndexMap" mapping each nodeTag to its global index.
// filling the local stiffness matrices stored in the elementData structure contained in "FEM_kinematicsMap".
void fillElementalKandF(const Eigen::Matrix<double, 3, 3> H_matrix, const double & rho, const double & bx, 
                const double & by, const int & dim, const std::vector<int> & tags, std::map<int, int> & nodeIndexMap,
                std::vector<double> & full_f_vector, std::map<int, elementData> & FEM_kinematicsMap)
{
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;
    std::vector<double> basisFunctionsGradient;

    for (std::size_t k = 0; k < tags.size(); ++k)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types
        for (std::size_t i = 0; i < elementTypes.size(); ++i)
        {
            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, elDim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);

            // adapting the composite integration rule as a function of the order of the elements
            std::string adapting_integration_rule = "CompositeGauss14";
            if(order < 8)
            {
                adapting_integration_rule = "CompositeGauss" + std::to_string(2*order);
            }
            
            // get Gauss points coordinates and weights
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    adapting_integration_rule, 
                                                    localCoords, weights);

            // get determinants and Jacobians
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords, tags[k]);

            // values and derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);

            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctionsGradient, 
                                                 numOrientations);

            /*---------computing the K matrix (and f vector, only volumic part for now)------*/
            // looping over the elements (computing elemental K matrices)
            #ifdef _MSC_VER
            // https://github.com/tesseract-ocr/tesseract/issues/3285
            // OpenMP 2.0 standard is enforced - we may use /openmp:llvm and keep size_t
            #pragma omp parallel for
            for (size_t el = 0; el < elementTags[i].size(); el++)
            {
            #else
            #pragma omp parallel for
            for (std::size_t el = 0; el < elementTags[i].size(); el++)
            {
            #endif
                // elemental K matrix and f vector initialization
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K_matrix(2*numNodes,2*numNodes);
                std::vector<double> f_vector(2*numNodes);
                K_matrix.setZero();
                for (int j = 0; j < 2*numNodes; ++j)f_vector[j] = 0.;

                // looping over the Gauss points
                for (std::size_t gaussIndex = 0; gaussIndex < weights.size(); ++gaussIndex)
                {
                    // converting jacobian to Eigen format
                    Eigen::Matrix<double, 2, 2> jacob2D;
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*gaussIndex+9*weights.size()*el], 3, 3);
                    // 3D to 2D (plane stress case)
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    Eigen::Matrix<double, 2, 2> jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    Eigen::Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix (formula slide 19) and the f vector (also using gauss integration)
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                    computeBmatrix(gaussIndex, numNodes, basisFunctionsGradient, determinants, jacobinvtrans, B_matrix);
                    
                    for (int l = 0; l < numNodes; ++l)
                    { // looping over the number of shape functions
                        f_vector[2*l] = f_vector[2*l] + rho*basisFunctions[l+numNodes*gaussIndex]*bx*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                        f_vector[2*l+1] = f_vector[2*l+1] + rho*basisFunctions[l+numNodes*gaussIndex]*by*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                    }

                    K_matrix = K_matrix + B_matrix.transpose()*H_matrix*B_matrix*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                FEM_kinematicsMap[elementTags[i][el]].Kl = K_matrix;
                assembleFlocal(el, nodeTags[i], numNodes, nodeIndexMap, full_f_vector, f_vector);
            }
        }
    }
}

// computing the strain displacement "B_matrix" at the "gaussIndex"-th gaussian point of the element
// consisting of "numNodes" nodes, based on the "basisFunctionsGradient", "determinants" and the inverse of the jacobian
// matrix "jacobinvtrans".
void computeBmatrix(const int &gaussIndex, const int &numNodes, const std::vector<double> & basisFunctionsGradient,
                      const std::vector<double> & determinants, 
                      const Eigen::Matrix<double, 2, 2> & jacobinvtrans,
                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & B_matrix)
{
    for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
        B_matrix(0,2*l) = jacobinvtrans(0,0)*basisFunctionsGradient[3*l+numNodes*3*gaussIndex]+jacobinvtrans(0,1)*basisFunctionsGradient[3*l+numNodes*3*gaussIndex+1];
        B_matrix(1,2*l) = 0.;
        B_matrix(2,2*l) = jacobinvtrans(1,0)*basisFunctionsGradient[3*l+numNodes*3*gaussIndex]+jacobinvtrans(1,1)*basisFunctionsGradient[3*l+numNodes*3*gaussIndex+1];
        B_matrix(0,2*l+1) = 0.;
        B_matrix(1,2*l+1) = B_matrix(2,2*l);
        B_matrix(2,2*l+1) = B_matrix(0,2*l);
    }
}

// assembles the local "f_vector" into the full f vector, "elIndex" is the element index with respect to "nodeTags"
// having "numNodes" nodes, "nodeTags" gathering all the nodeTags of the current element type in the current entity.
// "nodeIndexMap" associates a global index to each nodetag.
void assembleFlocal(const int &elIndex, const std::vector<std::size_t> & nodeTags, const int &numNodes, 
                std::map<int, int> & nodeIndexMap, std::vector<double> & full_f_vector, const std::vector<double> & f_vector)
{
    int first_node_index;
    int first_node;
    int my_row;
    for (int j = 0; j < 2*numNodes; ++j)
    {
        first_node_index = (int) j/2;
        first_node = nodeTags[numNodes*elIndex+first_node_index];
        my_row = nodeIndexMap[first_node] + j % 2; 
        #pragma omp atomic update
        full_f_vector[my_row] +=  f_vector[j];
    }
}

// fills the "full_f_vector" with the surfacic contribution of the neumann boundary conditions based on the 
// "nodeIndexMap" associating a global index to each node tag, the local integration of the surfacic forces is 
// performed using Gauss integration based on "integration rule".
// "BCkeys" gathers the keys associated to boundary conditions in the .geo file (using SetNumber(...)).
// "groups" associates a dimension and a tag to each physical group.
void includeNeumannBC(const std::vector<std::string> & BCkeys, const std::string & integration_rule,
                std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups,
                std::vector<double> & full_f_vector)
{
    double tx = 0; // surface traction along x axis
    double ty = 0; // surface traction along y axis

    int dim;
    int tag;
    std::vector<int> tags;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;

    for (auto &key : BCkeys)
    {
        // get corresponding value
        std::vector<double> value;
        gmsh::onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            std::string groupname = ""; // will remain empty only if tx == ty == 0
            if(words[2].compare("tx") == 0){
                tx = value[0];
                std::string new_key = "Boundary Conditions/" + words[1] + "/ty";
                gmsh::onelab::getNumber(new_key, value);
                ty = value[0];
                if(tx || ty) // there is no need the consider Neumann BC if it is homogeneous
                    groupname = words[1];
            }
            //here, if we are dealing with neumann BC's, we should have retrieved a non-empty groupname
            if(groupname.compare("") != 0){
                dim = groups[groupname].first;
                tag = groups[groupname].second;
                if (tag == 0)
                {
                    std::cerr << "Group '" << groupname << "' does not exist!\n";
                    return;
                }
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags); // entities of particular groupname

                // loop over entities
                for (std::size_t k = 0; k < tags.size(); ++k)
                {
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (std::size_t i = 0; i < elementTypes.size(); ++i)
                    {
                        // getting element properties
                        gmsh::model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                    
                        // get Gauss points coordinates and weights
                        gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                                integration_rule,
                                                                localCoords, weights);

                        // get determinants and Jacobians (only determinants will be useful)
                        std::vector<double> jacobians, determinants, coords;
                        gmsh::model::mesh::getJacobians(elementTypes[i], 
                                                        localCoords, jacobians, 
                                                        determinants, coords, tags[k]);

                        gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                            "Lagrange", numComponents, 
                                                            basisFunctions, 
                                                            numOrientations);

                        //looping over the element types
                        #pragma omp parallel for
                        for (std::size_t el = 0; el < elementTags[i].size(); el++){
                            std::vector<double> f_vector(2*numNodes);
                            for (int j = 0; j < 2*numNodes; ++j){
                                f_vector[j] = 0.;
                            }
                            // looping over the Gauss points
                            for (std::size_t j = 0; j < weights.size(); ++j){
                                // looping over the number of shape functions
                                for (int l = 0; l < numNodes; ++l){ 
                                    f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*j]*tx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
                                    f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*j]*ty*determinants[j+weights.size()*el]*weights[j]; // le weight non plus !
                                }
                            }
                            /*------------- ASSEMBLY PROCESS ------------*/
                            for (int j = 0; j < 2*numNodes; ++j){
                                int first_node_index = (int) j/2;
                                int first_node = nodeTags[i][numNodes*el+first_node_index];
                                int my_row = nodeIndexMap[first_node] + j % 2;
                                #pragma omp atomic update
                                full_f_vector[my_row] += f_vector[j];
                            }
                        }    
                    }
                }
            }
        }
    }
}

// for the coupled solver, fills the "full_f_vector" with the electrostatic pressure contained in the 
// "electrostaticPressure" map for different element tags, using Gauss integration "integration_rule", the "nodeIndexMap" 
// associating a global nodal index to each nodeTag, "nodeCoordsMap" gathering the coordinates of each node and "groups"
// containing the different physical groups of the domain.
void applyElecPressure(const std::string & integration_rule, std::map<int,std::pair<double,double>> & nodeCoordsMap, 
                    std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups, 
                    std::map<int, double> &electrostaticPressure, std::vector<double> & full_f_vector)
{
    // electrostatic pressure must only be applied on the BEM_FEM_boundary
    std::string groupname = "BEM_FEM_boundary";
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

    int first_nodeTag;
    int second_nodeTag;
    double nx;
    double ny;
    double p;
    double tx;
    double ty;
    double norm;
    int first_node_index;
    int first_node;
    int my_row;

    for(std::size_t i=0; i< tags.size(); ++i)
    {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> elnodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(std::size_t j=0; j< elementTypes.size(); ++j)
        {

            // getting element properties
            std::string elementName;
            int elDim, order, numNodes, numPrimaryNodes;
            std::vector<double> localNodeCoord;
            gmsh::model::mesh::getElementProperties(elementTypes[j],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                    
            // get Gauss points coordinates and weights
            std::vector<double> localCoords, weights;
            gmsh::model::mesh::getIntegrationPoints(elementTypes[j], 
                                                    integration_rule,
                                                    localCoords, weights);

            // get determinants and Jacobians (only determinants will be useful)
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[j], 
                                            localCoords, jacobians, 
                                            determinants, coords,tags[i]);

            int numComponents, numOrientations;
            std::vector<double> basisFunctions;
            gmsh::model::mesh::getBasisFunctions(elementTypes[j], localCoords, 
                                                "Lagrange", numComponents, 
                                                basisFunctions, 
                                                numOrientations);

            #pragma omp parallel for private(first_nodeTag, second_nodeTag, nx, ny, norm, p, tx, ty, first_node_index, first_node)
            for(std::size_t k=0; k < elementTags[j].size(); ++k)
            {
                std::vector<double> f_vector(2*numNodes);
                // computing the outward normal to the element, assuming the BEM domain is traveled "aire à gauche", hence
                // the FEM domain "aire à droite"
                first_nodeTag = elnodeTags[j][numNodes*k];
                second_nodeTag = elnodeTags[j][numNodes*k+1];
                nx = -(nodeCoordsMap[second_nodeTag].second - nodeCoordsMap[first_nodeTag].second);
                ny = (nodeCoordsMap[second_nodeTag].first - nodeCoordsMap[first_nodeTag].first);
                norm = sqrt(nx*nx+ny*ny);
                nx = nx / norm; //normalisation
                ny = ny / norm;
                p = electrostaticPressure[elementTags[j][k]];
                tx = p*nx; // projecting the pressure onto the outward normal
                ty = p*ny;

                for (int n = 0; n < 2*numNodes; ++n){
                    f_vector[n] = 0.;
                }
                // looping over the Gauss points
                for (std::size_t n = 0; n < weights.size(); ++n){
                    // looping over the number of shape functions
                    for (int l = 0; l < numNodes; ++l){ 
                        f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*n]*tx*determinants[n+weights.size()*k]*weights[n];
                        f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*n]*ty*determinants[n+weights.size()*k]*weights[n];
                    }
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                for (int n = 0; n < 2*numNodes; ++n){
                    first_node_index = (int) n/2;
                    first_node = elnodeTags[j][numNodes*k+first_node_index];
                    my_row = nodeIndexMap[first_node] + n % 2;
                    #pragma omp atomic update
                    full_f_vector[my_row] += f_vector[n];
                }
            }
        }
    }
}

// looping over the boundary conditions imposed in the .geo file and gathered in "BCkeys" related to the
// physical groups "groups", the goal is to fill the two vectors "new_DOFindices" and "dir_BC":
// "new_DOFindices" vector will contain the new indices of the K matrix after the "number_removed_lines" rows/columns
// corresponding to dirichlet boundary conditions have been removed. A removed index will correspond to "new_DOFindices"= -1.
// dir_BC[index] will contain the dirichlet boundary condition prescribed to the node whose global index is retrieved using
// "nodeIndexMap".
void fillDOFindicesDirBC(const std::vector<std::string> & BCkeys, std::map<int, int> & nodeIndexMap,
                        std::map<std::string, std::pair<int, int>> & groups, 
                        std::vector<int> & new_DOFindices, std::vector<double> & dir_BC, int & number_removed_lines)
{

    for (std::size_t j = 0; j < new_DOFindices.size(); ++j){
        new_DOFindices[j] = (int) j;
    } 
    for (std::size_t j = 0; j < dir_BC.size(); ++j){
        dir_BC[j] = 0;
    }

    double ux; // imposed displacement along x axis
    double uy; // imposed displacement along y axis

    int dim;
    int tag;
    std::vector<int> tags;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;

    int first_node_index;
    int first_node;
    int my_row;

    for (auto &key : BCkeys)
    {
        // get corresponding value
        std::vector<double> value;
        gmsh::onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("ux") == 0 || words[2].compare("uy") == 0){
                ux = 0;
                uy = 0;
                if(words[2].compare("ux") == 0)
                    ux = value[0];
                if(words[2].compare("uy") == 0)
                    uy = value[0];
                std::string groupname = words[1];
                
                //loop over elements of group_name
                dim = groups[groupname].first;
                tag = groups[groupname].second;
                if (tag == 0)
                {
                    std::cerr << "Group '" << groupname << "' does not exist!\n";
                    return;
                }
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

                for (std::size_t k = 0; k < tags.size(); ++k)
                {
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (std::size_t i = 0; i < elementTypes.size(); ++i)
                    {
                        gmsh::model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                        //looping over the elements
                        for (std::size_t el = 0; el < elementTags[i].size(); el++){
                            for (int j = 0; j < 2*numNodes; ++j){
                                first_node_index = (int) j/2;
                                first_node = nodeTags[i][numNodes*el+first_node_index];
                                my_row = nodeIndexMap[first_node] + j % 2; // j % 2 == 1 corresponds to y nodal displacement
                                if(words[2].compare("uy") == 0 && j % 2){ // we are dealing with a uy type of node
                                    dir_BC[my_row] = uy;
                                    new_DOFindices[my_row] = -1;
                                }
                                if(words[2].compare("ux") == 0 && j % 2 == 0){ // ux type of node
                                    dir_BC[my_row] = ux;
                                    new_DOFindices[my_row] = -1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (std::size_t j = 0; j < new_DOFindices.size(); ++j){
        if(new_DOFindices[j] < 0){ // line has to be removed
            number_removed_lines++;
        }
        else{
            new_DOFindices[j] -= number_removed_lines;
        }
    }
}

// updating the kinematics of the "FEM_elementTags" contained in the "FEM_kinematicsMap",
// based on the "full_d_vector" gathering the nodal displacements associated to the global
// node indices contained in "nodeIndexMap".
void updateKinematics(std::map<int, elementData> & FEM_kinematicsMap, const std::vector<double> & full_d_vector, 
                std::map<int, int> & nodeIndexMap, const std::vector<std::size_t> & FEM_elementTags)
{
    size_t numNodes;
    size_t tmp_nodeTag;
    int node_index;
    double numerator;
    double denominator;
    double angle1;
    double angle2;
    double norm1;
    double norm2;
    #pragma omp parallel for private(numNodes, tmp_nodeTag, node_index, numerator, denominator, angle1, angle2, norm1, norm2)
    for(std::size_t i=0; i< FEM_elementTags.size(); ++i)
    {
        elementData* element = & FEM_kinematicsMap[FEM_elementTags[i]];
        numNodes = element->nodes.size();
        element->uc = 0;
        element->vc = 0;
        for(size_t l = 0; l < numNodes; l++)
        {
            tmp_nodeTag = element->nodes[l];
            node_index = nodeIndexMap[tmp_nodeTag];
            element->u[l] = full_d_vector[node_index];
            element->v[l] = full_d_vector[node_index + 1];

            element->uc += element->u[l];
            element->vc += element->v[l];
        }
        element->uc = element->uc / numNodes;
        element->vc = element->vc / numNodes;

        numerator = 0;
        denominator = 0;
        for(size_t l = 0; l < numNodes; l++)
        {
            numerator += element->xi[l]*(element->yi[l]+element->v[l]-element->vc);
            numerator -= element->yi[l]*(element->xi[l]+element->u[l]-element->uc);
            denominator += element->xi[l]*(element->xi[l]+element->u[l]-element->uc);
            denominator += element->yi[l]*(element->yi[l]+element->v[l]-element->vc);
        }
        angle1 = atan(numerator/denominator); // returns value between -PI/2 and PI/2
        // we want angle1 and angle2 between -PI and PI
        if(angle1 < 0)
        {
            angle2 = angle1 + PI;
        }
        else
        {
            angle2 = angle1 - PI;
        }
        // filling ui and vi while finding the angle minimizing equation (70) in the report (VERIFIER EQUATION)
        std::vector<double> tmp_ui1(numNodes);
        std::vector<double> tmp_ui2(numNodes);
        std::vector<double> tmp_vi1(numNodes);
        std::vector<double> tmp_vi2(numNodes);
        norm1 = 0;
        norm2 = 0;
        for(size_t l = 0; l < numNodes; l++)
        {
            tmp_ui1[l] = cos(angle1)*(element->xi[l]+element->u[l]-element->uc) + sin(angle1)*(element->yi[l]+element->v[l]-element->vc) - element->xi[l];
            tmp_ui2[l] = cos(angle2)*(element->xi[l]+element->u[l]-element->uc) + sin(angle2)*(element->yi[l]+element->v[l]-element->vc) - element->xi[l];
            tmp_vi1[l] = -sin(angle1)*(element->xi[l]+element->u[l]-element->uc) + cos(angle1)*(element->yi[l]+element->v[l]-element->vc) - element->yi[l];
            tmp_vi2[l] = -sin(angle2)*(element->xi[l]+element->u[l]-element->uc) + cos(angle2)*(element->yi[l]+element->v[l]-element->vc) - element->yi[l];
            norm1 += tmp_ui1[l]*tmp_ui1[l] + tmp_vi1[l]*tmp_vi1[l];
            norm2 += tmp_ui2[l]*tmp_ui2[l] + tmp_vi2[l]*tmp_vi2[l];
        }
        if(norm1 < norm2)
        {
            element->theta = angle1;
            element->ui = tmp_ui1;
            element->vi = tmp_vi1;
        }
        else
        {
            element->theta = angle2;
            element->ui = tmp_ui2;
            element->vi = tmp_vi2;
        }
        for(size_t l = 0; l < numNodes; l++)
        {
            element->pl(2*l) = element->ui[l];
            element->pl(2*l+1) = element->vi[l];
        }
        updateKinMatrices(element, true);
    }
}

// retrieving the global nodal forces based on local nodal forces of the 
// "FEM_elementTags" contained in the "FEM_kinematicsMap", 
// each node Tag corresponds to a global index given by "nodeIndexMap".
void retrieveTotalNodalForces(std::vector<double> & final_nodal_forces, std::map<int, elementData> & FEM_kinematicsMap, 
                        std::map<int, int> & nodeIndexMap, const std::vector<std::size_t> & FEM_elementTags)
{
    for(size_t i = 0; i < final_nodal_forces.size(); i++)
    {
        final_nodal_forces[i] = 0;
    }

    #pragma omp parallel for
    for(size_t i = 0; i < FEM_elementTags.size(); i++)
    {
        elementData* element = & FEM_kinematicsMap[FEM_elementTags[i]];
        int numNodes = (int) element->nodes.size();
        Eigen::Matrix<double, Eigen::Dynamic, 1> fg = element->C.transpose()*element->fl;
        for(int j = 0; j < numNodes; j++)
        {
            int nodeTag = element->nodes[j];
            #pragma omp atomic update
            final_nodal_forces[nodeIndexMap[nodeTag]] += fg(2*j);
            #pragma omp atomic update
            final_nodal_forces[nodeIndexMap[nodeTag] + 1] += fg(2*j + 1);
        }
    }
}

// assembles the "global_f_vector" based on the local nodal forces vector of the 
// "FEM_elementTags" contained in the "FEM_kinematicsMap", 
// without taking into account the DOFs fixed by a dirichlet boundary condition,
// the nodes are renumbered using "new_DOFindices".
// each node Tag corresponds to a global index given by "nodeIndexMap".
void assembleGlobalF(std::vector<double> & global_f_vector, std::map<int, elementData> & FEM_kinematicsMap, 
                std::map<int, int> & nodeIndexMap, const std::vector<int> & new_DOFindices, 
                const std::vector<std::size_t> & FEM_elementTags)
{
    #pragma omp parallel for
    for(size_t i = 0; i < FEM_elementTags.size(); i++)
    {
        elementData* element = & FEM_kinematicsMap[FEM_elementTags[i]];
        int numNodes = (int) element->nodes.size();
        Eigen::Matrix<double, Eigen::Dynamic, 1> fg = element->C.transpose()*element->fl;
        for(int j = 0; j < numNodes; j++)
        {
            int nodeTag = element->nodes[j];
            int index1 = new_DOFindices[nodeIndexMap[nodeTag]]; // current index related to x DOF of the node (= -1 if dirichlet BC)
            int index2 = new_DOFindices[nodeIndexMap[nodeTag] + 1]; // current index related to y DOF of the node (= -1 if dirichlet BC)
            if(index1 >= 0) // need to assemble in global f vector
            {
                #pragma omp atomic update
                global_f_vector[index1] += fg(2*j);
            }
            if(index2 >= 0)
            {
                #pragma omp atomic update
                global_f_vector[index2] += fg(2*j + 1);
            }
        }
    }
}

// assembling the global tangent stiffness matrix based on the local elemental stiffness matrices contained
// in the "FEM_kinematicsMap" associated to each "FEM_elementTags".
// "new_DOFindices" maps the global initial node index (from "nodeIndexMap") to the index
// obtained without considering the fixed DOFs.
void assembleGlobalKg(Eigen::SparseMatrix<double> & global_tangent_matrix, std::map<int, elementData> & FEM_kinematicsMap, 
                    std::map<int, int> & nodeIndexMap, const std::vector<int> & new_DOFindices, 
                    const std::vector<std::size_t> & FEM_elementTags)
{
    // to allow efficient parallelization, one triplet list is created for every thread, then the different lists are assembled
    // to further increase efficiency, space is reserved for each list
    int estimated_number_of_entries = (int) new_DOFindices.size()/2 *8*8; // accurate in the case of Q4 elements, otherwise it provides a good estimation
    std::vector<Eigen::Triplet <double>> Kg_tripletList;
    Kg_tripletList.reserve(estimated_number_of_entries);
    
    int max_threads = omp_get_max_threads();
    int estimated_number_of_entries_per_thread = (int) estimated_number_of_entries / max_threads;
    std::vector<std::vector<Eigen::Triplet <double>>> thread_Kg_list(max_threads);
    for(int i = 0; i < max_threads; i++)
    {
        thread_Kg_list[i].reserve(estimated_number_of_entries_per_thread);
    }

    std::vector<Eigen::Triplet <double>>* my_Kg_list;

    #pragma omp parallel for private(my_Kg_list)
    for(size_t i = 0; i < FEM_elementTags.size(); i++)
    {
        elementData* element = & FEM_kinematicsMap[FEM_elementTags[i]];
        int numNodes = (int) element->nodes.size();
        Eigen::Matrix<double, Eigen::Dynamic, 1> tmp_F = element->P.transpose()*element->fl;

        Eigen::Matrix<double, 1, Eigen::Dynamic> F(1,2*numNodes);
        for(int j = 0; j < numNodes; j++)
        {
            F(0,2*j) = tmp_F(2*j+1, 0);
            F(0,2*j+1) = -tmp_F(2*j, 0);
        }
        
        // see report for the formula
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Kh = element->E*(-F.transpose()*element->G - element->G.transpose()*F*element->P)*element->E.transpose();

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Kg = element->C.transpose()*element->Kl*element->C + Kh;

        int threadNum = omp_get_thread_num();
        my_Kg_list = & thread_Kg_list[threadNum];

        // assemble Kg into tripletlist
        for (int j = 0; j < 2*numNodes; ++j)
        {
            int first_node_index = (int) j/2;
            int first_node = element->nodes[first_node_index];
            int my_row = new_DOFindices[nodeIndexMap[first_node] + j % 2]; 
            if(my_row >= 0) // else, the DOF corresponds to a Dirichlet BC
            {
                for(int l = 0; l <= j; ++l)
                {   // to exploit symmetry
                    int second_node_index = (int) l/2;
                    int second_node = element->nodes[second_node_index];
                    int my_col = new_DOFindices[nodeIndexMap[second_node] + l % 2];

                    if(my_col >= 0) // else, the DOF corresponds to a Dirichlet BC
                    {
                        (*my_Kg_list).push_back(Eigen::Triplet<double>(my_row,my_col,Kg(j,l)));
                        if(j != l)
                        {
                            (*my_Kg_list).push_back(Eigen::Triplet<double>(my_col,my_row,Kg(j,l)));
                        }
                    }
                }
            }
        }
    }
    // gathering the different lists in one single list
    Kg_tripletList = thread_Kg_list[0];
    if(max_threads > 1)
    {
        for (int i = 1; i < max_threads; i++)
        {
            Kg_tripletList.insert(Kg_tripletList.end(), thread_Kg_list[i].begin(), thread_Kg_list[i].end());
        }
    }
    global_tangent_matrix.setFromTriplets(Kg_tripletList.begin(), Kg_tripletList.end());
}

// updating the "full_d_vector" based on the "global_p_vector", using the global and reduced index mapping:
// "new_DOFindices".
void updateNodalDispVector(std::vector<double> & full_d_vector, const std::vector<double> & global_p_vector, 
                        const std::vector<int> & new_DOFindices)
{
    for(size_t i = 0; i < full_d_vector.size(); i++)
    {
        if(new_DOFindices[i] >= 0)
        {
            full_d_vector[i] = global_p_vector[new_DOFindices[i]];
        }
    }
}

// retrieving the nodal displacements on the BEM_FEM_boundary in order to communicate them to the
// BEM solver in the two-way coupled solver. Using the "groups" to find the physical group associated
// to the boundary, all nodes of the boundary are associated through the "boundaryDisplacementMap" their
// nodal displacement (u,v) retrieved from "full_d_vector" using the indexation of "nodeIndexMap".
void retrieveBoundaryDisplacements(std::map<int,std::pair<double,double>> & boundaryDisplacementMap, 
                            const std::vector<double> & full_d_vector, std::map<int, int> & nodeIndexMap, 
                            std::map<std::string, std::pair<int, int>> & groups)
{
    std::string groupname = "BEM_FEM_boundary";
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;

    int nodeTag;
    double u;
    double v;

    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;

    for(std::size_t i=0; i< tags.size(); ++i)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(std::size_t j=0; j< elementTypes.size(); ++j)
        {

            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[j],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);

            for(std::size_t k=0; k < elementTags[j].size(); ++k)
            {
                for(int l = 0; l < numNodes; l++)
                {
                    nodeTag = elnodeTags[j][numNodes*k + l];
                    u = full_d_vector[nodeIndexMap[nodeTag]];
                    v = full_d_vector[nodeIndexMap[nodeTag] + 1];
                    boundaryDisplacementMap[nodeTag] = std::pair<double, double>(u, v);
                }
            }
        }
    }
}

// computing the reaction forces on all the physical groups "groups" for which a dirichlet boundary
// condition has been prescribed (using the "BCkeys" to retrieve the groups from the .geo file),
// based on the "final_nodal_forces" and the "nodeIndexMap" to retrieve the global index associated
// to one given node tag.
void computeReactionForces(std::map<std::string, std::pair<int, int>> & groups,
                     std::map<int, int> & nodeIndexMap, const std::vector<double> & final_nodal_forces,
                     const std::vector<std::string> & BCkeys)
{
    double Rx;
    double Ry;

    int dim;
    int tag;
    std::vector<int> tags;

    std::vector<std::size_t> nodeTags;
    std::vector<double> coord; //useless but mandatory
    std::vector<double> parametricCoord; //useless but mandatory

    for (auto &key : BCkeys)
    {
        // get corresponding value
        std::vector<double> value;
        gmsh::onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("ux") == 0 || words[2].compare("uy") == 0){
                Rx = 0;
                Ry = 0;
                std::string groupname = words[1];
                std::cout << "-------------------- " << groupname << " ------------------- \n";
                
                //loop over elements of group_name
                dim = groups[groupname].first;
                tag = groups[groupname].second;
                if (tag == 0)
                {
                    std::cerr << "Group '" << groupname << "' does not exist!\n";
                    return;
                }
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

                // loop over entities
                for (std::size_t k = 0; k < tags.size(); ++k)
                {
                    // retrieving all the nodes in this particular entity
                    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tags[k], true);
                    
                    
                    for (std::size_t j = 0; j < nodeTags.size(); ++j){
                        if(words[2].compare("ux") == 0){    // computing Rx
                            int my_new_col = nodeIndexMap[nodeTags[j]]; //nodeIndexMap[nodeTags[j]] corresponds to line/column index in K, f
                            Rx += final_nodal_forces[my_new_col];
                        }
                        if(words[2].compare("uy") == 0){    // computing Ry
                            int my_new_col = nodeIndexMap[nodeTags[j]] + 1;
                            Ry += final_nodal_forces[my_new_col];
                        }
                    }
                }
                std::cout << "Rx: " << Rx << " [N/m], Ry: " << Ry << " [N/m] \n";
            }
        }
    }
}

// plotting the nodal displacement values associated to "FEM_nodeTags" in the "names[0]" model,
// based on the "full_d_vector".
void plotUxUy(const std::vector<std::size_t> & FEM_nodeTags, const std::vector<double> & full_d_vector, 
            const std::vector<std::string> & names)
{
    int viewtagUx = gmsh::view::add("ux");
    int viewtagUy = gmsh::view::add("uy");
    std::vector<std::vector<double>> dataUx(FEM_nodeTags.size());
    std::vector<std::vector<double>> dataUy(FEM_nodeTags.size());
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        dataUx[i].resize(1);
        dataUx[i][0] = full_d_vector[2*i]; //ux // à modifier ici
        dataUy[i].resize(1);
        dataUy[i][0] = full_d_vector[2*i+1]; //uy
    }

    gmsh::view::addModelData(viewtagUx, 0, names[0], "NodeData",
                            FEM_nodeTags, dataUx); //independent of time here
    gmsh::view::addModelData(viewtagUy, 0, names[0], "NodeData",
                            FEM_nodeTags, dataUy); //independent of time here
}

// computing the elemental-nodal strains and stresses in the FEM domain (dimension "dim" and entity tags "tags"),
// based on Hooke's matrix in plane stress "H_matrix", Young's modulus "E", Poisson's ratio "nu", looping over the
// "nbelems" elements contained in "FEM_kinematicsMap".
// returning the "max_VM_stress" and the corresponding nodetag "max_VM_nodeTag".
void computeStressStrainElNodal(const int & dim, const std::vector<int> & tags, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           const double nu, const double E, std::map<int, elementData> & FEM_kinematicsMap,
                           const int & nbelems, double & max_VM_stress, int & max_VM_nodeTag)
{
    std::vector<std::string> names;
    gmsh::model::list(names);

    int viewTagEpsX = gmsh::view::add("nodal_eps_x");
    int viewTagEpsY = gmsh::view::add("nodal_eps_y");
    int viewTagEpsXY = gmsh::view::add("nodal_2eps_xy");
    int viewTagEpsZ = gmsh::view::add("nodal_eps_z");
    int viewTagSigX = gmsh::view::add("nodal_sig_x");
    int viewTagSigY = gmsh::view::add("nodal_sig_y");
    int viewTagSigXY = gmsh::view::add("nodal_sig_xy");
    int viewTagSigVM = gmsh::view::add("nodal_sig_VM");

    std::vector<std::vector<std::vector<double>>> dataEps(4);
    dataEps[0].resize(nbelems); // eps X
    dataEps[1].resize(nbelems); // eps Y
    dataEps[2].resize(nbelems); // eps XY
    dataEps[3].resize(nbelems); // eps Z

    std::vector<std::vector<std::vector<double>>> dataSig(4);
    dataSig[0].resize(nbelems); // Sig X
    dataSig[1].resize(nbelems); // Sig Y
    dataSig[2].resize(nbelems); // Sig XY
    dataSig[3].resize(nbelems); // Sig VM
    std::vector<std::size_t> elems2D(nbelems);

    int data_index = 0;

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;
    for(std::size_t entityTag = 0; entityTag < tags.size(); ++entityTag)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[entityTag]);

        std::string name;
        int dim, order, numNodes, numPrimaryNodes;
        std::vector<double> paramCoord2D;
        int numComponents, numOrientations;
        std::vector<double> basisFunctionsGradient;

        Eigen::Vector3d epsilon; // will contain local values of the strains (Voigt notation)
        Eigen::Vector3d sigma; // will contain local values of the stresses (Voigt notation)

        Eigen::Matrix<double, 2, 2> jacob2D;
        Eigen::Matrix<double, 2, 2> jacobinv;
        Eigen::Matrix<double, 2, 2> jacobinvtrans;

        elementData* element;
        double theta;

        Eigen::Matrix<double, 2, 2> strain_tensor;
        Eigen::Matrix<double, 2, 2> stress_tensor;
        Eigen::Matrix<double, 2, 2> rotation_tensor;

        for (std::size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
        {
            // get info about this type of element
            gmsh::model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                                    numNodes, paramCoord2D, numPrimaryNodes);
            //paramCoord contains the coordinates of the nodes in 2D
            // passage en 3D pour utiliser getJacobians
            std::vector<double> paramCoord(3*paramCoord2D.size()/2);
            for(std::size_t j = 0; j < paramCoord2D.size()/2; ++j){
                paramCoord[3*j] = paramCoord2D[2*j];
                paramCoord[3*j+1] = paramCoord2D[2*j+1];
                paramCoord[3*j+2] = 0.;
            }
            
            // select element/node tags for this type of element
            auto &etags = elementTags[ityp];

            // get determinants and Jacobians (only jacob useful here) AT NODES AND NOT AT GAUSS POINTS
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[ityp], 
                                            paramCoord, jacobians, 
                                            determinants, coords, tags[entityTag]);

            // derivatives of the shape functions at the nodes of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            gmsh::model::mesh::getBasisFunctions(elementTypes[ityp], paramCoord, 
                                                    "GradLagrange", numComponents, 
                                                    basisFunctionsGradient, 
                                                    numOrientations);

            // loop over elements of type "ityp"
            for (std::size_t i = 0; i < etags.size(); ++i)
            {
                element = & FEM_kinematicsMap[etags[i]];
                elems2D[data_index] = etags[i];

                theta = element->theta;
                rotation_tensor(0,0) = cos(theta);
                rotation_tensor(0,1) = sin(theta);
                rotation_tensor(1,0) = -sin(theta);
                rotation_tensor(1,1) = cos(theta);

                dataEps[0][data_index].resize(numNodes);
                dataEps[1][data_index].resize(numNodes);
                dataEps[2][data_index].resize(numNodes);
                dataEps[3][data_index].resize(numNodes);
                dataSig[0][data_index].resize(numNodes);
                dataSig[1][data_index].resize(numNodes);
                dataSig[2][data_index].resize(numNodes);
                dataSig[3][data_index].resize(numNodes);

                for (int j = 0; j < numNodes; ++j){
                    // converting jacobian to Eigen format
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*j+9*numNodes*i], 3, 3);
                    // 3D to 2D (plane stress case)
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix (formula slide 19) 
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                    computeBmatrix(j, numNodes, basisFunctionsGradient, determinants, jacobinvtrans, B_matrix);
                    
                    // expressed in local axes (and in Voigt notation)
                    epsilon = B_matrix*element->pl; 
                    sigma = H_matrix*epsilon;

                    // converting to tensor notation (still in local axes)
                    strain_tensor(0,0) = epsilon(0);
                    strain_tensor(1,0) = epsilon(2)/2; // divide by 2 because engineering shear strain -> pure shear strain
                    strain_tensor(0,1) = epsilon(2)/2;
                    strain_tensor(1,1) = epsilon(1);

                    stress_tensor(0,0) = sigma(0);
                    stress_tensor(1,0) = sigma(2);
                    stress_tensor(0,1) = sigma(2);
                    stress_tensor(1,1) = sigma(1);

                    strain_tensor = rotation_tensor.transpose()*strain_tensor*rotation_tensor;
                    stress_tensor = rotation_tensor.transpose()*stress_tensor*rotation_tensor;

                    dataEps[0][data_index][j] = strain_tensor(0,0);
                    dataEps[1][data_index][j] = strain_tensor(1,1);
                    dataEps[2][data_index][j] = strain_tensor(1,0);
                    dataEps[3][data_index][j] = -nu/E*(stress_tensor(0,0)+stress_tensor(1,1));// from plane stress hypothesis
                    dataSig[0][data_index][j] = stress_tensor(0,0);
                    dataSig[1][data_index][j] = stress_tensor(1,1);
                    dataSig[2][data_index][j] = stress_tensor(1,0);
                    // VM stress in plane stress
                    dataSig[3][data_index][j] = sqrt(stress_tensor(0,0)*stress_tensor(0,0) + stress_tensor(1,1)*stress_tensor(1,1) - stress_tensor(0,0)*stress_tensor(1,1) + 3*stress_tensor(1,0)*stress_tensor(1,0));
                    
                    //retrieving maximal nodal von mises stress
                    if(dataSig[3][data_index][j] > max_VM_stress)
                    {
                        max_VM_stress = dataSig[3][data_index][j];
                        max_VM_nodeTag = (int) elnodeTags[ityp][numNodes*i+j];
                    }
                }
                data_index++;
            }
        }
    }

    gmsh::view::addModelData(viewTagEpsX, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[0], 1); 
    gmsh::view::addModelData(viewTagEpsY, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[1], 1); 
    gmsh::view::addModelData(viewTagEpsXY, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[2], 1); 
    gmsh::view::addModelData(viewTagEpsZ, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[3], 1); 
    gmsh::view::addModelData(viewTagSigX, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[0], 1); 
    gmsh::view::addModelData(viewTagSigY, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[1], 1); 
    gmsh::view::addModelData(viewTagSigXY, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[2], 1); 
    gmsh::view::addModelData(viewTagSigVM, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[3], 1); 
}

// computing the work done by external forces, as the nodal forces are energy consistent, it
// is sufficient to compute the "externalWork" as the dot product of the nodal displacements
// "full_d_vector" and the nodal forces "final_nodal_forces".
void computeWorkDoneByExternalForces(double &externalWork, const std::vector<double> & full_d_vector, 
                                const std::vector<double> & final_nodal_forces)
{
    for(size_t i = 0; i < full_d_vector.size(); i++)
    {
        externalWork += full_d_vector[i]*final_nodal_forces[i];
    }
}

// computing the total (internal) strain energy "strain_energy" in the whole FEM domain (dimension "dim" and 
// entity tags "tags"), using the elementData stored in "FEM_kinematicsMap" gathering the informations for all
// FEM elements. The Hooke's matrix in plane stress "H_matrix" is used to compute the stresses from the strains.
// Computed as a sum of integrals over the finite elements, using Gauss integration "integration_rule".
void computeStrainEnergy(double & strainEnergy, const int & dim, const std::vector<int> & tags,
                           std::map<int, elementData> & FEM_kinematicsMap, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           const std::string & integration_rule)
{
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;
    std::vector<double> basisFunctionsGradient;

    for (std::size_t k = 0; k < tags.size(); ++k)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types
        for (std::size_t i = 0; i < elementTypes.size(); ++i)
        {
            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, elDim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);
            
            // get Gauss points coordinates and weights
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    integration_rule,
                                                    localCoords, weights);

            // get determinants and Jacobians
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords, tags[k]);

            // values and derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);

            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctionsGradient, 
                                                 numOrientations);

            #pragma omp parallel for
            for (std::size_t el = 0; el < elementTags[i].size(); el++){

                Eigen::Matrix<double, Eigen::Dynamic, 1> local_nodal_displacement = FEM_kinematicsMap[elementTags[i][el]].pl;

                double elementalStrainEnergy = 0;

                // looping over the Gauss points: Gauss integration
                for (std::size_t gaussIndex = 0; gaussIndex < weights.size(); ++gaussIndex){
                    // converting jacobian to Eigen format
                    Eigen::Matrix<double, 2, 2> jacob2D;
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*gaussIndex+9*weights.size()*el], 3, 3);
                    // 3D to 2D (plane stress case)
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    Eigen::Matrix<double, 2, 2> jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    Eigen::Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                    computeBmatrix(gaussIndex, numNodes, basisFunctionsGradient, determinants, jacobinvtrans, B_matrix);

                    Eigen::Vector3d epsilon = B_matrix*local_nodal_displacement;
                    Eigen::Vector3d sigma = H_matrix*epsilon;

                    double dot_product = 0;
                    for(int z = 0; z < 3; z++)
                    {
                        dot_product += epsilon(z)*sigma(z);
                    }

                    elementalStrainEnergy += 0.5*dot_product*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                #pragma omp atomic update
                strainEnergy += elementalStrainEnergy;
            }
        }
    }
}

//displays the time of the desired code section "function_name", if "display_time" is set as true.
void displayTime(const double &start, const double &end, const std::string &function_name, const bool &display_time)
{
    if(display_time)
    {
        std::cout << std::endl;
        std::cout << function_name << std::endl;
        std::cout << "Elapsed time in nanoseconds: "
        << end - start << " ns" << std::endl;

        std::cout << "Elapsed time in microseconds: "
        << end - start << " µs" << std::endl;

        std::cout << "Elapsed time in milliseconds: "
        << end - start << " ms" << std::endl;

        std::cout << "Elapsed time in seconds: "
        << end - start << " sec";
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

/*-------------- FUNCTIONS SPECIFIC TO LINEAR FEM ------------------*/

// assembling the full stiffness matrix "full_K_matrix" using first a triplet vector, assembling the volumic part 
// of the nodal forces "full_f_vector", in the FEM_domain (dimension "dim" and entity tag "tags") using the global
// node indexation defined by "nodeIndexMap". "H_matrix" is the plane stress Hooke's matrix, "rho" the density of the
// material, "bx" the volumic force along x, "by" the volumic force along y.
void assembleKVolumicF(const Eigen::Matrix<double, 3, 3> & H_matrix, const double & rho, const double & bx, const double & by,
                const int & dim, const std::vector<int> & tags, std::map<int, int> & nodeIndexMap, 
                Eigen::SparseMatrix<double> & full_K_matrix, std::vector<double> & full_f_vector)
{
    // to allow efficient parallelization, one triplet list (as a vector) is created for every thread, then the different lists are assembled
    // to further increase efficiency, space is reserved for each list
    int estimated_number_of_entries = (int) full_f_vector.size()/2 *8*8; // accurate in the case of Q4 elements, otherwise it provides a good estimation
    std::vector<Eigen::Triplet <double>> K_tripletList;
    K_tripletList.reserve(estimated_number_of_entries);
    
    int max_threads = omp_get_max_threads();
    int estimated_number_of_entries_per_thread = (int) estimated_number_of_entries / max_threads;
    std::vector<std::vector<Eigen::Triplet <double>>> thread_K_list(max_threads);
    for(int i = 0; i < max_threads; i++)
    {
        thread_K_list[i].reserve(estimated_number_of_entries_per_thread);
    }

    std::vector<Eigen::Triplet <double>>* my_K_list;

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;
    std::vector<double> basisFunctionsGradient;

    for (std::size_t k = 0; k < tags.size(); ++k)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types
        for (std::size_t i = 0; i < elementTypes.size(); ++i)
        {
            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, elDim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);

            // adapting the composite integration rule as a function of the order of the elements
            std::string adapting_integration_rule = "CompositeGauss14";
            if(order < 8)
            {
                adapting_integration_rule = "CompositeGauss" + std::to_string(2*order);
            }
            
            // get Gauss points coordinates and weights
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    adapting_integration_rule, 
                                                    localCoords, weights);

            // get determinants and Jacobians
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords, tags[k]);

            // values and derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);

            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctionsGradient, 
                                                 numOrientations);

            /*---------computing the K matrix (and f vector, only volumic part for now)------*/
            // looping over the elements (computing elemental K matrices)
            #ifdef _MSC_VER
            // https://github.com/tesseract-ocr/tesseract/issues/3285
            // OpenMP 2.0 standard is enforced - we may use /openmp:llvm and keep size_t
            #pragma omp parallel for
            for (int el = 0; el < elementTags[i].size(); el++){
            #else
            #pragma omp parallel for private (my_K_list)
            for (std::size_t el = 0; el < elementTags[i].size(); el++)
            {
            #endif
                // elemental K matrix and f vector initialization
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K_matrix(2*numNodes,2*numNodes);
                std::vector<double> f_vector(2*numNodes);
                K_matrix.setZero();
                for (int j = 0; j < 2*numNodes; ++j)f_vector[j] = 0.;

                int threadNum = omp_get_thread_num();
                my_K_list = & thread_K_list[threadNum];

                // looping over the Gauss points
                for (std::size_t gaussIndex = 0; gaussIndex < weights.size(); ++gaussIndex){
                    // converting jacobian to Eigen format
                        Eigen::Matrix<double, 2, 2> jacob2D;
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*gaussIndex+9*weights.size()*el], 3, 3);
                    // 3D to 2D (plane stress case)
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    Eigen::Matrix<double, 2, 2> jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    Eigen::Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix and the f vector (also using gauss integration)
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                    computeBmatrix(gaussIndex, numNodes, basisFunctionsGradient, determinants, jacobinvtrans, B_matrix);
                    
                    for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
                        f_vector[2*l] = f_vector[2*l] + rho*basisFunctions[l+numNodes*gaussIndex]*bx*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                        f_vector[2*l+1] = f_vector[2*l+1] + rho*basisFunctions[l+numNodes*gaussIndex]*by*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex]; 
                    }

                    K_matrix = K_matrix + B_matrix.transpose()*H_matrix*B_matrix*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                // USE LIST OF TRIPLETS TO PUSH ONLY ONCE -> more rapid and easier to parallelize
                assembleKtriplet(el, nodeTags[i], numNodes, nodeIndexMap, *my_K_list, K_matrix);
                assembleFlocal(el, nodeTags[i], numNodes, nodeIndexMap, full_f_vector, f_vector);
            }
        }
    }
    // gathering all triplet vectors in one single vector
    K_tripletList = thread_K_list[0];
    if(max_threads > 1)
    {
        for (int i = 1; i < max_threads; i++)
        {
            K_tripletList.insert(K_tripletList.end(), thread_K_list[i].begin(), thread_K_list[i].end());
        }
    }
    // assembling the full sparse K matrix
    full_K_matrix.setFromTriplets(K_tripletList.begin(), K_tripletList.end());
}

// assembles the elemental "K_matrix" into a vector of triplets "k_tripletList" using 
// the node tags provided by "nodeTags" vector (index of the element is "elIndex").
// The element has "numNodes" nodes. "nodeIndexMap" provides the global nodal index corresponding to a given node tag.
void assembleKtriplet(const int &elIndex, const std::vector<std::size_t> & nodeTags, const int numNodes, 
                    std::map<int, int> & nodeIndexMap, std::vector<Eigen::Triplet <double>> & k_tripletList,
                    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & K_matrix)
{
    //we go through the lower half of the elemental K_matrix before assembling in the triplet vector
    for (int j = 0; j < 2*numNodes; ++j){
        int first_node_index = (int) j/2;
        int first_node = nodeTags[numNodes*elIndex+first_node_index];
        int my_row = nodeIndexMap[first_node] + j % 2; 
        for(int l = 0; l <= j; ++l){ // to exploit symmetry
            int second_node_index = (int) l/2;
            int second_node = nodeTags[numNodes*elIndex+second_node_index];
            int my_col = nodeIndexMap[second_node] + l % 2;

            k_tripletList.push_back(Eigen::Triplet<double>(my_row,my_col,K_matrix(j,l)));
            if(j != l){
                k_tripletList.push_back(Eigen::Triplet<double>(my_col,my_row,K_matrix(j,l)));
            }
        }
    }
}

// the "full_f_vector" and "full_K_matrix" are reduced in "reduced_f_vector" and "reduced_K_matrix" respectively,
// taking into account the dirichlet boundary conditions contained in "dir_BC".
// "new_DOFindices" gathers the new indices of the DOFs when taking the boundary conditions into account.
void fillReducedKf(const std::vector<double> & full_f_vector, const Eigen::SparseMatrix<double> & full_K_matrix,
                const std::vector<int> & new_DOFindices, const std::vector<double> & dir_BC,
                Eigen::SparseMatrix<double> & reduced_K_matrix, std::vector<double> & reduced_f_vector)
{
    //the correction term on the f vector comes from the "Finite Element Method" theoretical course (JP Ponthot, ULiege), Chap2 Slide 38/50
    std::list<Eigen::Triplet <double>> reduced_k_tripletList;

    // initializing the reduced f vector
    for (std::size_t j = 0; j < reduced_f_vector.size(); ++j){
        reduced_f_vector[j] = 0.;
    }

    // looping over the non-zero elements of full_K_matrix
    int row;
    int new_row;
    int new_col;
    double k_value;
    for (int col=0; col < full_K_matrix.outerSize(); ++col){
        new_col = new_DOFindices[col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(full_K_matrix,col); it; ++it){
            k_value = it.value();
            row = it.row();   // row index
            new_row = new_DOFindices[row];
            if(new_row >= 0 && new_col >= 0){ // if the given element must not be removed
                if(new_row < new_col){ // taking advantage of the symmetry
                    reduced_k_tripletList.push_back(Eigen::Triplet<double>(new_row,new_col,k_value));
                    reduced_k_tripletList.push_back(Eigen::Triplet<double>(new_col,new_row,k_value));
                }
                if(new_row == new_col){ // diagonal element
                    reduced_k_tripletList.push_back(Eigen::Triplet<double>(new_row,new_col,k_value));
                }
            }
            if(new_col < 0 && new_row >= 0 && dir_BC[col]){ // correction term RHS due to non_homogeneous boundary condition
                reduced_f_vector[new_row] -= k_value*dir_BC[col];
            }
        }
        if(new_col >= 0 ){ // line must not be removed
            reduced_f_vector[new_col] += full_f_vector[col];
        }
    }
    reduced_K_matrix.setFromTriplets(reduced_k_tripletList.begin(), reduced_k_tripletList.end());
}


// computing the elemental-nodal strains and stresses in the FEM domain (dimension "dim" and entity tags "tags"),
// based on Hooke's matrix in plane stress "H_matrix", Young's modulus "E", Poisson's ratio "nu", and the full nodal
// displacement vector "full_d_vector" and the global nodal index map "nodeIndexMap".
// "nbelems" is the number of 2D elements contained in the FEM domain.
// returning the "max_VM_stress" and the corresponding nodetag "max_VM_nodeTag".
void computeStressStrainLinearWay(const int & dim, const std::vector<int> & tags, const std::vector<double> &full_d_vector,
                           const Eigen::Matrix<double, 3, 3> &H_matrix, const double &nu, const double &E,
                           std::map<int, int> & nodeIndexMap, const int &nbelems, double &max_VM_stress, int &max_VM_nodeTag)
{
    std::vector<std::string> names;
    gmsh::model::list(names);

    int viewTagEpsX = gmsh::view::add("nodal_eps_x");
    int viewTagEpsY = gmsh::view::add("nodal_eps_y");
    int viewTagEpsXY = gmsh::view::add("nodal_2eps_xy");
    int viewTagEpsZ = gmsh::view::add("nodal_eps_z");
    int viewTagSigX = gmsh::view::add("nodal_sig_x");
    int viewTagSigY = gmsh::view::add("nodal_sig_y");
    int viewTagSigXY = gmsh::view::add("nodal_sig_xy");
    int viewTagSigVM = gmsh::view::add("nodal_sig_VM");

    std::vector<std::vector<std::vector<double>>> dataEps(4);
    dataEps[0].resize(nbelems); // eps X
    dataEps[1].resize(nbelems); // eps Y
    dataEps[2].resize(nbelems); // eps XY
    dataEps[3].resize(nbelems); // eps Z

    std::vector<std::vector<std::vector<double>>> dataSig(4);
    dataSig[0].resize(nbelems); // Sig X
    dataSig[1].resize(nbelems); // Sig Y
    dataSig[2].resize(nbelems); // Sig XY
    dataSig[3].resize(nbelems); // Sig VM
    std::vector<std::size_t> elems2D(nbelems);


    int data_index = 0;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;

    for(std::size_t entityTag = 0; entityTag< tags.size(); ++entityTag)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[entityTag]);

        std::string name;
        int dim, order, numNodes, numPrimaryNodes;
        std::vector<double> paramCoord2D;
        int numComponents, numOrientations;
        std::vector<double> basisFunctionsGradient;
        Eigen::Vector3d epsilon; // will contain local values of the strains
        Eigen::Vector3d sigma; // will contain local values of the stresses
        int nodeTag;
        Eigen::Matrix<double, 2, 2> jacob2D;
        Eigen::Matrix<double, 2, 2> jacobinv;
        Eigen::Matrix<double, 2, 2> jacobinvtrans;

        for (std::size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
        {
            // get info about this type of element
            gmsh::model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                                    numNodes, paramCoord2D, numPrimaryNodes);
            //paramCoord contains the coordinates of the nodes in 2D
            // passage en 3D pour utiliser getJacobians
            std::vector<double> paramCoord(3*paramCoord2D.size()/2);
            for(std::size_t j = 0; j < paramCoord2D.size()/2; ++j){
                paramCoord[3*j] = paramCoord2D[2*j];
                paramCoord[3*j+1] = paramCoord2D[2*j+1];
                paramCoord[3*j+2] = 0.;
            }
            

            // select element/node tags for this type of element
            auto &etags = elementTags[ityp];
            auto &enods = elnodeTags[ityp];

            // will contain the d vector inside each element
            Eigen::VectorXd nodal_displacement(2*numNodes); //use eigen vector to use Eigen matric product later

            // get determinants and Jacobians (only jacob useful here) AT NODES AND NOT AT GAUSS POINTS
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[ityp], 
                                            paramCoord, jacobians, 
                                            determinants, coords, tags[entityTag]);

            // derivatives of the shape functions at the nodes of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            gmsh::model::mesh::getBasisFunctions(elementTypes[ityp], paramCoord, 
                                                    "GradLagrange", numComponents, 
                                                    basisFunctionsGradient, 
                                                    numOrientations);

            // loop over elements of type "ityp"
            for (std::size_t i = 0; i < etags.size(); ++i)
            {
                elems2D[data_index] = etags[i];

                dataEps[0][data_index].resize(numNodes);
                dataEps[1][data_index].resize(numNodes);
                dataEps[2][data_index].resize(numNodes);
                dataEps[3][data_index].resize(numNodes);
                dataSig[0][data_index].resize(numNodes);
                dataSig[1][data_index].resize(numNodes);
                dataSig[2][data_index].resize(numNodes);
                dataSig[3][data_index].resize(numNodes);

                // loop over nodes of the element, filling the local nodal displacement vector "d"
                for (int j = 0; j < numNodes; ++j)
                {
                    nodeTag = enods[numNodes * i + j];
                    nodal_displacement(2*j) = full_d_vector[nodeIndexMap[nodeTag]];
                    nodal_displacement(2*j+1) = full_d_vector[nodeIndexMap[nodeTag]+1];
                }
                // looping over the nodes
                for (int j = 0; j < numNodes; ++j){
                    // converting jacobian to Eigen format
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*j+9*numNodes*i], 3, 3);
                    // 3D to 2D (plane stress case)
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix  
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                    computeBmatrix(j, numNodes, basisFunctionsGradient, determinants, jacobinvtrans, B_matrix);

                    epsilon = B_matrix*nodal_displacement;
                    sigma = H_matrix*epsilon;

                    dataEps[0][data_index][j] = epsilon(0);
                    dataEps[1][data_index][j] = epsilon(1);
                    dataEps[2][data_index][j] = epsilon(2);
                    dataEps[3][data_index][j] = -nu/E*(sigma(0)+sigma(1));
                    dataSig[0][data_index][j] = sigma(0);
                    dataSig[1][data_index][j] = sigma(1);
                    dataSig[2][data_index][j] = sigma(2);
                    dataSig[3][data_index][j] = sqrt(sigma(0)*sigma(0) + sigma(1)*sigma(1) - sigma(0)*sigma(1) + 3*sigma(2)*sigma(2));

                    //retrieving maximal nodal von mises stress
                    if(dataSig[3][data_index][j] > max_VM_stress)
                    {
                        max_VM_stress = dataSig[3][data_index][j];
                        max_VM_nodeTag = (int) elnodeTags[ityp][numNodes*i+j];
                    }
                }
                data_index++;
            }
        }
    }
    
    gmsh::view::addModelData(viewTagEpsX, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[0], 1);
    gmsh::view::addModelData(viewTagEpsY, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[1], 1);
    gmsh::view::addModelData(viewTagEpsXY, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[2], 1);
    gmsh::view::addModelData(viewTagEpsZ, 0, names[0], "ElementNodeData",
                                elems2D, dataEps[3], 1);
    gmsh::view::addModelData(viewTagSigX, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[0], 1);
    gmsh::view::addModelData(viewTagSigY, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[1], 1);
    gmsh::view::addModelData(viewTagSigXY, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[2], 1);
    gmsh::view::addModelData(viewTagSigVM, 0, names[0], "ElementNodeData",
                                elems2D, dataSig[3], 1);
}

// computing the total (internal) strain energy "strain_energy" in the whole FEM domain (dimension "dim" and 
// entity tags "tags"), using the nodal displacements in "full_d_vector". The nodes are indexed using the "nodeIndexMap".
// The Hooke's matrix in plane stress "H_matrix" is used to compute the stresses from the strains.
// Computed as a sum of integrals over the finite elements, using Gauss integration "integration_rule".
void computeStrainEnergyLinearWay(double & strainEnergy, const int & dim, const std::vector<int> & tags,
                           const std::vector<double> &full_d_vector, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           std::map<int, int> & nodeIndexMap, const std::string & integration_rule)
{
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;
    std::vector<double> basisFunctionsGradient;
    int nodeTag;

    Eigen::Vector3d epsilon; // will contain local values of the strains
    Eigen::Vector3d sigma; // will contain local values of the stresses

    for (std::size_t k = 0; k < tags.size(); ++k)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types
        for (std::size_t i = 0; i < elementTypes.size(); ++i)
        {
            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, elDim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);
            
            // get Gauss points coordinates and weights
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    integration_rule,
                                                    localCoords, weights);

            // get determinants and Jacobians
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords, tags[k]);

            // values and derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);

            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctionsGradient, 
                                                 numOrientations);
                                                
            Eigen::VectorXd nodal_displacement(2*numNodes); //use eigen vector to use Eigen matric product later

            for (std::size_t el = 0; el < elementTags[i].size(); el++){

                for (int j = 0; j < numNodes; ++j)
                {
                    nodeTag = nodeTags[i][numNodes * el + j];
                    nodal_displacement(2*j) = full_d_vector[nodeIndexMap[nodeTag]];
                    nodal_displacement(2*j+1) = full_d_vector[nodeIndexMap[nodeTag]+1];
                }

                double elementalStrainEnergy = 0;

                // looping over the Gauss points
                for (std::size_t gaussIndex = 0; gaussIndex < weights.size(); ++gaussIndex){
                    // converting jacobian to Eigen format
                    Eigen::Matrix<double, 2, 2> jacob2D;
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*gaussIndex+9*weights.size()*el], 3, 3);
                    // 3D to 2D (plane stress case)
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    Eigen::Matrix<double, 2, 2> jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    Eigen::Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix 
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                    computeBmatrix(gaussIndex, numNodes, basisFunctionsGradient, determinants, jacobinvtrans, B_matrix);

                    epsilon = B_matrix*nodal_displacement;
                    sigma = H_matrix*epsilon;

                    double dot_product = 0;
                    for(int z = 0; z < 3; z++)
                    {
                        dot_product += epsilon(z)*sigma(z);
                    }

                    elementalStrainEnergy += 0.5*dot_product*determinants[gaussIndex+weights.size()*el]*weights[gaussIndex];
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                strainEnergy += elementalStrainEnergy;
            }
        }
    }
}