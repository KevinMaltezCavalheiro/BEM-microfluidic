#ifdef _MSC_VER
#include <gmsh.h_cwrap>       // gmsh main header
#else
#include <gmsh.h>       // gmsh main header
#endif

#include "functionsFEM.hpp"

#include <iostream>     // for std::cout
#include <map>
#include <sstream>  // for std::stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices
#include <Eigen/SparseCholesky> // solving sparse linear systems
#include <Eigen/Core>
#include <cmath>
#include <algorithm> // to sort and merge the FEM node vectors
#include <omp.h>

#define PI 3.141592654

//Compute the matrix H
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

//Compute the full matrix "f" of one given thread
void assembleFthread(const int i, const int el, 
                          const std::vector<std::vector<std::size_t>> & nodeTags, 
                          const int numNodes, std::map<int, int> & nodeIndexMap,
                          std::list<std::pair<int,double>> & thread_f_vector, 
                          std::vector<double> & f_vector)
{
    int first_node_index;
    int first_node;
    int my_row;
    for (int j = 0; j < 2*numNodes; ++j)
    {
        first_node_index = (int) j/2;
        first_node = nodeTags[i][numNodes*el+first_node_index];
        my_row = nodeIndexMap[first_node] + j % 2; 
        thread_f_vector.push_back(std::pair<int, double>(my_row, f_vector[j]));
    }
}


void computeBmatrix(const int j, const int numNodes, const std::vector<double> & basisFunctionsGradient,
                      const std::vector<double> & determinants, 
                      const Eigen::Matrix<double, 2, 2> & jacobinvtrans,
                      const Eigen::Matrix<double, 3, 3> & H_matrix,
                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & B_matrix)
{
    for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
        B_matrix(0,2*l) = jacobinvtrans(0,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(0,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
        B_matrix(1,2*l) = 0.;
        B_matrix(2,2*l) = jacobinvtrans(1,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(1,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
        B_matrix(0,2*l+1) = 0.;
        B_matrix(1,2*l+1) = B_matrix(2,2*l);
        B_matrix(2,2*l+1) = B_matrix(0,2*l);
    }
}

//Get physical properties of domain
void getPhysicalProperties(const std::vector<std::string> & keys, double & E, double & nu, double & rho,
                              double & bx, double & by)
{
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

void fillElementalKandF(const Eigen::Matrix<double, 3, 3> H_matrix, const double rho, const double bx, const double by,
              const std::string & integration_rule, 
              int & dim, std::vector<int> & tags, std::map<int, int>  & nodeIndexMap,
              std::vector<std::list<Eigen::Triplet <double>>> & my_k_list, 
              std::vector<std::list<std::pair<int,double>>> & thread_f_vector,
              std::map<int, elementData> & FEM_kinematicsMap)
{
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    std::vector<double> jacobians, determinants, coords;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;
    std::vector<double> basisFunctionsGradient;

    for (std::size_t k = 0; k < tags.size(); ++k)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types, enlever les déclarations des boucles pr aller plus vite (cfr Boman)
        for (std::size_t i = 0; i < elementTypes.size(); ++i)
        {
            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, elDim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);
            
            // get Gauss points coordinates and weights
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    integration_rule, // "Gauss2"
                                                    localCoords, weights);

            // get determinants and Jacobians
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords);

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
            //double test1 = omp_get_wtime();
            //Eigen::setNbThreads(1);
            #pragma omp parallel for
            #ifdef _MSC_VER
            // https://github.com/tesseract-ocr/tesseract/issues/3285
            // OpenMP 2.0 standard is enforced - we may use /openmp:llvm and keep size_t
            for (int el = 0; el < elementTags[i].size(); el++){
            #else
            for (std::size_t el = 0; el < elementTags[i].size(); el++){
            #endif
                // elemental K matrix and f vector initialization
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K_matrix(2*numNodes,2*numNodes);
                std::vector<double> f_vector(2*numNodes);
                K_matrix.setZero();
                for (int j = 0; j < 2*numNodes; ++j)f_vector[j] = 0.;

                // looping over the Gauss points (formula of slide 17)
                for (std::size_t j = 0; j < weights.size(); ++j){
                    // converting jacobian to Eigen format
                    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                        jacob(&jacobians[9*j+9*weights.size()*el], 3, 3);
                    // 3D to 2D (plane stress case)
                    Eigen::Matrix<double, 2, 2> jacob2D;
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


                    computeBmatrix(j, numNodes, basisFunctionsGradient, determinants, 
                                     jacobinvtrans, H_matrix, B_matrix);
                    
                    for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
                        f_vector[2*l] = f_vector[2*l] + rho*basisFunctions[l+numNodes*j]*bx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
                        f_vector[2*l+1] = f_vector[2*l+1] + rho*basisFunctions[l+numNodes*j]*by*determinants[j+weights.size()*el]*weights[j]; // pas oublier le weights[j] non plus !
                    }

                    K_matrix = K_matrix + B_matrix.transpose()*H_matrix*B_matrix*determinants[j+weights.size()*el]*weights[j];
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                // USE LIST OF TRIPLETS TO PUSH ONLY ONCE -> more rapid and easier to parallelize
                FEM_kinematicsMap[elementTags[i][el]].Kl = K_matrix;
                assembleFthread(i, el, nodeTags, numNodes, nodeIndexMap, thread_f_vector[omp_get_thread_num()], f_vector);
                //std::cout << omp_get_thread_num() << "\n";
            }
            //double test2 = omp_get_wtime() - test1;
            //std::cout << "loop time: " << test2 << "\n";
        }
    }
}

/*-----------computing the f vector components related to this Neumann BC---------*/
void computeNeumannBC(const int numNodes, const int i, const std::vector<std::vector<std::size_t>> & elementTags,
                        const double tx, const double ty, std::map<int, int> & nodeIndexMap,
                        std::vector<double> & f_vector,
                        const std::vector<double> & weights, const std::vector<double> & basisFunctions, 
                        const std::vector<double> & determinants, const std::vector<std::vector<std::size_t>> & nodeTags,
                        std::vector<double> & full_f_vector)
{
    int first_node_index;
    int first_node;
    int my_row;

    for (std::size_t el = 0; el < elementTags[i].size(); el++){
        for (int j = 0; j < 2*numNodes; ++j){
            f_vector[j] = 0.;
        }
    // looping over the Gauss points (formula of slide 17)
        for (std::size_t j = 0; j < weights.size(); ++j){
            // looping over the number of shape functions
            for (int l = 0; l < numNodes; ++l){ 
                f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*j]*tx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
                f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*j]*ty*determinants[j+weights.size()*el]*weights[j]; // le weight non plus !
            }
        }
    /*------------- ASSEMBLY PROCESS ------------*/
        for (int j = 0; j < 2*numNodes; ++j){
            first_node_index = (int) j/2;
            first_node = nodeTags[i][numNodes*el+first_node_index];
            my_row = nodeIndexMap[first_node] + j % 2;
            full_f_vector[my_row] += f_vector[j];
        }
    }    
}

/*-----------focus on neumann BC's (surface traction) into full_f_vector------------*/
void includeNeumannBC(const std::vector<std::string> & keys, const std::string & integration_rule,
                std::map<int, int> & nodeIndexMap,
                std::map<std::string, std::pair<int, int>> & groups,
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
    std::vector<double> jacobians, determinants, coords;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;

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
            std::string groupname = ""; // will remain empty only if tx == ty == 0
            if(words[2].compare("tx") == 0){
                tx = value[0];
                std::string new_key = "Boundary Conditions/" + words[1] + "/ty";
                gmsh::onelab::getNumber(new_key, value);
                ty = value[0];
                if(tx || ty) // there is no need the consider Neumann BC if it is homogeneous
                    groupname = words[1];
            }
            //here, if we are dealing with neumann BC's, we should have retrieved a non-empty group_name
            if(groupname.compare("") != 0){
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
                                                                integration_rule, //maybe use different integration method on boundary than in domain
                                                                localCoords, weights);

                        // get determinants and Jacobians (only determinants will be useful)
                        gmsh::model::mesh::getJacobians(elementTypes[i], 
                                                        localCoords, jacobians, 
                                                        determinants, coords,tags[k]); // ai dû rajouter mon tags[k] pr focus on this particular entity

                        gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                            "Lagrange", numComponents, 
                                                            basisFunctions, 
                                                            numOrientations);

                        /*-----------computing the f vector components related to this Neumann BC---------*/
                        //looping over the elements
                        std::vector<double> f_vector(2*numNodes);
                        computeNeumannBC(numNodes, i, elementTags, tx, ty, nodeIndexMap,
                                           f_vector, weights, basisFunctions, 
                                           determinants, nodeTags, full_f_vector);
                    }
                }
            }
        }
    }
}

// looping over the boundary conditions imposed in the .geo file, the purpose of the whole loop is
// to fill the two vectors defined just above
void fillDOFindicesDirBC(const std::vector<std::string> & keys, 
                             std::map<int, int> & nodeIndexMap,
                             std::map<std::string, std::pair<int, int>> & groups, 
                             std::vector<int> & new_DOFindices,
                             std::vector<double> & dir_BC)
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
}

void computeReactionForces(std::map<std::string, std::pair<int, int>> & groups,
                     std::map<int, int> & nodeIndexMap, std::vector<double> & final_nodal_forces,
                     const std::vector<std::string> & keys)
{
    double Rx;
    double Ry;

    int dim;
    int tag;
    std::vector<int> tags;

    std::vector<std::size_t> nodeTags;
    std::vector<double> coord; //useless but mandatory
    std::vector<double> parametricCoord; //useless but mandatory

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
                            int my_new_col = nodeIndexMap[nodeTags[j]]; // le noeud numéro nodeIndexMap[nodeTags[j]] correspond à la ligne/colonne dans K, f
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

void plotUxUy(const std::vector<std::size_t> & FEM_nodeTags, const std::vector<double> & full_d_vector, const int nbViews, std::vector<std::string> & names)
{
    int viewtag1 = gmsh::view::add("ux");
    int viewtag2 = gmsh::view::add("uy");
    std::vector<std::vector<double>> data1(FEM_nodeTags.size());
    std::vector<std::vector<double>> data2(FEM_nodeTags.size());
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        data1[i].resize(1);
        data1[i][0] = full_d_vector[2*i]; //ux
        data2[i].resize(1);
        data2[i][0] = full_d_vector[2*i+1]; //uy
    }

    gmsh::view::addModelData(viewtag1, 0, names[0], "NodeData",
                            FEM_nodeTags, data1); //independent of time here
    gmsh::view::addModelData(viewtag2, 0, names[0], "NodeData",
                            FEM_nodeTags, data2); //independent of time here

    // save the results to disk
    //gmsh::view::write(viewtag1, "ux.msh");
    //gmsh::view::write(viewtag2, "uy.msh");

    // set display options & view results
    gmsh::option::setNumber("View[" + std::to_string(nbViews) + "].IntervalsType", 3);
    gmsh::option::setNumber("View[" + std::to_string(nbViews) + "].AdaptVisualizationGrid", 1);
    gmsh::option::setNumber("View[" + std::to_string(nbViews) + "].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[" + std::to_string(nbViews) + "].TargetError", -0.0001);
    gmsh::option::setNumber("View[" + std::to_string(nbViews +1) + "].IntervalsType", 3);
    gmsh::option::setNumber("View[" + std::to_string(nbViews +1) + "].AdaptVisualizationGrid", 1);
    gmsh::option::setNumber("View[" + std::to_string(nbViews +1) + "].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[" + std::to_string(nbViews +1) + "].TargetError", -0.0001);
}

void computeStressStrainElNodal(const std::vector<int> & elementTypes, const std::vector<std::vector<std::size_t>> &elementTags,
                           const std::vector<std::vector<std::size_t>> &elnodeTags,
                           std::vector<std::size_t> & elems2D, const Eigen::Matrix<double, 3, 3> &H_matrix, 
                           const double nu, const double E, std::map<int, elementData> & FEM_kinematicsMap,
                           std::vector<std::vector<double>> & data3, std::vector<std::vector<double>> & data4,
                           std::vector<std::vector<double>> & data5, std::vector<std::vector<double>> & data6,
                           std::vector<std::vector<double>> & data7, std::vector<std::vector<double>> & data8,
                           std::vector<std::vector<double>> & data9, std::vector<std::vector<double>> & data10,
                           int & k)
{
    // The strains and stresses are computed in local axes
    // They are then converted back into global axes, using a simple rotation matrix

    std::string name;
    int dim, order, numNodes, numPrimaryNodes;
    std::vector<double> paramCoord2D;
    std::vector<double> jacobians, determinants, coords;
    int numComponents, numOrientations;
    std::vector<double> basisFunctionsGradient;

    Eigen::Vector3d epsilon; // will contain local values of the strains (Voigt notation)
    Eigen::Vector3d sigma; // will contain local values of the stresses (Voigt notation)

    Eigen::Matrix<double, 2, 2> jacob2D;
    Eigen::Matrix<double, 2, 2> jacobinv;
    Eigen::Matrix<double, 2, 2> jacobinvtrans;

    elementData* element;
    double theta;
    double dN_dx_local;
    double dN_dy_local;

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
        gmsh::model::mesh::getJacobians(elementTypes[ityp], 
                                        paramCoord, jacobians, 
                                        determinants, coords);

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
            elems2D[k] = etags[i];

            theta = element->theta;
            rotation_tensor(0,0) = cos(theta);
            rotation_tensor(0,1) = sin(theta);
            rotation_tensor(1,0) = -sin(theta);
            rotation_tensor(1,1) = cos(theta);

            data3[k].resize(numNodes);
            data4[k].resize(numNodes);
            data5[k].resize(numNodes);
            data6[k].resize(numNodes);
            data7[k].resize(numNodes);
            data8[k].resize(numNodes);
            data9[k].resize(numNodes);
            data10[k].resize(numNodes);

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

                computeBmatrix(j, numNodes, basisFunctionsGradient, determinants, 
                                 jacobinvtrans, H_matrix, B_matrix);
                
                // WARNING, the B_matrix is expressed in terms of global axes -> convert back to local axes
                for(int k = 0; k < numNodes; k++)
                {
                    dN_dx_local = cos(theta)*B_matrix(0,2*k) + sin(theta)*B_matrix(2,2*k);
                    dN_dy_local = -sin(theta)*B_matrix(0,2*k) + cos(theta)*B_matrix(2,2*k);
                    B_matrix(0,2*k) = dN_dx_local;
                    B_matrix(2,2*k) = dN_dy_local;
                    B_matrix(1,2*k+1) = dN_dy_local;
                    B_matrix(2,2*k+1) = dN_dx_local;
                }
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

                //converting back to global axes
                strain_tensor = rotation_tensor.transpose()*strain_tensor*rotation_tensor;
                stress_tensor = rotation_tensor.transpose()*stress_tensor*rotation_tensor;

                data3[k][j] = strain_tensor(0,0);
                data4[k][j] = strain_tensor(1,1);
                data5[k][j] = strain_tensor(1,0);
                data6[k][j] = -nu/E*(stress_tensor(0,0)+stress_tensor(1,1));// formula slide 11/20 elasticity
                data7[k][j] = stress_tensor(0,0);
                data8[k][j] = stress_tensor(1,1);
                data9[k][j] = stress_tensor(1,0);
                data10[k][j] = sqrt(stress_tensor(0,0)*stress_tensor(0,0) + stress_tensor(1,1)*stress_tensor(1,1) - stress_tensor(0,0)*stress_tensor(1,1) + 3*stress_tensor(1,0)*stress_tensor(1,0));
            }
            k++;
        }
    }
}

void applyElecPressure(const std::string & integration_rule, std::map<int,std::pair<double,double>> & nodeCoordsMap, std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups, std::map<int, double> &electrostaticPressure, std::vector<double> & full_f_vector)
{
    std::string groupname = "BEM_FEM_boundary";
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;

    std::string elementName;
    int elDim, order, numNodes, numPrimaryNodes;
    std::vector<double> localNodeCoord;
    std::vector<double> localCoords, weights;
    std::vector<double> jacobians, determinants, coords;
    int numComponents, numOrientations;
    std::vector<double> basisFunctions;

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
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(std::size_t j=0; j< elementTypes.size(); ++j)
        {

            // getting element properties
            gmsh::model::mesh::getElementProperties(elementTypes[j],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                    
            // get Gauss points coordinates and weights
            gmsh::model::mesh::getIntegrationPoints(elementTypes[j], 
                                                    integration_rule, //maybe use different integration method on boundary than in domain
                                                    localCoords, weights);

            // get determinants and Jacobians (only determinants will be useful)
            gmsh::model::mesh::getJacobians(elementTypes[j], 
                                            localCoords, jacobians, 
                                            determinants, coords,tags[i]); // ai dû rajouter mon tags[k] pr focus on this particular entity

            gmsh::model::mesh::getBasisFunctions(elementTypes[j], localCoords, 
                                                "Lagrange", numComponents, 
                                                basisFunctions, 
                                                numOrientations);

            /*-----------computing the f vector components related to this Neumann BC---------*/
            std::vector<double> f_vector(2*numNodes);

            for(std::size_t k=0; k < elementTags[j].size(); ++k)
            {
                first_nodeTag = elnodeTags[j][numNodes*k];
                second_nodeTag = elnodeTags[j][numNodes*k+1];
                nx = -(nodeCoordsMap[second_nodeTag].second - nodeCoordsMap[first_nodeTag].second);
                ny = (nodeCoordsMap[second_nodeTag].first - nodeCoordsMap[first_nodeTag].first);
                norm = sqrt(nx*nx+ny*ny);
                nx = nx / norm; //normalisation
                ny = ny / norm;
                p = electrostaticPressure[elementTags[j][k]];
                tx = p*nx;
                ty = p*ny;

                for (int n = 0; n < 2*numNodes; ++n){
                    f_vector[n] = 0.;
                }
                // looping over the Gauss points (formula of slide 17)
                for (std::size_t n = 0; n < weights.size(); ++n){
                    // looping over the number of shape functions
                    for (int l = 0; l < numNodes; ++l){ 
                        f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*n]*tx*determinants[n+weights.size()*k]*weights[n]; //pas oublier le dét !
                        f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*n]*ty*determinants[n+weights.size()*k]*weights[n]; // le weight non plus !
                    }
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                for (int n = 0; n < 2*numNodes; ++n){
                    first_node_index = (int) n/2;
                    first_node = elnodeTags[j][numNodes*k+first_node_index];
                    my_row = nodeIndexMap[first_node] + n % 2;
                    full_f_vector[my_row] += f_vector[n];
                }
            }
        }
    }
}

int initKinematics(std::vector<std::size_t> & FEM_elementTags, std::map<int, elementData> & FEM_kinematicsMap, std::map<int,std::pair<double,double>> & FEM_nodeCoordsMap, int & dim, std::vector<int> & tags)
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

void updateKinematics(std::map<int, elementData> & FEM_kinematicsMap, std::vector<double> & full_d_vector, std::map<int, int> & nodeIndexMap, std::vector<std::size_t> & FEM_elementTags)
{
    elementData* element;
    size_t numNodes;
    size_t tmp_nodeTag;
    int node_index;
    double numerator;
    double denominator;
    double angle1;
    double angle2;
    double norm1;
    double norm2;
    for(std::size_t i=0; i< FEM_elementTags.size(); ++i)
    {
        element = & FEM_kinematicsMap[FEM_elementTags[i]];
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
        denominator =0;
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
        // filling ui and vi while finding the angle minimizing equation (4) in the paper (what i call the norm)
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

void updateKinMatrices(elementData* element, bool computeFl)
{
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

    //update B matrix
    element->B = element->P*element->E.transpose();

    //update fl vector, only if Kl has been initialized
    if(computeFl)
    {
        element->fl = element->Kl*element->pl;
    }
}

void Assembly_global_f(std::vector<double> & global_f_vector, std::map<int, elementData> & FEM_kinematicsMap, std::map<int, int> &  nodeIndexMap, std::vector<int> & new_DOFindices, std::vector<std::size_t> & FEM_elementTags)
{
    elementData* element;
    Eigen::Matrix<double, Eigen::Dynamic, 1> fg;
    int numNodes;
    int index1;
    int index2;
    int nodeTag;
    for(size_t i = 0; i < FEM_elementTags.size(); i++)
    {
        element = & FEM_kinematicsMap[FEM_elementTags[i]];
        numNodes = (int) element->nodes.size();
        fg = element->B.transpose()*element->fl;
        for(int j = 0; j < numNodes; j++)
        {
            nodeTag = element->nodes[j];
            index1 = new_DOFindices[nodeIndexMap[nodeTag]]; // current index related to x DOF of the node (= -1 if dirichlet BC)
            index2 = new_DOFindices[nodeIndexMap[nodeTag] + 1]; // current index related to y DOF of the node (= -1 if dirichlet BC)
            if(index1 >= 0) // need to assemble in global f vector
            {
                global_f_vector[index1] += fg(2*j);
            }
            if(index2 >= 0)
            {
                global_f_vector[index2] += fg(2*j + 1);
            }
        }
    }
}

double computeResidualNorm(std::vector<double> & global_f_vector, std::vector<double> & f_app_vector, Eigen::Matrix<double, Eigen::Dynamic, 1> & residual)
{
    // norm computed as a geometrical norm (one could also use an other norm)
    double res = 0;
    for(size_t i = 0; i < global_f_vector.size(); i++)
    {
        residual(i,0) = global_f_vector[i] - f_app_vector[i];
        res += residual(i,0)*residual(i,0);
    }
    res = sqrt(res);
    return res;
}

void Assembly_global_Kg(Eigen::SparseMatrix<double> & global_tangent_matrix, std::map<int, elementData> & FEM_kinematicsMap, std::map<int, int> &  nodeIndexMap, std::vector<int> & new_DOFindices, std::vector<std::size_t> & FEM_elementTags)
{
    std::list<Eigen::Triplet <double>> Kg_tripletList;

    elementData* element;
    Eigen::Matrix<double, Eigen::Dynamic, 1> fg;
    int numNodes;
    Eigen::Matrix<double, Eigen::Dynamic, 1> tmp_F;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Kh;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Kg;

    int first_node_index;
    int first_node;
    int my_row;
    int second_node_index;
    int second_node;
    int my_col;
    for(size_t i = 0; i < FEM_elementTags.size(); i++)
    {
        element = & FEM_kinematicsMap[FEM_elementTags[i]];
        numNodes = (int) element->nodes.size();
        tmp_F = element->P.transpose()*element->fl;
        /*if(FEM_elementTags[i] == 16)
            std::cout << "this is the first tmp F matrix: " << tmp_F << "\n";*/

        Eigen::Matrix<double, 1, Eigen::Dynamic> F(1,2*numNodes);
        for(int j = 0; j < numNodes; j++)
        {
            F(0,2*j) = tmp_F(2*j+1, 0);
            F(0,2*j+1) = -tmp_F(2*j, 0);
        }
        /*if(FEM_elementTags[i] == 16)
            std::cout << "this is the first F matrix: " << F << "\n";*/
        
        Kh = element->E*(-F.transpose()*element->G - element->G.transpose()*F*element->P)*element->E.transpose();

        /*if(FEM_elementTags[i] == 16)
            std::cout << "this is the first Kh matrix: " << Kh << "\n";*/

        Kg = element->B.transpose()*element->Kl*element->B + Kh;
        /*if(FEM_elementTags[i] == 16)
            std::cout << "this is the first Kg matrix: " << Kg << "\n";*/

        // assemble Kg into tripletlist
        for (int j = 0; j < 2*numNodes; ++j)
        {
            first_node_index = (int) j/2;
            first_node = element->nodes[first_node_index];
            my_row = new_DOFindices[nodeIndexMap[first_node] + j % 2]; // node list begins at 1, '-1' correction to begin at 0
            if(my_row >= 0) // else, the DOF corresponds to a Dirichlet BC
            {
                for(int l = 0; l <= j; ++l)
                {   // to exploit symmetry
                    second_node_index = (int) l/2;
                    second_node = element->nodes[second_node_index];
                    my_col = new_DOFindices[nodeIndexMap[second_node] + l % 2];

                    if(my_col >= 0) // else, the DOF corresponds to a Dirichlet BC
                    {
                        Kg_tripletList.push_back(Eigen::Triplet<double>(my_row,my_col,Kg(j,l)));
                        if(j != l)
                        {
                            Kg_tripletList.push_back(Eigen::Triplet<double>(my_col,my_row,Kg(j,l)));
                        }
                    }
                }
            }
        }
    }
    global_tangent_matrix.setFromTriplets(Kg_tripletList.begin(), Kg_tripletList.end());
}
void updateNodalDispVector(std::vector<double> & full_d_vector, std::vector<double> & global_p_vector, std::vector<int> & new_DOFindices)
{
    for(size_t i = 0; i < full_d_vector.size(); i++)
    {
        if(new_DOFindices[i] >= 0)
        {
            full_d_vector[i] = global_p_vector[new_DOFindices[i]];
        }
    }
}

void retrieveFinalNodalForces(std::vector<double> & final_nodal_forces, std::map<int, elementData> & FEM_kinematicsMap, std::map<int, int> & nodeIndexMap, std::vector<std::size_t> & FEM_elementTags)
{
    elementData* element;
    Eigen::Matrix<double, Eigen::Dynamic, 1> fg;
    int numNodes;
    int nodeTag;
    for(size_t i = 0; i < FEM_elementTags.size(); i++)
    {
        element = & FEM_kinematicsMap[FEM_elementTags[i]];
        numNodes = (int) element->nodes.size();
        fg = element->B.transpose()*element->fl;
        for(int j = 0; j < numNodes; j++)
        {
            nodeTag = element->nodes[j];
            final_nodal_forces[nodeIndexMap[nodeTag]] += fg(2*j);
            final_nodal_forces[nodeIndexMap[nodeTag] + 1] += fg(2*j + 1);
        }
    }
}

void retrieveBoundaryDisplacements(std::map<int,std::pair<double,double>> & boundaryDisplacementMap, const std::vector<double> & full_d_vector, std::map<int, int> & nodeIndexMap, std::map<std::string, std::pair<int, int>> & groups)
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
                    //std::cout << "nodeTag : " << nodeTag << ": u = " << u << ", v = " << v << "\n";
                }
            }
        }
    }
}
