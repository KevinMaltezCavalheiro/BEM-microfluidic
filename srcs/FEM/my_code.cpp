#include <gmsh.h>       // gmsh main header
#include <iostream>     // for std::cout
#include <map>
#include <sstream>  // for std::stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices
#include <Eigen/SparseCholesky> // solving sparse linear systems
#include <cmath>

double compute_f_value(const int, const Eigen::SparseMatrix<double> &, const std::vector<double> &);

//Compute the matrix H
Eigen::Matrix<double, 3, 3> Compute_H_matrix(const double E, const double nu){
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

//Compute the full matrix "K" and "f".
void Assembly_K_triplet_F(const int i, const int el, 
                          const std::vector<std::vector<std::size_t>> nodeTags, 
                          const int numNodes, std::vector<Eigen::Triplet <double>> & k_tripletList,
                          std::vector<double> & full_f_vector, 
                          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & K_matrix,
                          std::vector<double> & f_vector){

    for (int j = 0; j < 2*numNodes; ++j){
        int first_node_index = (int) j/2;
        int first_node = nodeTags[i][numNodes*el+first_node_index];
        int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
        for(int l = 0; l <= j; ++l){ // to exploit symmetry
            int second_node_index = (int) l/2;
            int second_node = nodeTags[i][numNodes*el+second_node_index];
            int my_col = 2*(second_node-1) + l % 2;

            //full_K_matrix.coeffRef(my_row,my_col) += K_matrix(j,l); //utiliser triplets puis reessayer avec noeuds ordre 2
            k_tripletList.push_back(Eigen::Triplet<double>(my_row,my_col,K_matrix(j,l)));
            if(j != l){
                //full_K_matrix.coeffRef(my_col,my_row) += K_matrix(l,j); // use symmetry
                k_tripletList.push_back(Eigen::Triplet<double>(my_col,my_row,K_matrix(j,l)));
            }
        }
        full_f_vector[my_row] += f_vector[j];
    }
}


void Compute_B_matrix(const int j, const int numNodes, const std::vector<double> basisFunctionsGradient,
                      const std::vector<double> determinants, 
                      const Eigen::Matrix<double, 2, 2> jacobinvtrans,
                      const Eigen::Matrix<double, 3, 3> H_matrix,
                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & B_matrix){

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
void Get_physical_Prop_Domain(const std::vector<std::string> keys, double & E, double & nu, double & rho,
                              double & bx, double & by){
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
/*------------computing K matrix and f vector------------------------*/
// we first consider the K matrix and the volumic part of the f vector

void Loop_K_F(const Eigen::Matrix<double, 3, 3> H_matrix, const double rho, const double bx, const double by,
              const std::string integration_rule, 
              int & dim, int & tag, std::vector<int> & tags,
              Eigen::SparseMatrix<double> & full_K_matrix, 
              std::vector<double> & full_f_vector,
              std::vector<Eigen::Triplet <double>> & k_tripletList){

    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
    for (std::size_t k = 0; k < tags.size(); ++k)
    {
        //std::cout << "Entity (" << dim << "," << tags[k] << "):\n";
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types
        for (std::size_t i = 0; i < elementTypes.size(); ++i)
        {

            // getting element properties
            std::string elementName;
            int dim, order, numNodes, numPrimaryNodes;
            std::vector<double> localNodeCoord;
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, dim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);
            
            // get Gauss points coordinates and weights
            std::vector<double> localCoords, weights;
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    integration_rule, // "Gauss2"
                                                    localCoords, weights);

            // get determinants and Jacobians
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords);

            // values and derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            int numComponents, numOrientations;
            std::vector<double> basisFunctions;
            std::vector<double> basisFunctionsGradient;
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
            for (std::size_t el = 0; el < elementTags[i].size(); el++){
                // elemental K matrix and f vector initialization
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K_matrix(2*numNodes,2*numNodes);
                std::vector<double> f_vector(2*numNodes);
                for (int j = 0; j < 2*numNodes; ++j){
                    for (int l = 0; l < 2*numNodes; ++l){
                        K_matrix(j,l)=0.;
                    }
                    f_vector[j] = 0.;
                }
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


                    Compute_B_matrix(j, numNodes, basisFunctionsGradient, determinants, 
                                     jacobinvtrans, H_matrix, B_matrix);
                    
                    for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
                        f_vector[2*l] = f_vector[2*l] + rho*basisFunctions[l+numNodes*j]*bx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
                        f_vector[2*l+1] = f_vector[2*l+1] + rho*basisFunctions[l+numNodes*j]*by*determinants[j+weights.size()*el]*weights[j]; // pas oublier le weights[j] non plus !
                    }

                    K_matrix = K_matrix + B_matrix.transpose()*H_matrix*B_matrix*determinants[j+weights.size()*el]*weights[j];
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                // USE LIST OF TRIPLETS TO PUSH ONLY ONCE -> more rapid and easier to parallelize
                Assembly_K_triplet_F(i, el, nodeTags, numNodes, k_tripletList, full_f_vector, K_matrix, f_vector);
            }
        }
    }
    full_K_matrix.setFromTriplets(k_tripletList.begin(), k_tripletList.end());
}

/*-----------computing the f vector components related to this Neumann BC---------*/
void Compute_Neumann_BC(const int numNodes, const int i, const std::vector<std::vector<std::size_t>> elementTags,
                        const double tx, const double ty, std::map<int, int> & node_index,
                        std::vector<double> & f_vector,
                        const std::vector<double> weights, const std::vector<double> & basisFunctions, 
                        const std::vector<double> determinants, const std::vector<std::vector<std::size_t>> nodeTags,
                        std::vector<double> & full_f_vector){

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
            int first_node_index = (int) j/2;
            int first_node = nodeTags[i][numNodes*el+first_node_index];
            int my_row = node_index[first_node] + j % 2;
            full_f_vector[my_row] += f_vector[j];
        }
    }    
}

/*-----------focus on neumann BC's (surface traction) into full_f_vector------------*/
void Neumann_BC(const std::vector<std::string> keys, const std::string integration_rule,
                int & dim, int & tag, std::vector<int> & tags, std::map<int, int> & node_index,
                std::map<std::string, std::pair<int, int>> & groups,
                std::string & groupname,
                std::vector<double> & full_f_vector){
    double tx = 0; // surface traction along x axis
    double ty = 0; // surface traction along y axis
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
            groupname = ""; // will remain empty only if tx == ty == 0
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
                    std::vector<int> elementTypes;
                    std::vector<std::vector<std::size_t>> elementTags;
                    std::vector<std::vector<std::size_t>> nodeTags;
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (std::size_t i = 0; i < elementTypes.size(); ++i)
                    {
                        // getting element properties
                        std::string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        std::vector<double> localNodeCoord;
                        gmsh::model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, dim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                    
                        // get Gauss points coordinates and weights
                        std::vector<double> localCoords, weights;
                        gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                                integration_rule, //maybe use different integration method on boundary than in domain
                                                                localCoords, weights);

                        // get determinants and Jacobians (only determinants will be useful)
                        std::vector<double> jacobians, determinants, coords;
                        gmsh::model::mesh::getJacobians(elementTypes[i], 
                                                        localCoords, jacobians, 
                                                        determinants, coords,tags[k]); // ai dû rajouter mon tags[k] pr focus on this particular entity

                        int numComponents, numOrientations;
                        std::vector<double> basisFunctions;
                        gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                            "Lagrange", numComponents, 
                                                            basisFunctions, 
                                                            numOrientations);

                        /*-----------computing the f vector components related to this Neumann BC---------*/
                        //looping over the elements
                        std::vector<double> f_vector(2*numNodes);
                        Compute_Neumann_BC(numNodes, i, elementTags, tx, ty, node_index,
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
void Fill_New_indices_Dir_BC(const std::vector<std::string> keys, 
                             int & dim, int & tag, std::vector<int> & tags, std::map<int, int> & node_index,
                             std::map<std::string, std::pair<int, int>> & groups, 
                             std::string & groupname, std::vector<int> & new_indices,
                             std::vector<double> & dir_BC){

    for (std::size_t j = 0; j < new_indices.size(); ++j){
        new_indices[j] = (int) j;
    } 
    for (std::size_t j = 0; j < dir_BC.size(); ++j){
        dir_BC[j] = 0;
    }

    double ux = 0; // imposed displacement along x axis
    double uy = 0; // imposed displacement along y axis
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
                if(words[2].compare("ux") == 0)
                    ux = value[0];
                if(words[2].compare("uy") == 0)
                    uy = value[0];
                groupname = words[1];
                
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
                    std::vector<int> elementTypes;
                    std::vector<std::vector<std::size_t>> elementTags;
                    std::vector<std::vector<std::size_t>> nodeTags;
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (std::size_t i = 0; i < elementTypes.size(); ++i)
                    {

                        // getting element properties
                        std::string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        std::vector<double> localNodeCoord;
                        gmsh::model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, dim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                        //looping over the elements
                        for (std::size_t el = 0; el < elementTags[i].size(); el++){
                            for (int j = 0; j < 2*numNodes; ++j){
                                int first_node_index = (int) j/2;
                                int first_node = nodeTags[i][numNodes*el+first_node_index];
                                int my_row = node_index[first_node] + j % 2; // j % 2 == 1 corresponds to y nodal displacement
                                if(words[2].compare("uy") == 0 && j % 2){ // we are dealing with a uy type of node
                                    dir_BC[my_row] = uy;
                                    new_indices[my_row] = -1;
                                }
                                if(words[2].compare("ux") == 0 && j % 2 == 0){ // ux type of node
                                    dir_BC[my_row] = ux;
                                    new_indices[my_row] = -1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Fill_K_F_Reduced(const std::vector<double> & full_f_vector, const Eigen::SparseMatrix<double> full_K_matrix,
                      const std::vector<int> new_indices, const std::vector<double> dir_BC,
                      Eigen::SparseMatrix<double> & reduced_K_matrix, 
                      std::vector<double> & reduced_f_vector,
                      std::vector<Eigen::Triplet <double>> & reduced_k_tripletList){
    for (std::size_t j = 0; j < reduced_f_vector.size(); ++j){
        reduced_f_vector[j] = 0.;
    }

    // looping over the non-zero elements of full_K_matrix
    int row;
    int new_row;
    int new_col;
    double k_value;
    for (int col=0; col<full_K_matrix.outerSize(); ++col){
        new_col = new_indices[col];
        for (Eigen::SparseMatrix<double>::InnerIterator it(full_K_matrix,col); it; ++it){
            k_value = it.value();
            row = it.row();   // row index
            new_row = new_indices[row];
            if(new_row >= 0 && new_col >= 0){ // if the given element must not be removed
                if(new_row < new_col){ // taking advantage of the symmetry
                    //reduced_K_matrix.coeffRef(new_row,new_col) = k_value;
                    //reduced_K_matrix.coeffRef(new_col,new_row) = k_value;
                    reduced_k_tripletList.push_back(Eigen::Triplet<double>(new_row,new_col,k_value));
                    reduced_k_tripletList.push_back(Eigen::Triplet<double>(new_col,new_row,k_value));
                }
                if(new_row == new_col){ // diagonal element
                    //reduced_K_matrix.coeffRef(new_row,new_col) = k_value;
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

void Reaction_Forces(double & Rx, double & Ry, int & dim, int & tag, 
                     std::vector<int> & tags, std::string & groupname,
                     std::map<std::string, std::pair<int, int>> & groups,
                     std::map<int, int> & node_index, Eigen::VectorXd nodal_forces,
                     const std::vector<std::string> keys){

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
                groupname = words[1];
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
                    std::vector<std::size_t> nodeTags;
                    std::vector<double> coord; //useless
                    std::vector<double> parametricCoord; //useless

                    // retrieving all the nodes in this particular entity
                    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tags[k], true);
                    
                    
                    for (std::size_t j = 0; j < nodeTags.size(); ++j){
                        if(words[2].compare("ux") == 0){    // computing Rx
                            int my_new_col = node_index[nodeTags[j]]; // le noeud numéro node_index[nodeTags[j]] correspond à la ligne/colonne dans K, f
                            Rx += nodal_forces(my_new_col);
                        }
                        if(words[2].compare("uy") == 0){    // computing Ry
                            int my_new_col = node_index[nodeTags[j]] + 1;
                            Ry += nodal_forces(my_new_col);
                        }
                    }
                }
                std::cout << "Rx: " << Rx << " [N/m], Ry: " << Ry << " [N/m] \n";
            }
        }
    }
}

std::vector<std::string> Plot_Save_ux_uy(const std::vector<std::size_t> all_nodeTags, const std::vector<double> full_d_vector){
    int viewtag1 = gmsh::view::add("ux");
    int viewtag2 = gmsh::view::add("uy");
    std::vector<std::vector<double>> data1(all_nodeTags.size());
    std::vector<std::vector<double>> data2(all_nodeTags.size());
    for (std::size_t i = 0; i < all_nodeTags.size(); ++i)
    {
        data1[i].resize(1);
        data1[i][0] = full_d_vector[2*i]; //ux
        data2[i].resize(1);
        data2[i][0] = full_d_vector[2*i+1]; //uy
    }

    std::vector<std::string> names;
    gmsh::model::list(names);
    gmsh::view::addModelData(viewtag1, 0, names[0], "NodeData",
                            all_nodeTags, data1); //independent of time here
    gmsh::view::addModelData(viewtag2, 0, names[0], "NodeData",
                            all_nodeTags, data2); //independent of time here

    // save the results to disk
    gmsh::view::write(viewtag1, "ux.msh");
    gmsh::view::write(viewtag2, "uy.msh");

    // set display options & view results
    gmsh::option::setNumber("View[0].IntervalsType", 3);
    gmsh::option::setNumber("View[0].AdaptVisualizationGrid", 1);
    gmsh::option::setNumber("View[0].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[0].TargetError", -0.0001);
    gmsh::option::setNumber("View[1].IntervalsType", 3);
    gmsh::option::setNumber("View[1].AdaptVisualizationGrid", 1);
    gmsh::option::setNumber("View[1].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[1].TargetError", -0.0001);
    return names;
}

void Compute_Stress_Strain_Elemental_Nodal(const std::vector<int> elementTypes, const std::vector<std::vector<std::size_t>> elementTags,
                           const std::vector<std::vector<std::size_t>> elnodeTags,
                           const std::vector<double> full_d_vector, std::vector<std::size_t> & elems2D,
                           const Eigen::Matrix<double, 3, 3> H_matrix, const double nu, const double E,
                           std::vector<int> & nb_adjacent_elements, std::vector<double> & nodal_sigma_x, 
                           std::vector<double> & nodal_sigma_y, std::vector<double> & nodal_sigma_xy,
                           std::map<int, int> & node_index,
                           std::vector<std::vector<double>> & data3, std::vector<std::vector<double>> & data4,
                           std::vector<std::vector<double>> & data5, std::vector<std::vector<double>> & data6,
                           std::vector<std::vector<double>> & data7, std::vector<std::vector<double>> & data8,
                           std::vector<std::vector<double>> & data9, std::vector<std::vector<double>> & data10){

    int k = 0;

    for (std::size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
    {
        // get info about this type of element
        std::string name;
        int dim, order, numNodes, numPrimaryNodes;
        std::vector<double> paramCoord2D;
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
                                        determinants, coords);

        // derivatives of the shape functions at the nodes of the reference element
        // we are here working with derivatives in the reference space, same for all elements
        int numComponents, numOrientations;
        std::vector<double> basisFunctionsGradient;

        gmsh::model::mesh::getBasisFunctions(elementTypes[ityp], paramCoord, 
                                                "GradLagrange", numComponents, 
                                                basisFunctionsGradient, 
                                                numOrientations);

        Eigen::Vector3d epsilon; // will contain local values of the strains
        Eigen::Vector3d sigma; // will contain local values of the stresses

        // loop over elements of type "ityp"
        for (std::size_t i = 0; i < etags.size(); ++i)
        {
            elems2D[k] = etags[i];

            data3[k].resize(numNodes);
            data4[k].resize(numNodes);
            data5[k].resize(numNodes);
            data6[k].resize(numNodes);
            data7[k].resize(numNodes);
            data8[k].resize(numNodes);
            data9[k].resize(numNodes);
            data10[k].resize(numNodes);


            // loop over nodes of the element, filling the local nodal displacement vector "d"
            for (int j = 0; j < numNodes; ++j)
            {
                int nodeTag = enods[numNodes * i + j];
                nodal_displacement(2*j) = full_d_vector[node_index[nodeTag]];
                nodal_displacement(2*j+1) = full_d_vector[node_index[nodeTag]+1];

                nb_adjacent_elements[nodeTag-1] +=1; // -1 is because nodeTags begin at 1
            }
            // looping over the nodes
            for (int j = 0; j < numNodes; ++j){
                int nodeTag = enods[numNodes * i + j];
                // converting jacobian to Eigen format
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                    jacob(&jacobians[9*j+9*numNodes*i], 3, 3);
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

                // computing the B matrix (formula slide 19) 
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);

                Compute_B_matrix(j, numNodes, basisFunctionsGradient, determinants, 
                                 jacobinvtrans, H_matrix, B_matrix);

                epsilon = B_matrix*nodal_displacement;
                sigma = H_matrix*epsilon;

                data3[k][j] = epsilon(0);
                data4[k][j] = epsilon(1);
                data5[k][j] = epsilon(2);
                data6[k][j] = -nu/E*(sigma(0)+sigma(1));// formula slide 11/20 elasticity
                data7[k][j] = sigma(0);
                data8[k][j] = sigma(1);
                data9[k][j] = sigma(2);
                data10[k][j] = sqrt(sigma(0)*sigma(0) + sigma(1)*sigma(1) - sigma(0)*sigma(1) + 3*sigma(2)*sigma(2));


                nodal_sigma_x[nodeTag - 1] += sigma(0);
                nodal_sigma_y[nodeTag - 1] += sigma(1);
                nodal_sigma_xy[nodeTag - 1] += sigma(2);
            }
            k++;
        }
    }
    // computing mean nodal values (not sure if needed)
    for (std::size_t i = 0; i < nodal_sigma_x.size(); ++i){
        nodal_sigma_x[i] = nodal_sigma_x[i] / nb_adjacent_elements[i];
        nodal_sigma_y[i] = nodal_sigma_y[i] / nb_adjacent_elements[i];
        nodal_sigma_xy[i] = nodal_sigma_xy[i] / nb_adjacent_elements[i];
    }
}

double compute_f_value(const int col, const Eigen::SparseMatrix<double> & full_K_matrix, const std::vector<double> & full_d_vector){
    int row;
    double k_value;
    double f_value = 0.;
    for (Eigen::SparseMatrix<double>::InnerIterator it(full_K_matrix,col); it; ++it){
        k_value = it.value();
        row = it.row();   // row index
        f_value += k_value*full_d_vector[row];
    }
    return f_value;
}

int main(int argc, char **argv)
{
    /*--------------INIT----------------*/
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    gmsh::initialize();

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);



    /*--------Integration rule used in the code----------*/
    // could be useful to pass it as an argument in the .geo file
    std::string integration_rule = "CompositeGauss4"; // tested with other methods, OK too.

    /*--------get physical groups--------------*/
    std::map<std::string, std::pair<int, int>> groups;
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);
    for (std::size_t i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    /*--------get number of nodes-----------------*/
    std::vector<std::size_t> all_nodeTags;
    std::vector<double> all_nodecoord;
    std::vector<double> all_nodeparametricCoord;

    gmsh::model::mesh::getNodes(all_nodeTags, all_nodecoord, all_nodeparametricCoord);

    /*-------map between nodeTag and index of first nodal displacement associated to the node (second one is the same +1)------*/
    std::map<int, int> node_index;
    for (std::size_t i = 0; i < all_nodeTags.size(); ++i)
    {
        int nodeTag = all_nodeTags[i];// use i à la place de nodeTag pour traiter cas ou le nodeTag n'est pas une suite continue
        int node_first_index = 2*i; 
        node_index[nodeTag] = node_first_index;
        //std::cout << node_index[nodeTag] << " ";
    }

    /*-----------get physical properties of domain----------------------*/
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Volumic Forces|Materials).+");
    double E = 0; // Young modulus
    double nu = 0; // Poisson ratio
    double rho = 0; // volumic mass
    double bx = 0; // volumic force along x axis
    double by = 0; // volumic force along y axis

    Get_physical_Prop_Domain(keys, E, nu, rho, bx, by);

    /*-----------------H matrix involved in stress computation------------------*/
    Eigen::Matrix<double, 3, 3> H_matrix = Compute_H_matrix(E, nu);

    /*------------computing K matrix and f vector------------------------*/
    // we first consider the K matrix and the volumic part of the f vector
    Eigen::SparseMatrix<double> full_K_matrix(2*all_nodeTags.size(),2*all_nodeTags.size());
    std::vector<double> full_f_vector(2*all_nodeTags.size());

    std::vector<Eigen::Triplet <double>> k_tripletList;

    std::string groupname = "domain"; // we consider the domain as containing the whole 2d model

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        std::cerr << "Group '" << groupname << "' does not exist!\n";
        return 1;
    }

    // Getting entities of the domain
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

    // Looping over entities and compute K and F
    Loop_K_F(H_matrix, rho, bx, by, integration_rule, dim, tag, tags, full_K_matrix, 
              full_f_vector, k_tripletList);

    full_K_matrix.setFromTriplets(k_tripletList.begin(), k_tripletList.end());


    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    gmsh::onelab::getNames(keys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    Neumann_BC(keys, integration_rule, dim, tag, tags, node_index, groups, groupname, full_f_vector);

    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    // new_indices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not
    // if the given index has been removed, it will be set to -1 in new_indices
    std::vector<int> new_indices(2*all_nodeTags.size());
    for (std::size_t j = 0; j < new_indices.size(); ++j){
        new_indices[j] = (int) j;
    } 

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    std::vector<double> dir_BC(2*all_nodeTags.size());
    for (std::size_t j = 0; j < dir_BC.size(); ++j){
        dir_BC[j] = 0;
    }

    // Fill the two vectors defined just above
    Fill_New_indices_Dir_BC(keys, dim, tag, tags, node_index, groups, groupname, new_indices, dir_BC);

    int number_removed_lines = 0;
    for (std::size_t j = 0; j < new_indices.size(); ++j){
        if(new_indices[j] < 0){ // line has to be removed
            number_removed_lines++;
        }
        else{
            new_indices[j] -= number_removed_lines;
        }
    }


    /*------------FILLING THE REDUCED K MATRIX AND F VECTOR--------------*/
    //the correction term on the f vector comes from the "FEM" theoretical course (JP Ponthot), Chap2 Slide 38/50
    Eigen::SparseMatrix<double> reduced_K_matrix(2*all_nodeTags.size() - number_removed_lines,2*all_nodeTags.size() - number_removed_lines);
    std::vector<double> reduced_f_vector(2*all_nodeTags.size() - number_removed_lines);

    std::vector<Eigen::Triplet <double>> reduced_k_tripletList;

    Fill_K_F_Reduced(full_f_vector, full_K_matrix, new_indices, dir_BC, reduced_K_matrix, reduced_f_vector, reduced_k_tripletList);

    /*--------------SOLVING THE REDUCED SYSTEM--------------------*/
    Eigen::VectorXd final_reduced_d(2*all_nodeTags.size() - number_removed_lines);

    Eigen::VectorXd final_reduced_f(2*all_nodeTags.size() - number_removed_lines);
    for (std::size_t j = 0; j < reduced_f_vector.size(); ++j){
        final_reduced_f(j) = reduced_f_vector[j];
    }

    // solve Kd=f
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(reduced_K_matrix);

    final_reduced_d = solver.solve(final_reduced_f);


    /*-----------RECONSTRUCTING FULL NODAL VALUES--------------*/
    std::vector<double> full_d_vector(2*all_nodeTags.size());
    for (std::size_t j = 0; j < full_d_vector.size(); ++j){
        if(new_indices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = final_reduced_d(new_indices[j]);
        }
    }

    /*--------------COMPUTING NODAL FORCES-------------*/
    Eigen::VectorXd eigen_d(2*all_nodeTags.size());
    for(size_t i = 0 ; i < 2*all_nodeTags.size(); ++i){
        eigen_d(i) = full_d_vector[i];
    }
    Eigen::VectorXd nodal_forces(2*all_nodeTags.size());
    nodal_forces = full_K_matrix*eigen_d;

    /*--------------REACTION FORCES COMPUTATION------------------*/

    double Rx = 0.;
    double Ry = 0.;
    Reaction_Forces(Rx, Ry, dim, tag, tags, groupname, groups, node_index, nodal_forces, keys);

    /*--------------PLOTTING AND SAVING NODAL VALUES------------------*/
    // loop over the nodes
    std::vector<std::string> names = Plot_Save_ux_uy(all_nodeTags, full_d_vector);

    /*-------------POST PROCESSING---------------*/
    /*-------------Computing and saving strains/stresses-----------------*/

    // get all the elements for dim = 2
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, 2);
    
    // compute total number of 2D elements
    int nbelems=0;
    for(std::size_t i=0; i< elementTags.size(); ++i)
        nbelems+=elementTags[i].size();


    /*-------------------ELEMENTAL-NODAL VALUES-----------------------*/
    // éventuellement faire un tableau ici qui regroupe tout le bazar pour ne pas avoir la même ligne 7 fois ...
    int viewtag3 = gmsh::view::add("nodal_eps_x");
    int viewtag4 = gmsh::view::add("nodal_eps_y");
    int viewtag5 = gmsh::view::add("nodal_2eps_xy");
    int viewtag6 = gmsh::view::add("nodal_eps_z");
    int viewtag7 = gmsh::view::add("nodal_sig_x");
    int viewtag8 = gmsh::view::add("nodal_sig_y");
    int viewtag9 = gmsh::view::add("nodal_sig_xy");
    int viewtag10 = gmsh::view::add("nodal_sig_VM");
    std::vector<std::vector<double>> data3(nbelems);
    std::vector<std::vector<double>> data4(nbelems);
    std::vector<std::vector<double>> data5(nbelems);
    std::vector<std::vector<double>> data6(nbelems);
    std::vector<std::vector<double>> data7(nbelems);
    std::vector<std::vector<double>> data8(nbelems);
    std::vector<std::vector<double>> data9(nbelems);
    std::vector<std::vector<double>> data10(nbelems);
    std::vector<std::size_t> elems2D(nbelems);


    std::vector<int> nb_adjacent_elements(all_nodeTags.size()); // number of 2d elements adjacent to node (i+1)
    std::vector<double> nodal_sigma_x(all_nodeTags.size()); // will store nodal values (may be helpful later)
    std::vector<double> nodal_sigma_y(all_nodeTags.size());
    std::vector<double> nodal_sigma_xy(all_nodeTags.size());

    Compute_Stress_Strain_Elemental_Nodal(elementTypes, elementTags, elnodeTags, full_d_vector, elems2D, H_matrix, 
                                    nu, E, nb_adjacent_elements, nodal_sigma_x, nodal_sigma_y, nodal_sigma_xy, node_index,
                                    data3, data4, data5, data6, data7, data8, data9, data10);

    gmsh::view::addModelData(viewtag3, 0, names[0], "ElementNodeData",
                                elems2D, data3, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag4, 0, names[0], "ElementNodeData",
                                elems2D, data4, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag5, 0, names[0], "ElementNodeData",
                                elems2D, data5, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag6, 0, names[0], "ElementNodeData",
                                elems2D, data6, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag7, 0, names[0], "ElementNodeData",
                                elems2D, data7, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag8, 0, names[0], "ElementNodeData",
                                elems2D, data8, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag9, 0, names[0], "ElementNodeData",
                                elems2D, data9, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag10, 0, names[0], "ElementNodeData",
                                elems2D, data10, 1);  // the last ,1 is important!

    // representing vector displacement field and nodal forces
    int viewtag11 = gmsh::view::add("u");
    int viewtag12 = gmsh::view::add("f");
    std::vector<std::vector<double>> data11(all_nodeTags.size());
    std::vector<std::vector<double>> data12(all_nodeTags.size());
    for (std::size_t i = 0; i < all_nodeTags.size(); ++i)
    {
        data11[i].resize(3);
        data12[i].resize(3);
        data11[i][0] = full_d_vector[2*i]; //ux
        data11[i][1] = full_d_vector[2*i+1]; //uy
        data11[i][2] = 0.; //uz
        data12[i][0] = nodal_forces(2*i); //fx
        data12[i][1] = nodal_forces(2*i+1); //fy
        data12[i][2] = 0.; //fz
    }

    gmsh::view::addModelData(viewtag11, 0, names[0], "NodeData",
                             all_nodeTags, data11); //independent of time here
    gmsh::view::addModelData(viewtag12, 0, names[0], "NodeData",
                             all_nodeTags, data12); //independent of time here



    // empeche que tout se plot l'un sur l'autre
    int nb_views = 12; //hardcoded
    for( int i = 0; i < nb_views; ++i){
        std::string my_string = "View[" + std::to_string(i) + "].Visible";
        gmsh::option::setNumber(my_string, 0);
    }


    
    gmsh::fltk::run(); // opens the gmsh window
    
    gmsh::finalize();
    return 0;
}