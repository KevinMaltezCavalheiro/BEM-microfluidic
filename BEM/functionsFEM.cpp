#ifdef _MSC_VER
#include <gmsh.h_cwrap>       // gmsh main header
#else
#include <gmsh.h>       // gmsh main header
#endif

#include <iostream>     // for cout
#include <map>
#include <sstream>  // for stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices
#include <Eigen/SparseCholesky> // solving sparse linear systems
#include <Eigen/Core>
#include <cmath>
#include <algorithm> // to sort and merge the FEM node vectors
#include <omp.h>

using namespace gmsh;
using namespace std;
using namespace Eigen;

//Compute the matrix H
Matrix<double, 3, 3> Compute_H_matrix(const double E, const double nu){
    Matrix<double, 3, 3> H_matrix;
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

//Compute the full matrix "K".
void Assembly_K_triplet(const int i, const int el, 
                          const vector<vector<size_t>> & nodeTags,
                          const int numNodes, map<int, int> & node_map,
                          list<Triplet <double>> & k_tripletList,
                          Matrix<double, Dynamic, Dynamic> & K_matrix,
                          vector<double> & f_vector){

    for (int j = 0; j < 2*numNodes; ++j){
        int first_node_index = (int) j/2;
        int first_node = nodeTags[i][numNodes*el+first_node_index];
        int my_row = node_map[first_node] + j % 2; // node list begins at 1, '-1' correction to begin at 0
        for(int l = 0; l <= j; ++l){ // to exploit symmetry
            int second_node_index = (int) l/2;
            int second_node = nodeTags[i][numNodes*el+second_node_index];
            int my_col = node_map[second_node] + l % 2;

            k_tripletList.push_back(Triplet<double>(my_row,my_col,K_matrix(j,l)));
            if(j != l){
                k_tripletList.push_back(Triplet<double>(my_col,my_row,K_matrix(j,l)));
            }
        }
    }
}

//Compute the full matrix "f" of one given thread
void Assembly_f_thread(const int i, const int el, 
                          const vector<vector<size_t>> & nodeTags,
                          const int numNodes, map<int, int> & node_map,
                          list<pair<int,double>> & thread_f_vector,
                          vector<double> & f_vector){

    for (int j = 0; j < 2*numNodes; ++j){
        int first_node_index = (int) j/2;
        int first_node = nodeTags[i][numNodes*el+first_node_index];
        int my_row = node_map[first_node] + j % 2; 
        thread_f_vector.push_back(pair<int, double>(my_row, f_vector[j]));
    }
}


void Compute_B_matrix(const int j, const int numNodes, const vector<double> & basisFunctionsGradient,
                      const vector<double> & determinants,
                      const Matrix<double, 2, 2> & jacobinvtrans,
                      const Matrix<double, 3, 3> & H_matrix,
                      Matrix<double, Dynamic, Dynamic> & B_matrix){

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
void Get_physical_Prop_Domain(const vector<string> & keys, double & E, double & nu, double & rho,
                              double & bx, double & by){
    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
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

void Loop_K_F(const Matrix<double, 3, 3> H_matrix, const double rho, const double bx, const double by,
              const string & integration_rule,
              int & dim, vector<int> & tags, map<int, int>  & node_map,
              vector<list<Triplet <double>>> & my_k_list,
              vector<list<pair<int,double>>> & thread_f_vector){

    for (size_t k = 0; k < tags.size(); ++k)
    {
        //cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types, enlever les déclarations des boucles pr aller plus vite (cfr Boman)
        for (size_t i = 0; i < elementTypes.size(); ++i)
        {

            // getting element properties
            string elementName;
            int dim, order, numNodes, numPrimaryNodes;
            vector<double> localNodeCoord;
            model::mesh::getElementProperties(elementTypes[i],
                                            elementName, dim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);
            
            // get Gauss points coordinates and weights
            vector<double> localCoords, weights;
            model::mesh::getIntegrationPoints(elementTypes[i],
                                                    integration_rule, // "Gauss2"
                                                    localCoords, weights);

            // get determinants and Jacobians
            vector<double> jacobians, determinants, coords;
            model::mesh::getJacobians(elementTypes[i],
                                            localCoords, jacobians, 
                                            determinants, coords);

            // values and derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            int numComponents, numOrientations;
            vector<double> basisFunctions;
            vector<double> basisFunctionsGradient;
            model::mesh::getBasisFunctions(elementTypes[i], localCoords,
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);

            model::mesh::getBasisFunctions(elementTypes[i], localCoords,
                                                 "GradLagrange", numComponents, 
                                                 basisFunctionsGradient, 
                                                 numOrientations);

            /*---------computing the K matrix (and f vector, only volumic part for now)------*/
            // looping over the elements (computing elemental K matrices)
            //double test1 = omp_get_wtime();
            //setNbThreads(1);
            #pragma omp parallel for
            #ifdef _MSC_VER
            // https://github.com/tesseract-ocr/tesseract/issues/3285
            // OpenMP 2.0 standard is enforced - we may use /openmp:llvm and keep size_t
            for (int el = 0; el < elementTags[i].size(); el++){
            #else
            for (size_t el = 0; el < elementTags[i].size(); el++){
            #endif
                // elemental K matrix and f vector initialization
                Matrix<double, Dynamic, Dynamic> K_matrix(2*numNodes,2*numNodes);
                vector<double> f_vector(2*numNodes);
                K_matrix.setZero();
                for (int j = 0; j < 2*numNodes; ++j)f_vector[j] = 0.;

                // looping over the Gauss points (formula of slide 17)
                for (size_t j = 0; j < weights.size(); ++j){
                    // converting jacobian to Eigen format
                    Map<Matrix<double, Dynamic, Dynamic> >
                        jacob(&jacobians[9*j+9*weights.size()*el], 3, 3);
                    // 3D to 2D (plane stress case)
                    Matrix<double, 2, 2> jacob2D;
                    jacob2D(0,0) = jacob(0,0);
                    jacob2D(1,0) = jacob(1,0);
                    jacob2D(0,1) = jacob(0,1);
                    jacob2D(1,1) = jacob(1,1);
                            
                    // computing the inverse
                    Matrix<double, 2, 2> jacobinv = jacob2D.inverse();

                    // computing the transpose of the inverse
                    Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();

                    // computing the B matrix (formula slide 19) and the f vector (also using gauss integration)
                    Matrix<double, Dynamic, Dynamic> B_matrix(3,2*numNodes);


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
                //Assembly_K_triplet_F(i, el, nodeTags, numNodes, k_tripletList, full_f_vector, K_matrix, f_vector);
                Assembly_K_triplet(i, el, nodeTags, numNodes, node_map, my_k_list[omp_get_thread_num()], K_matrix, f_vector);
                Assembly_f_thread(i, el, nodeTags, numNodes, node_map, thread_f_vector[omp_get_thread_num()], f_vector);
                //cout << omp_get_thread_num() << "\n";
            }
            //double test2 = omp_get_wtime() - test1;
            //cout << "loop time: " << test2 << "\n";
        }
    }
}

/*-----------computing the f vector components related to this Neumann BC---------*/
void Compute_Neumann_BC(const int numNodes, const int i, const vector<vector<size_t>> & elementTags,
                        const double tx, const double ty, map<int, int> & node_map,
                        vector<double> & f_vector,
                        const vector<double> & weights, const vector<double> & basisFunctions,
                        const vector<double> & determinants, const vector<vector<size_t>> & nodeTags,
                        vector<double> & full_f_vector){

    for (size_t el = 0; el < elementTags[i].size(); el++){
        for (int j = 0; j < 2*numNodes; ++j){
            f_vector[j] = 0.;
        }
    // looping over the Gauss points (formula of slide 17)
        for (size_t j = 0; j < weights.size(); ++j){
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
            int my_row = node_map[first_node] + j % 2;
            full_f_vector[my_row] += f_vector[j];
        }
    }    
}

/*-----------focus on neumann BC's (surface traction) into full_f_vector------------*/
void Neumann_BC(const vector<string> & keys, const string & integration_rule,
                map<int, int> & node_map,
                map<string, pair<int, int>> & groups,
                vector<double> & full_f_vector){
    double tx = 0; // surface traction along x axis
    double ty = 0; // surface traction along y axis
    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            string groupname = ""; // will remain empty only if tx == ty == 0
            if(words[2].compare("tx") == 0){
                tx = value[0];
                string new_key = "Boundary Conditions/" + words[1] + "/ty";
                onelab::getNumber(new_key, value);
                ty = value[0];
                if(tx || ty) // there is no need the consider Neumann BC if it is homogeneous
                    groupname = words[1];
            }
            //here, if we are dealing with neumann BC's, we should have retrieved a non-empty group_name
            if(groupname.compare("") != 0){
                //loop over elements of group_name
                int dim = groups[groupname].first;
                int tag = groups[groupname].second;
                if (tag == 0)
                {
                    cerr << "Group '" << groupname << "' does not exist!\n";
                    return;
                }
                vector<int> tags;
                model::getEntitiesForPhysicalGroup(dim, tag, tags);

                // loop over entities
                for (size_t k = 0; k < tags.size(); ++k)
                {
                    vector<int> elementTypes;
                    vector<vector<size_t>> elementTags;
                    vector<vector<size_t>> nodeTags;
                    model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (size_t i = 0; i < elementTypes.size(); ++i)
                    {
                        // getting element properties
                        string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        vector<double> localNodeCoord;
                        model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, dim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                    
                        // get Gauss points coordinates and weights
                        vector<double> localCoords, weights;
                        model::mesh::getIntegrationPoints(elementTypes[i],
                                                                integration_rule, //maybe use different integration method on boundary than in domain
                                                                localCoords, weights);

                        // get determinants and Jacobians (only determinants will be useful)
                        vector<double> jacobians, determinants, coords;
                        model::mesh::getJacobians(elementTypes[i],
                                                        localCoords, jacobians, 
                                                        determinants, coords,tags[k]); // ai dû rajouter mon tags[k] pr focus on this particular entity

                        int numComponents, numOrientations;
                        vector<double> basisFunctions;
                        model::mesh::getBasisFunctions(elementTypes[i], localCoords,
                                                            "Lagrange", numComponents, 
                                                            basisFunctions, 
                                                            numOrientations);

                        /*-----------computing the f vector components related to this Neumann BC---------*/
                        //looping over the elements
                        vector<double> f_vector(2*numNodes);
                        Compute_Neumann_BC(numNodes, i, elementTags, tx, ty, node_map,
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
void Fill_New_indices_Dir_BC(const vector<string> & keys,
                             map<int, int> & node_map,
                             map<string, pair<int, int>> & groups,
                             vector<int> & new_indices,
                             vector<double> & dir_BC){

    for (size_t j = 0; j < new_indices.size(); ++j){
        new_indices[j] = (int) j;
    } 
    for (size_t j = 0; j < dir_BC.size(); ++j){
        dir_BC[j] = 0;
    }

    double ux = 0; // imposed displacement along x axis
    double uy = 0; // imposed displacement along y axis
    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("ux") == 0 || words[2].compare("uy") == 0){
                if(words[2].compare("ux") == 0)
                    ux = value[0];
                if(words[2].compare("uy") == 0)
                    uy = value[0];
                string groupname = words[1];
                
                //loop over elements of group_name
                int dim = groups[groupname].first;
                int tag = groups[groupname].second;
                if (tag == 0)
                {
                    cerr << "Group '" << groupname << "' does not exist!\n";
                    return;
                }
                vector<int> tags;
                model::getEntitiesForPhysicalGroup(dim, tag, tags);

                // loop over entities
                for (size_t k = 0; k < tags.size(); ++k)
                {
                    vector<int> elementTypes;
                    vector<vector<size_t>> elementTags;
                    vector<vector<size_t>> nodeTags;
                    model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (size_t i = 0; i < elementTypes.size(); ++i)
                    {

                        // getting element properties
                        string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        vector<double> localNodeCoord;
                        model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, dim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                        //looping over the elements
                        for (size_t el = 0; el < elementTags[i].size(); el++){
                            for (int j = 0; j < 2*numNodes; ++j){
                                int first_node_index = (int) j/2;
                                int first_node = nodeTags[i][numNodes*el+first_node_index];
                                int my_row = node_map[first_node] + j % 2; // j % 2 == 1 corresponds to y nodal displacement
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

void Fill_K_F_Reduced(const vector<double> & full_f_vector, const SparseMatrix<double> & full_K_matrix,
                      const vector<int> & new_indices, const vector<double> & dir_BC,
                      SparseMatrix<double> & reduced_K_matrix,
                      vector<double> & reduced_f_vector,
                      list<Triplet <double>> & reduced_k_tripletList){
    for (size_t j = 0; j < reduced_f_vector.size(); ++j){
        reduced_f_vector[j] = 0.;
    }

    // looping over the non-zero elements of full_K_matrix
    int row;
    int new_row;
    int new_col;
    double k_value;
    for (int col=0; col<full_K_matrix.outerSize(); ++col){
        new_col = new_indices[col];
        for (SparseMatrix<double>::InnerIterator it(full_K_matrix,col); it; ++it){
            k_value = it.value();
            row = it.row();   // row index
            new_row = new_indices[row];
            if(new_row >= 0 && new_col >= 0){ // if the given element must not be removed
                if(new_row < new_col){ // taking advantage of the symmetry
                    reduced_k_tripletList.push_back(Triplet<double>(new_row,new_col,k_value));
                    reduced_k_tripletList.push_back(Triplet<double>(new_col,new_row,k_value));
                }
                if(new_row == new_col){ // diagonal element
                    reduced_k_tripletList.push_back(Triplet<double>(new_row,new_col,k_value));
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

void Reaction_Forces(double & Rx, double & Ry,
                     map<string, pair<int, int>> & groups,
                     map<int, int> & node_map, VectorXd & nodal_forces,
                     const vector<string> keys){

    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("ux") == 0 || words[2].compare("uy") == 0){
                Rx = 0;
                Ry = 0;
                string groupname = words[1];
                cout << "-------------------- " << groupname << " ------------------- \n";
                
                //loop over elements of group_name
                int dim = groups[groupname].first;
                int tag = groups[groupname].second;
                if (tag == 0)
                {
                    cerr << "Group '" << groupname << "' does not exist!\n";
                    return;
                }
                vector<int> tags;
                model::getEntitiesForPhysicalGroup(dim, tag, tags);

                // loop over entities
                for (size_t k = 0; k < tags.size(); ++k)
                {
                    vector<size_t> nodeTags;
                    vector<double> coord; //useless
                    vector<double> parametricCoord; //useless

                    // retrieving all the nodes in this particular entity
                    model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tags[k], true);
                    
                    
                    for (size_t j = 0; j < nodeTags.size(); ++j){
                        if(words[2].compare("ux") == 0){    // computing Rx
                            int my_new_col = node_map[nodeTags[j]]; // le noeud numéro node_index[nodeTags[j]] correspond à la ligne/colonne dans K, f
                            Rx += nodal_forces(my_new_col);
                        }
                        if(words[2].compare("uy") == 0){    // computing Ry
                            int my_new_col = node_map[nodeTags[j]] + 1;
                            Ry += nodal_forces(my_new_col);
                        }
                    }
                }
                cout << "Rx: " << Rx << " [N/m], Ry: " << Ry << " [N/m] \n";
            }
        }
    }
}

vector<string> Plot_Save_ux_uy(const vector<size_t> & FEM_nodeTags, const vector<double> & full_d_vector, const int nbViews){
    int viewtag1 = view::add("ux");
    int viewtag2 = view::add("uy");
    vector<vector<double>> data1(FEM_nodeTags.size());
    vector<vector<double>> data2(FEM_nodeTags.size());
    for (size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        data1[i].resize(1);
        data1[i][0] = full_d_vector[2*i]; //ux
        data2[i].resize(1);
        data2[i][0] = full_d_vector[2*i+1]; //uy
    }

    vector<string> names;
    model::list(names);
    view::addModelData(viewtag1, 0, names[0], "NodeData",
                            FEM_nodeTags, data1); //independent of time here
    view::addModelData(viewtag2, 0, names[0], "NodeData",
                            FEM_nodeTags, data2); //independent of time here

    // save the results to disk
    //view::write(viewtag1, "ux.msh");
    //view::write(viewtag2, "uy.msh");

    // set display options & view results
    option::setNumber("View[" + to_string(nbViews) + "].IntervalsType", 3);
    option::setNumber("View[" + to_string(nbViews) + "].AdaptVisualizationGrid", 1);
    option::setNumber("View[" + to_string(nbViews) + "].MaxRecursionLevel", 3);
    option::setNumber("View[" + to_string(nbViews) + "].TargetError", -0.0001);
    option::setNumber("View[" + to_string(nbViews +1) + "].IntervalsType", 3);
    option::setNumber("View[" + to_string(nbViews +1) + "].AdaptVisualizationGrid", 1);
    option::setNumber("View[" + to_string(nbViews +1) + "].MaxRecursionLevel", 3);
    option::setNumber("View[" + to_string(nbViews +1) + "].TargetError", -0.0001);
    return names;
}

void Compute_Stress_Strain_Elemental_Nodal(const vector<int> & elementTypes, const vector<vector<size_t>> &elementTags,
                           const vector<vector<size_t>> &elnodeTags,
                           const vector<double> &full_d_vector, vector<size_t> & elems2D,
                           const Matrix<double, 3, 3> &H_matrix, const double nu, const double E,
                           map<int, int> & node_map,
                           vector<vector<double>> & data3, vector<vector<double>> & data4,
                           vector<vector<double>> & data5, vector<vector<double>> & data6,
                           vector<vector<double>> & data7, vector<vector<double>> & data8,
                           vector<vector<double>> & data9, vector<vector<double>> & data10,
                           int & k){

    for (size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
    {
        // get info about this type of element
        string name;
        int dim, order, numNodes, numPrimaryNodes;
        vector<double> paramCoord2D;
        model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                                numNodes, paramCoord2D, numPrimaryNodes);
        //paramCoord contains the coordinates of the nodes in 2D
        // passage en 3D pour utiliser getJacobians
        vector<double> paramCoord(3*paramCoord2D.size()/2);
        for(size_t j = 0; j < paramCoord2D.size()/2; ++j){
            paramCoord[3*j] = paramCoord2D[2*j];
            paramCoord[3*j+1] = paramCoord2D[2*j+1];
            paramCoord[3*j+2] = 0.;
        }
        

        // select element/node tags for this type of element
        auto &etags = elementTags[ityp];
        auto &enods = elnodeTags[ityp];

        // will contain the d vector inside each element
        VectorXd nodal_displacement(2*numNodes); //use eigen vector to use Eigen matric product later

        // get determinants and Jacobians (only jacob useful here) AT NODES AND NOT AT GAUSS POINTS
        vector<double> jacobians, determinants, coords;
        model::mesh::getJacobians(elementTypes[ityp],
                                        paramCoord, jacobians, 
                                        determinants, coords);

        // derivatives of the shape functions at the nodes of the reference element
        // we are here working with derivatives in the reference space, same for all elements
        int numComponents, numOrientations;
        vector<double> basisFunctionsGradient;

        model::mesh::getBasisFunctions(elementTypes[ityp], paramCoord,
                                                "GradLagrange", numComponents, 
                                                basisFunctionsGradient, 
                                                numOrientations);

        Vector3d epsilon; // will contain local values of the strains
        Vector3d sigma; // will contain local values of the stresses

        // loop over elements of type "ityp"
        for (size_t i = 0; i < etags.size(); ++i)
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
                nodal_displacement(2*j) = full_d_vector[node_map[nodeTag]];
                nodal_displacement(2*j+1) = full_d_vector[node_map[nodeTag]+1];
            }
            // looping over the nodes
            for (int j = 0; j < numNodes; ++j){
                // converting jacobian to Eigen format
                Map<Matrix<double, Dynamic, Dynamic> >
                    jacob(&jacobians[9*j+9*numNodes*i], 3, 3);
                // 3D to 2D (plane stress case)
                Matrix<double, 2, 2> jacob2D;
                jacob2D(0,0) = jacob(0,0);
                jacob2D(1,0) = jacob(1,0);
                jacob2D(0,1) = jacob(0,1);
                jacob2D(1,1) = jacob(1,1);
                        
                // computing the inverse
                Matrix<double, 2, 2> jacobinv = jacob2D.inverse();

                // computing the transpose of the inverse
                Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();

                // computing the B matrix (formula slide 19) 
                Matrix<double, Dynamic, Dynamic> B_matrix(3,2*numNodes);

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
            }
            k++;
        }
    }
}
void applyElecPressure(const string & integration_rule, map<int, int> & node_map, map<string, pair<int, int>> & groups, map<int, double> &electrostaticPressure, vector<double> & full_f_vector){

    string groupname = "BEM_FEM_boundary";
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elnodeTags;

    map<int,pair<double,double>> nodeCoordsMap;
    //filling the map with : nodeTag -> (x,y)
    for(size_t i=0; i< tags.size(); ++i)
    {
        vector<size_t> entity_nodeTags;
        vector<double> entity_nodecoord;
        vector<double> entity_nodeparametricCoord;
        model::mesh::getNodes(entity_nodeTags, entity_nodecoord, entity_nodeparametricCoord, dim, tags[i], true);
        for ( size_t j=0; j < entity_nodeTags.size(); j++){
            nodeCoordsMap[entity_nodeTags[j]] = pair<double, double>(entity_nodecoord[3*j], entity_nodecoord[3*j+1]);
        }
    }

    for(size_t i=0; i< tags.size(); ++i)
    {
        model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(size_t j=0; j< elementTypes.size(); ++j)
        {

            // getting element properties
            string elementName;
            int elDim, order, numNodes, numPrimaryNodes;
            vector<double> localNodeCoord;
            model::mesh::getElementProperties(elementTypes[j],
                                                        elementName, elDim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                    
            // get Gauss points coordinates and weights
            vector<double> localCoords, weights;
            model::mesh::getIntegrationPoints(elementTypes[j],
                                                    integration_rule, //maybe use different integration method on boundary than in domain
                                                    localCoords, weights);

            // get determinants and Jacobians (only determinants will be useful)
            vector<double> jacobians, determinants, coords;
            model::mesh::getJacobians(elementTypes[j],
                                            localCoords, jacobians, 
                                            determinants, coords,tags[i]); // ai dû rajouter mon tags[k] pr focus on this particular entity

            int numComponents, numOrientations;
            vector<double> basisFunctions;
            model::mesh::getBasisFunctions(elementTypes[j], localCoords, 
                                                "Lagrange", numComponents, 
                                                basisFunctions, 
                                                numOrientations);

            /*-----------computing the f vector components related to this Neumann BC---------*/
            //looping over the elements
            vector<double> f_vector(2*numNodes);


            for(size_t k=0; k < elementTags[j].size(); ++k)
            {
                int first_nodeTag = elnodeTags[j][numNodes*k];
                int second_nodeTag = elnodeTags[j][numNodes*k+1];
                double nx = -(nodeCoordsMap[second_nodeTag].second - nodeCoordsMap[first_nodeTag].second);
                double ny = (nodeCoordsMap[second_nodeTag].first - nodeCoordsMap[first_nodeTag].first);
                nx = nx / sqrt(nx*nx+ny*ny); //normalisation
                ny = ny / sqrt(nx*nx+ny*ny);
                double p = electrostaticPressure[elementTags[j][k]];
                double tx = p*nx;
                double ty = p*ny;

                for (int n = 0; n < 2*numNodes; ++n){
                    f_vector[n] = 0.;
                }
                // looping over the Gauss points (formula of slide 17)
                for (size_t n = 0; n < weights.size(); ++n){
                    // looping over the number of shape functions
                    for (int l = 0; l < numNodes; ++l){ 
                        f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*n]*tx*determinants[n+weights.size()*k]*weights[n]; //pas oublier le dét !
                        f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*n]*ty*determinants[n+weights.size()*k]*weights[n]; // le weight non plus !
                    }
                }
                /*------------- ASSEMBLY PROCESS ------------*/
                for (int n = 0; n < 2*numNodes; ++n){
                    int first_node_index = (int) n/2;
                    int first_node = elnodeTags[j][numNodes*k+first_node_index];
                    int my_row = node_map[first_node] + n % 2;
                    full_f_vector[my_row] += f_vector[n];
                }
            }
        }
    }
}
