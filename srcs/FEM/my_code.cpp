#include <gmsh.h>       // gmsh main header
#include <iostream>     // for std::cout
#include <map>
#include <sstream>  // for std::stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices (QUESTION: OK to use it ?)
#include <Eigen/SparseCholesky> // solving sparse linear systems (QUESTION: OK to use it ?)
#include <cmath>


// This code is only a functional DRAFT, it returns correct results, but can be improved
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
    std::string integration_rule = "CompositeGauss2"; // tested with other methods, OK too.

    /*--------get physical groups--------------*/
    std::map<std::string, std::pair<int, int>> groups;
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
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

    /*-----------get physical properties of domain----------------------*/
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Volumic Forces|Materials).+");
    double E; // Young modulus
    double nu; // Poisson ratio
    double rho; // volumic mass
    double bx; // volumic force along x axis
    double by; // volumic force along y axis
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

    /*-----------------H matrix involved in stress computation------------------*/
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

    /*------------computing K matrix and f vector------------------------*/
    // we first consider the K matrix and the volumic part of the f vector
    Eigen::SparseMatrix<double> full_K_matrix(2*all_nodeTags.size(),2*all_nodeTags.size());
    std::vector<double> full_f_vector(2*all_nodeTags.size());

    std::string groupname = "domain"; // we consider the domain as containing the whole 2d model

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        std::cerr << "Group '" << groupname << "' does not exist!\n";
        return 1;
    }

    // getting entities of the domain
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

    std::cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        std::cout << tags[i] << " ";
    std::cout << '\n';

    // looping over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        std::cout << "Entity (" << dim << "," << tags[k] << "):\n";
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // looping over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            /*std::cout << "\telement type " << elementTypes[i] << ":\n";
            std::cout << "\t\telements: ";
            for (int j = 0; j < elementTags[i].size(); ++j)
                std::cout << elementTags[i][j] << " ";
            std::cout << '\n';
            std::cout << "\t\tnodes: ";
            for (int j = 0; j < nodeTags[i].size(); ++j)
                std::cout << nodeTags[i][j] << " ";
            std::cout << '\n';*/

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
            std::cout << "\t" << weights.size() << " integration points\n";

            std::cout << "\tweights: ";
            for (int j = 0; j < weights.size(); ++j)
                std::cout << weights[j] << " ";
            std::cout << '\n';

            std::cout << "\tlocalCoords (gmsh): ";
            for (int j = 0; j < localCoords.size(); ++j)
                std::cout << localCoords[j] << " ";
            std::cout << '\n';

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
            for (int el = 0; el < elementTags[i].size(); el++){
            //for (int el = 0; el < 1; el++){
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
                for (int j = 0; j < weights.size(); ++j){
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
                    for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
                        B_matrix(0,2*l) = jacobinvtrans(0,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(0,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
                        B_matrix(1,2*l) = 0.;
                        B_matrix(2,2*l) = jacobinvtrans(1,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(1,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
                        B_matrix(0,2*l+1) = 0.;
                        B_matrix(1,2*l+1) = B_matrix(2,2*l);
                        B_matrix(2,2*l+1) = B_matrix(0,2*l);
                        f_vector[2*l] = f_vector[2*l] + rho*basisFunctions[l+numNodes*j]*bx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
                        f_vector[2*l+1] = f_vector[2*l+1] + rho*basisFunctions[l+numNodes*j]*by*determinants[j+weights.size()*el]*weights[j]; // pas oublier le weights[j] non plus !
                    }
                    
                    K_matrix = K_matrix + B_matrix.transpose()*H_matrix*B_matrix*determinants[j+weights.size()*el]*weights[j];
                }
                /*std::cout << "\telemental f vector of element " << elementTags[i][el] << ": ";
                for (int j = 0; j < f_vector.size(); ++j)
                    std::cout << f_vector[j] << " ";
                std::cout << '\n';*/
                //std::cout << "\t elemental K matrix of element " << elementTags[i][el] << ":\n" << K_matrix << "\n";
                /*------------- ASSEMBLY PROCESS ------------*/
                // USE LIST OF TRIPLETS TO PUSH ONLY ONCE -> more rapid and easier to parallelize
                for (int j = 0; j < 2*numNodes; ++j){
                    int first_node_index = (int) j/2;
                    int first_node = nodeTags[i][numNodes*el+first_node_index];
                    int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
                    for(int l = 0; l <= j; ++l){ // to exploit symmetry
                        int second_node_index = (int) l/2;
                        int second_node = nodeTags[i][numNodes*el+second_node_index];
                        int my_col = 2*(second_node-1) + l % 2;

                        full_K_matrix.coeffRef(my_row,my_col) += K_matrix(j,l);
                        if(j != l)
                            full_K_matrix.coeffRef(my_col,my_row) += K_matrix(l,j); // use symmetry
                    }
                    full_f_vector[my_row] += f_vector[j];
                }
            }
        }
    }
    /*std::cout << "\t full K matrix:\n" << full_K_matrix << "\n";
    std::cout << "\tfull f vector : ";
    for (int j = 0; j < full_f_vector.size(); ++j)
        std::cout << full_f_vector[j] << " ";
    std::cout << '\n';*/
    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    gmsh::onelab::getNames(keys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    double tx; // surface traction along x axis
    double ty; // surface traction along y axis
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
                std::cout << "------------------------------------" << "\n";
                std::cout << groupname << ", tx: " << tx << ", ty: " << ty << "\n";
                //loop over elements of group_name
                dim = groups[groupname].first;
                tag = groups[groupname].second;
                if (tag == 0)
                {
                    std::cerr << "Group '" << groupname << "' does not exist!\n";
                    return 1;
                }
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

                std::cout << "Entities in group named '" << groupname << "': ";
                for (int i = 0; i < tags.size(); ++i)
                    std::cout << tags[i] << " ";
                std::cout << '\n';
                // loop over entities
                for (int k = 0; k < tags.size(); ++k)
                {
                    std::cout << "Entity (" << dim << "," << tags[k] << "):\n";
                    std::vector<int> elementTypes;
                    std::vector<std::vector<std::size_t>> elementTags;
                    std::vector<std::vector<std::size_t>> nodeTags;
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (int i = 0; i < elementTypes.size(); ++i)
                    {
                        std::cout << "\telement type " << elementTypes[i] << ":\n";
                        std::cout << "\t\telements: ";
                        for (int j = 0; j < elementTags[i].size(); ++j)
                            std::cout << elementTags[i][j] << " ";
                        std::cout << '\n';
                        std::cout << "\t\tnodes: ";
                        for (int j = 0; j < nodeTags[i].size(); ++j)
                            std::cout << nodeTags[i][j] << " ";
                        std::cout << '\n';

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
                        std::cout << "\t" << weights.size() << " integration points\n";

                        std::cout << "\tweights: ";
                        for (int j = 0; j < weights.size(); ++j)
                            std::cout << weights[j] << " ";
                        std::cout << '\n';

                        std::cout << "\tlocalCoords (gmsh): ";
                        for (int j = 0; j < localCoords.size(); ++j)
                            std::cout << localCoords[j] << " ";
                        std::cout << '\n';

                        // get determinants and Jacobians (only determinants will be useful)
                        std::vector<double> jacobians, determinants, coords;
                        gmsh::model::mesh::getJacobians(elementTypes[i], 
                                                        localCoords, jacobians, 
                                                        determinants, coords,tags[k]); // ai dû rajouter mon tags[k] pr focus on this particular entity
                        std::cout << "\tGot " << determinants.size() 
                                << " Jacobians for type " << elementTypes[i] << "\n";

                        std::cout << "\tdeterminants (gmsh): ";
                        for (int j = 0; j < determinants.size(); ++j)
                            std::cout << determinants[j] << " ";
                        std::cout << '\n';
                        int numComponents, numOrientations;
                        std::vector<double> basisFunctions;
                        gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                            "Lagrange", numComponents, 
                                                            basisFunctions, 
                                                            numOrientations);
                        std::cout << "\tGot " << basisFunctions.size() 
                                << " basis functions for type "
                                << elementTypes[i] << '\n';
                        for (int j = 0; j < basisFunctions.size(); ++j)
                            std::cout << basisFunctions[j] << " ";
                        std::cout << '\n';
                        /*-----------computing the f vector components related to this Neumann BC---------*/
                        //looping over the elements
                        std::vector<double> f_vector(2*numNodes);
                        for (int el = 0; el < elementTags[i].size(); el++){
                            for (int j = 0; j < 2*numNodes; ++j){
                                f_vector[j] = 0.;
                            }
                            // looping over the Gauss points (formula of slide 17)
                            for (int j = 0; j < weights.size(); ++j){
                                // looping over the number of shape functions
                                for (int l = 0; l < numNodes; ++l){ 
                                    f_vector[2*l] = f_vector[2*l] + basisFunctions[l+numNodes*j]*tx*determinants[j+weights.size()*el]*weights[j]; //pas oublier le dét !
                                    f_vector[2*l+1] = f_vector[2*l+1] + basisFunctions[l+numNodes*j]*ty*determinants[j+weights.size()*el]*weights[j]; // le weight non plus !
                                }
                            }
                            std::cout << "\telemental f vector of element " << elementTags[i][el] << ": ";
                            for (int j = 0; j < f_vector.size(); ++j)
                                std::cout << f_vector[j] << " ";
                            std::cout << '\n';
                            /*------------- ASSEMBLY PROCESS ------------*/
                            for (int j = 0; j < 2*numNodes; ++j){
                                int first_node_index = (int) j/2;
                                int first_node = nodeTags[i][numNodes*el+first_node_index];
                                int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
                                full_f_vector[my_row] += f_vector[j];
                                //std::cout << my_row << "\n";
                            }
                        }
                    }
                }
            }
        }
    }
    /*std::cout << "\tfull f vector : ";
    for (int j = 0; j < full_f_vector.size(); ++j)
        std::cout << full_f_vector[j] << " ";
    std::cout << '\n';*/
    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    // new_indices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not
    // if the given index has been removed, it will set to -1 in new_indices
    std::vector<int> new_indices(2*all_nodeTags.size());
    for (int j = 0; j < new_indices.size(); ++j){
        new_indices[j] = j;
    } 

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    std::vector<double> dir_BC(2*all_nodeTags.size());
    for (int j = 0; j < dir_BC.size(); ++j){
        dir_BC[j] = 0;
    } 

    // looping over the boundary conditions imposed in the .geo file, the purpose of the whole loop is
    // to fill the two vectors defined just above
    double ux; // imposed displacement along x axis
    double uy; // imposed displacement along y axis
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
                    return 1;
                }
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

                std::cout << "Entities in group named '" << groupname << "': ";
                for (int i = 0; i < tags.size(); ++i)
                    std::cout << tags[i] << " ";
                std::cout << '\n';
                // loop over entities
                for (int k = 0; k < tags.size(); ++k)
                {
                    std::cout << "Entity (" << dim << "," << tags[k] << "):\n";
                    std::vector<int> elementTypes;
                    std::vector<std::vector<std::size_t>> elementTags;
                    std::vector<std::vector<std::size_t>> nodeTags;
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (int i = 0; i < elementTypes.size(); ++i)
                    {
                        std::cout << "\telement type " << elementTypes[i] << ":\n";
                        std::cout << "\t\telements: ";
                        for (int j = 0; j < elementTags[i].size(); ++j)
                            std::cout << elementTags[i][j] << " ";
                        std::cout << '\n';
                        std::cout << "\t\tnodes: ";
                        for (int j = 0; j < nodeTags[i].size(); ++j)
                            std::cout << nodeTags[i][j] << " ";
                        std::cout << '\n';

                        // getting element properties
                        std::string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        std::vector<double> localNodeCoord;
                        gmsh::model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, dim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                        //looping over the elements
                        for (int el = 0; el < elementTags[i].size(); el++){
                            for (int j = 0; j < 2*numNodes; ++j){
                                int first_node_index = (int) j/2;
                                int first_node = nodeTags[i][numNodes*el+first_node_index];
                                int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
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

    /*std::cout << "\t\tdir_BC: ";
    for (int j = 0; j < dir_BC.size(); ++j)
        std::cout << dir_BC[j] << " ";
    std::cout << '\n';*/

    int number_removed_lines = 0;
    for (int j = 0; j < new_indices.size(); ++j){
        if(new_indices[j] < 0){ // line has to be removed
            number_removed_lines++;
        }
        else{
            new_indices[j] -= number_removed_lines;
        }
    }

    /*std::cout << "\t\tnew indices: ";
    for (int j = 0; j < new_indices.size(); ++j)
        std::cout << new_indices[j] << " ";
    std::cout << '\n';*/

    /*------------FILLING THE REDUCED K MATRIX AND F VECTOR--------------*/
    //the correction term on the f vector comes from the "FEM" theoretical course (JP Ponthot), Chap2 Slide 38/50
    Eigen::SparseMatrix<double> reduced_K_matrix(2*all_nodeTags.size() - number_removed_lines,2*all_nodeTags.size() - number_removed_lines);
    std::vector<double> reduced_f_vector(2*all_nodeTags.size() - number_removed_lines);

    for (int j = 0; j < reduced_f_vector.size(); ++j){
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
                    reduced_K_matrix.coeffRef(new_row,new_col) = k_value;
                    reduced_K_matrix.coeffRef(new_col,new_row) = k_value;
                }
                if(new_row == new_col){ // diagonal element
                    reduced_K_matrix.coeffRef(new_row,new_col) = k_value;
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
    //std::cout << "\t full K matrix:\n" << full_K_matrix << "\n";
    //std::cout << "\t reduced K matrix:\n" << reduced_K_matrix << "\n";
    /*std::cout << "\treduced f vector : ";
    for (int j = 0; j < reduced_f_vector.size(); ++j)
        std::cout << reduced_f_vector[j] << " ";
    std::cout << '\n';*/

    /*--------------SOLVING THE REDUCED SYSTEM--------------------*/
    Eigen::VectorXd final_reduced_d(2*all_nodeTags.size() - number_removed_lines);

    Eigen::VectorXd final_reduced_f(2*all_nodeTags.size() - number_removed_lines);
    for (int j = 0; j < reduced_f_vector.size(); ++j){
        final_reduced_f(j) = reduced_f_vector[j];
    }
    //std::cout << "final f vector: " << final_reduced_f << "\n";


    // solve Kd=f
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(reduced_K_matrix);
    /*if(solver.info()!=Success) {
    // decomposition failed
    return;
    }*/
    final_reduced_d = solver.solve(final_reduced_f);
    /*if(solver.info()!=Success) {
    // solving failed
    return;
    }*/

    //std::cout << "final d vector: " << final_reduced_d << "\n";

    /*-----------RECONSTRUCTING FULL NODAL VALUES--------------*/
    std::vector<double> full_d_vector(2*all_nodeTags.size());
    for (int j = 0; j < full_d_vector.size(); ++j){
        if(new_indices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = final_reduced_d(new_indices[j]);
        }
    }
    /*std::cout << "\t final d vector : ";
    for (int j = 0; j < full_d_vector.size(); ++j)
        std::cout << full_d_vector[j] << " ";
    std::cout << '\n';*/

    /*--------------PLOTTING AND SAVING THEM------------------*/
    // loop over the nodes
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

    /*-------------POST PROCESSING---------------*/
    /*-------------Computing and saving strains/stresses-----------------*/
    // ELEMENTAL STRAIN/STRESS COMPUTATION for now (mean value of values computed at GP's)

    // get all the elements for dim=2
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, 2);
    
    // compute total number of 2D elements
    int nbelems=0;
    for(std::size_t i=0; i< elementTags.size(); ++i)
        nbelems+=elementTags[i].size();
    std::cout << "the mesh has " << nbelems << " elements\n";

    int viewtag3 = gmsh::view::add("eps_x");
    int viewtag4 = gmsh::view::add("eps_y");
    int viewtag5 = gmsh::view::add("2eps_xy");
    int viewtag6 = gmsh::view::add("eps_z");
    int viewtag7 = gmsh::view::add("sigma_x");
    int viewtag8 = gmsh::view::add("sigma_y");
    int viewtag9 = gmsh::view::add("sigma_xy");
    int viewtag10 = gmsh::view::add("sigma_VM");


    std::vector<std::vector<double>> data3(nbelems);
    std::vector<std::vector<double>> data4(nbelems);
    std::vector<std::vector<double>> data5(nbelems);
    std::vector<std::vector<double>> data6(nbelems);
    std::vector<std::vector<double>> data7(nbelems);
    std::vector<std::vector<double>> data8(nbelems);
    std::vector<std::vector<double>> data9(nbelems);
    std::vector<std::vector<double>> data10(nbelems);
    std::vector<std::size_t> elems2D(nbelems);

    // for each element type
    int k=0;
    for (std::size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
    {
        // get info about this type of element
        std::string name;
        int dim, order, numNodes, numPrimaryNodes;
        std::vector<double> paramCoord;
        gmsh::model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                                numNodes, paramCoord, numPrimaryNodes);

        std::cout << "paramCoord: " << "\n";
        for(int j = 0; j < paramCoord.size(); ++j)
            std::cout << paramCoord[j] << " ";
        std::cout << "\n";

        // select element/node tags for this type of element
        auto &etags = elementTags[ityp];
        auto &enods = elnodeTags[ityp];
        std::cout << numNodes << "\n";

        // will contain the d vector inside each element
        Eigen::VectorXd nodal_displacement(2*numNodes); //use eigen vector to use Eigen matric product later

        // get Gauss points coordinates and weights (only coordinates useful here)
        std::vector<double> localCoords, weights;
        gmsh::model::mesh::getIntegrationPoints(elementTypes[ityp], 
                                                integration_rule, // "Gauss2", en vrai ne pas le hardcoder (on l'utilise plusieurs fois dans le code)
                                                localCoords, weights);

        std::cout << "localCoords: " << "\n";
        for(int j = 0; j < localCoords.size(); ++j)
            std::cout << localCoords[j] << " ";
        std::cout << "\n";

        // get determinants and Jacobians (only jacob useful here)
        std::vector<double> jacobians, determinants, coords;
        gmsh::model::mesh::getJacobians(elementTypes[ityp], 
                                        localCoords, jacobians, 
                                        determinants, coords);

        // derivatives of the shape functions at the Gauss points of the reference element
        // we are here working with derivatives in the reference space, same for all elements
        int numComponents, numOrientations;
        std::vector<double> basisFunctionsGradient;

        gmsh::model::mesh::getBasisFunctions(elementTypes[ityp], localCoords, 
                                                "GradLagrange", numComponents, 
                                                basisFunctionsGradient, 
                                                numOrientations);

        Eigen::Vector3d epsilon; // will contain mean elemental values of the strains
        Eigen::Vector3d sigma; // will contain mean elemental values of the stresses
        Eigen::Vector3d tmp_epsilon; // will store temporary values

        // loop over elements of type "ityp"
        for (std::size_t i = 0; i < etags.size(); ++i)
        {
            elems2D[k] = etags[i];

            data3[k].resize(1);
            data4[k].resize(1);
            data5[k].resize(1);
            data6[k].resize(1);
            data7[k].resize(1);
            data8[k].resize(1);
            data9[k].resize(1);
            data10[k].resize(1);
            epsilon(0) = 0.;
            epsilon(1) = 0.;
            epsilon(2) = 0.;
            sigma(0) = 0.;
            sigma(1) = 0.;
            sigma(2) = 0.;
            // loop over nodes of the element, filling the local nodal displacement vector "d"
            for (std::size_t j = 0; j < numNodes; ++j)
            {
                int nodeTag = enods[numNodes * i + j];
                nodal_displacement(2*j) = full_d_vector[2*(nodeTag-1)];
                nodal_displacement(2*j+1) = full_d_vector[2*(nodeTag-1)+1];
            }
            // looping over the Gauss points (formula of slide 17)
            for (int j = 0; j < weights.size(); ++j){
                // converting jacobian to Eigen format
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                    jacob(&jacobians[9*j+9*weights.size()*i], 3, 3);
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
                for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
                    B_matrix(0,2*l) = jacobinvtrans(0,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(0,1)*basisFunctionsGradient[3*l+numNodes*3*j+1]; //ici faudrait pas un weights.size() à la place de numNodes ?? pense pas, ça fonctionne
                    B_matrix(1,2*l) = 0.;
                    B_matrix(2,2*l) = jacobinvtrans(1,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(1,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
                    B_matrix(0,2*l+1) = 0.;
                    B_matrix(1,2*l+1) = B_matrix(2,2*l);
                    B_matrix(2,2*l+1) = B_matrix(0,2*l);
                }
                tmp_epsilon = B_matrix*nodal_displacement;
                epsilon += tmp_epsilon;
                sigma += H_matrix*tmp_epsilon;
            }
            epsilon = epsilon / weights.size(); // performing an arithmetic mean value
            sigma = sigma / weights.size();
            data3[k][0] = epsilon(0);
            data4[k][0] = epsilon(1);
            data5[k][0] = epsilon(2);
            data6[k][0] = -nu/E*(sigma(0)+sigma(1));// formula slide 11/20 elasticity
            data7[k][0] = sigma(0);
            data8[k][0] = sigma(1);
            data9[k][0] = sigma(2);
            data10[k][0] = sqrt(sigma(0)*sigma(0) + sigma(1)*sigma(1) - sigma(0)*sigma(1) + 3*sigma(2)*sigma(2)); // VM stress in plane stress configuration

            k++;
        }
    }
    gmsh::view::addModelData(viewtag3, 0, names[0], "ElementData",
                                elems2D, data3);
    gmsh::view::addModelData(viewtag4, 0, names[0], "ElementData",
                                elems2D, data4);
    gmsh::view::addModelData(viewtag5, 0, names[0], "ElementData",
                                elems2D, data5);
    gmsh::view::addModelData(viewtag6, 0, names[0], "ElementData",
                                elems2D, data6);
    gmsh::view::addModelData(viewtag7, 0, names[0], "ElementData",
                                elems2D, data7);
    gmsh::view::addModelData(viewtag8, 0, names[0], "ElementData",
                                elems2D, data8);
    gmsh::view::addModelData(viewtag9, 0, names[0], "ElementData",
                                elems2D, data9);
    gmsh::view::addModelData(viewtag10, 0, names[0], "ElementData",
                                elems2D, data10);

    /*-------------------ELEMENTAL-NODAL VALUES-----------------------*/
    int viewtag11 = gmsh::view::add("nodal_eps_x");
    int viewtag12 = gmsh::view::add("nodal_eps_y");
    int viewtag13 = gmsh::view::add("nodal_2eps_xy");
    int viewtag14 = gmsh::view::add("nodal_eps_z");
    int viewtag15 = gmsh::view::add("nodal_sig_x");
    int viewtag16 = gmsh::view::add("nodal_sig_y");
    int viewtag17 = gmsh::view::add("nodal_sig_xy");
    int viewtag18 = gmsh::view::add("nodal_sig_VM");
    std::vector<std::vector<double>> data11(nbelems);
    std::vector<std::vector<double>> data12(nbelems);
    std::vector<std::vector<double>> data13(nbelems);
    std::vector<std::vector<double>> data14(nbelems);
    std::vector<std::vector<double>> data15(nbelems);
    std::vector<std::vector<double>> data16(nbelems);
    std::vector<std::vector<double>> data17(nbelems);
    std::vector<std::vector<double>> data18(nbelems);
    std::vector<std::size_t> elems2D_bis(nbelems);

    std::vector<int> nb_adjacent_elements(all_nodeTags.size()); // number of 2d elements adjacent to node (i+1)
    std::vector<double> nodal_sigma_x(all_nodeTags.size()); // will store nodal values (may be helpful later)
    std::vector<double> nodal_sigma_y(all_nodeTags.size());
    std::vector<double> nodal_sigma_xy(all_nodeTags.size());

    for(int i = 0; i < all_nodeTags.size(); ++i){
        nb_adjacent_elements[i] = 0;
        nodal_sigma_x[i] = 0.;
        nodal_sigma_y[i] = 0.;
        nodal_sigma_xy[i] = 0.;
    }

    // for each element type
    k=0;
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
        for( int j = 0; j < paramCoord2D.size()/2; ++j){
            paramCoord[3*j] = paramCoord2D[2*j];
            paramCoord[3*j+1] = paramCoord2D[2*j+1];
            paramCoord[3*j+2] = 0.;
        }

        /*std::cout << "paramCoord: " << "\n";
        for(int j = 0; j < paramCoord.size(); ++j)
            std::cout << paramCoord[j] << " ";
        std::cout << "\n";*/
        

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

        /*std::cout << "jacobians: " << "\n";
        for(int j = 0; j < jacobians.size(); ++j)
            std::cout << jacobians[j] << " ";
        std::cout << "\n";*/

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
        Eigen::Vector3d tmp_epsilon; // will store temporary values

        // loop over elements of type "ityp"
        for (std::size_t i = 0; i < etags.size(); ++i)
        {
            elems2D_bis[k] = etags[i];

            data11[k].resize(numNodes);
            data12[k].resize(numNodes);
            data13[k].resize(numNodes);
            data14[k].resize(numNodes);
            data15[k].resize(numNodes);
            data16[k].resize(numNodes);
            data17[k].resize(numNodes);
            data18[k].resize(numNodes);

            // loop over nodes of the element, filling the local nodal displacement vector "d"
            for (std::size_t j = 0; j < numNodes; ++j)
            {
                int nodeTag = enods[numNodes * i + j];
                nodal_displacement(2*j) = full_d_vector[2*(nodeTag-1)];
                nodal_displacement(2*j+1) = full_d_vector[2*(nodeTag-1)+1];

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
                for (int l = 0; l < numNodes; ++l){ // looping over the number of shape functions
                    B_matrix(0,2*l) = jacobinvtrans(0,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(0,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
                    B_matrix(1,2*l) = 0.;
                    B_matrix(2,2*l) = jacobinvtrans(1,0)*basisFunctionsGradient[3*l+numNodes*3*j]+jacobinvtrans(1,1)*basisFunctionsGradient[3*l+numNodes*3*j+1];
                    B_matrix(0,2*l+1) = 0.;
                    B_matrix(1,2*l+1) = B_matrix(2,2*l);
                    B_matrix(2,2*l+1) = B_matrix(0,2*l);
                }
                epsilon = B_matrix*nodal_displacement;
                sigma = H_matrix*epsilon;

                data11[k][j] = epsilon(0);
                data12[k][j] = epsilon(1);
                data13[k][j] = epsilon(2);
                data14[k][j] = -nu/E*(sigma(0)+sigma(1));// formula slide 11/20 elasticity
                data15[k][j] = sigma(0);
                data16[k][j] = sigma(1);
                data17[k][j] = sigma(2);
                data18[k][j] = sqrt(sigma(0)*sigma(0) + sigma(1)*sigma(1) - sigma(0)*sigma(1) + 3*sigma(2)*sigma(2)); // VM stress in plane stress configuration

                nodal_sigma_x[nodeTag - 1] += sigma(0);
                nodal_sigma_y[nodeTag - 1] += sigma(1);
                nodal_sigma_xy[nodeTag - 1] += sigma(2);
            }

            k++;
        }
    }
    gmsh::view::addModelData(viewtag11, 0, names[0], "ElementNodeData",
                                elems2D_bis, data11, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag12, 0, names[0], "ElementNodeData",
                                elems2D_bis, data12, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag13, 0, names[0], "ElementNodeData",
                                elems2D_bis, data13, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag14, 0, names[0], "ElementNodeData",
                                elems2D_bis, data14, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag15, 0, names[0], "ElementNodeData",
                                elems2D_bis, data15, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag16, 0, names[0], "ElementNodeData",
                                elems2D_bis, data16, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag17, 0, names[0], "ElementNodeData",
                                elems2D_bis, data17, 1);  // the last ,1 is important!
    gmsh::view::addModelData(viewtag18, 0, names[0], "ElementNodeData",
                                elems2D_bis, data18, 1);  // the last ,1 is important!

    // empeche que tout se plot l'un sur l'autre
    int nb_views = 18; //hardcoded
    for( int i = 0; i < nb_views; ++i){
        std::string hello = "View[" + std::to_string(i) + "].Visible";
        gmsh::option::setNumber(hello, 0);
    }

    // computing mean nodal values
    for (int i = 0; i < all_nodeTags.size(); ++i){
        nodal_sigma_x[i] = nodal_sigma_x[i] / nb_adjacent_elements[i];
        nodal_sigma_y[i] = nodal_sigma_y[i] / nb_adjacent_elements[i];
        nodal_sigma_xy[i] = nodal_sigma_xy[i] / nb_adjacent_elements[i];
    }

    // computing reaction forces (if ux imposed, need to compute Rx (pareil pr y))
    // using INTEGRATION PAR TRAPEZE
    double Rx = 0.;
    double Ry = 0.;
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
                    return 1;
                }
                gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

                // loop over entities
                for (int k = 0; k < tags.size(); ++k)
                {
                    std::vector<int> elementTypes;
                    std::vector<std::vector<std::size_t>> elementTags;
                    std::vector<std::vector<std::size_t>> nodeTags;
                    gmsh::model::mesh::getElements(elementTypes, elementTags,
                                                nodeTags, dim, tags[k]);
                    
                    for (int i = 0; i < elementTypes.size(); ++i)
                    {

                        // getting element properties
                        std::string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        std::vector<double> localNodeCoord;
                        gmsh::model::mesh::getElementProperties(elementTypes[i],
                                                        elementName, dim, order,
                                                        numNodes, localNodeCoord, 
                                                        numPrimaryNodes);
                        // getting nodal coordinates
                        std::vector<std::size_t> my_nodeTags; // rien de nouveau % nodeTags[i]
                        std::vector<double> model_coord; // utile ici !! les coordonnées
                        std::vector<double> param_coord; // inutile
                        gmsh::model::mesh::getNodesByElementType(elementTypes[i], my_nodeTags, model_coord, param_coord, tags[k]);

                        //looping over the elements (they are 1D -> 2 nodes)
                        for (int el = 0; el < elementTags[i].size(); el++){
                            int first_node = nodeTags[i][numNodes*el];
                            int second_node = nodeTags[i][numNodes*el + 1];
                            
                            // retrieving the coordinates
                            double x1 = model_coord[3*numNodes*el];
                            double y1 = model_coord[3*numNodes*el + 1];
                            double x2 = model_coord[3*numNodes*el + 3];
                            double y2 = model_coord[3*numNodes*el + 4];
                            // computing the normal (how to deal with the signs ?????)
                            double dx = x1 - x2;
                            double dy = y1 - y2;
                            double length = sqrt(dx*dx+dy*dy);
                            double nx;
                            double ny;
                            if (dx == 0){
                                nx = 1;
                                ny = 0;
                            }
                            else if (dy == 0){
                                nx = 0;
                                ny = 1;
                            }
                            else{
                                double slope = dx/dy;
                                nx = 1/sqrt(slope*slope+1);
                                ny = slope/sqrt(slope*slope+1);
                            }

                            double mean_sigma_x = (nodal_sigma_x[first_node - 1] + nodal_sigma_x[second_node - 1])/2; //-1 because nodeTags begin at 1
                            double mean_sigma_y = (nodal_sigma_y[first_node - 1] + nodal_sigma_y[second_node - 1])/2;
                            double mean_sigma_xy = (nodal_sigma_x[first_node - 1] + nodal_sigma_xy[second_node - 1])/2;

                            if(words[2].compare("ux") == 0){
                                Rx += (mean_sigma_x*nx + mean_sigma_xy*ny)*length;
                            }
                            if(words[2].compare("uy") == 0){
                                Ry += (mean_sigma_xy*nx + mean_sigma_y*ny)*length;
                            }
                        }
                    }
                }
                std::cout << "Rx: " << Rx << " [N/m], Ry: " << Ry << " [N/m] \n";
            }
        }
    }

    gmsh::fltk::run(); // opens the gmsh window
    
    gmsh::finalize();
    return 0;
}
