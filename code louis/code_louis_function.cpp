#include <gmsh.h>       // gmsh main header
#include <iostream>     // for cout
#include <map>
#include <sstream>  // for stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices (QUESTION: OK to use it ?)
#include <Eigen/SparseCholesky> // solving sparse linear systems (QUESTION: OK to use it ?)
#include <cmath>
#include "code_louis_function.hpp"

using namespace gmsh;
using namespace std;
using namespace Eigen;

void initialisation(double &E,double &nu,double &rho,double &bx,double &by,map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys);

void initialise_H_matrix(Matrix<double, 3, 3> &H_matrix, double nu, double E);

int compute_K_matrix(vector<size_t> &all_nodeTags, map<string, pair<int, int>> &groups, string &integration_rule, double bx, double rho, double by, Matrix<double, 3, 3> &H_matrix, SparseMatrix<double> &full_K_matrix,vector<double> &full_f_vector, string &groupname);

int boundary_conditions(vector<string> &keys, string &groupname, double &tx, double &ty, map<string, pair<int, int>> &groups, vector<int> &tags, string &integration_rule, vector<double> &full_f_vector);

int code_simple(int argc, char **argv);

int code_louis_function(int argc, char **argv);

int dirichlet(vector<int> &new_indices, vector<double> &dir_BC, vector<string> &keys, double &ux, double &uy, string &groupname, map<string, pair<int, int>> &groups, vector<int> &tags, int &number_removed_lines);

void filling(vector<size_t> &all_nodeTags, int &number_removed_lines, vector<int> &new_indices, SparseMatrix<double> &full_K_matrix, vector<double> &full_f_vector, vector<double> &dir_BC,SparseMatrix<double> &reduced_K_matrix,vector<double> &reduced_f_vector);

void reduce_system_solver(vector<size_t> &all_nodeTags, int &number_removed_lines, vector<double> &reduced_f_vector, SparseMatrix<double> &reduced_K_matrix, VectorXd &final_reduced_d, VectorXd &final_reduced_f);

void nodal_values(vector<double> &full_d_vector, vector<int> &new_indices, vector<double> &dir_BC, VectorXd &final_reduced_d);

void plot_and_save(vector<size_t> &all_nodeTags, vector<double> &full_d_vector, vector<string> &names);

//void post_processing(string &integration_rule, vector<double> &full_d_vector, Matrix<double, 3, 3> &H_matrix, double &nu, double &E, vector<string> &names);

int code_louis_function(int argc, char **argv)
{
    /*--------------INIT----------------*/
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    initialize();

    open(argv[1]);
    model::mesh::generate(2);

    /*PARAMETERS*/

    //physical
    double E(0); // Young modulus
    double nu(0); // Poisson ratio
    double rho(0); // volumic mass
    double bx(0); // volumic force along x axis
    double by(0); // volumic force along y axis
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    double tx(0); // surface traction along x axis
    double ty(0); // surface traction along y axis
    double ux(0); // imposed displacement along x axis
    double uy(0); // imposed displacement along y axis

    //number of nodes
    vector<size_t> all_nodeTags;
    vector<double> all_nodecoord;
    vector<double> all_nodeparametricCoord;

    //physical group
    map<string, pair<int, int>> groups;
    vector<int> tags;
    int number_removed_lines = 0;
    vector<string> names;

    //physical properties of the domain
    vector<string> keys;
    string groupname = "domain"; // we consider the domain as containing the whole 2d model

    /*--------Integration rule used in the code----------*/
    // could be useful to pass it as an argument in the .geo file
    string integration_rule = "CompositeGauss4"; // tested with other methods, OK too.
    //H matrix
    Matrix<double, 3, 3> H_matrix;

    /*FONCTIONS*/
    //initialise tout le domaine physique
    initialisation(E,nu,rho,bx,by,groups,all_nodeTags,all_nodecoord,all_nodeparametricCoord,keys);

    /*-----------------H matrix involved in stress computation------------------*/
    initialise_H_matrix(H_matrix,nu,E);

    /*------------computing K matrix and f vector------------------------*/
    //K matrix
    SparseMatrix<double> full_K_matrix(2*all_nodeTags.size(),2*all_nodeTags.size());
    //f vector
    vector<double> full_f_vector(2*all_nodeTags.size());
    compute_K_matrix(all_nodeTags,groups,integration_rule,bx,rho,by,H_matrix,full_K_matrix,full_f_vector,groupname);

    /*cout << "\t full K matrix:\n" << full_K_matrix << "\n";
     cout << "\tfull f vector : ";
     for (int j = 0; j < full_f_vector.size(); ++j)
     cout << full_f_vector[j] << " ";
     cout << '\n';*/
    /*-----------------BOUNDARY CONDITIONS------------*/

    boundary_conditions(keys,groupname,tx,ty,groups,tags,integration_rule,full_f_vector);

    /*cout << "\tfull f vector : ";
     for (int j = 0; j < full_f_vector.size(); ++j)
     cout << full_f_vector[j] << " ";
     cout << '\n';
     */
    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    vector<int> new_indices(2*all_nodeTags.size());
    vector<double> dir_BC(2*all_nodeTags.size());
    dirichlet(new_indices, dir_BC, keys, ux, uy, groupname, groups, tags, number_removed_lines);


    /*------------FILLING THE REDUCED K MATRIX AND F VECTOR--------------*/
    SparseMatrix<double> reduced_K_matrix(2*all_nodeTags.size() - number_removed_lines,2*all_nodeTags.size() - number_removed_lines);
    vector<double> reduced_f_vector(2*all_nodeTags.size() - number_removed_lines);
    filling(all_nodeTags,number_removed_lines,new_indices,full_K_matrix,full_f_vector,dir_BC,reduced_K_matrix,reduced_f_vector);

    /*--------------SOLVING THE REDUCED SYSTEM--------------------*/

    VectorXd final_reduced_d(2*all_nodeTags.size() - number_removed_lines);
    VectorXd final_reduced_f(2*all_nodeTags.size() - number_removed_lines);
    reduce_system_solver(all_nodeTags, number_removed_lines, reduced_f_vector, reduced_K_matrix, final_reduced_d, final_reduced_f);

    /*-----------RECONSTRUCTING FULL NODAL VALUES--------------*/
    vector<double> full_d_vector(2*all_nodeTags.size());
    nodal_values(full_d_vector,new_indices,dir_BC,final_reduced_d);

    /*cout << "\t final d vector : ";
     for (int j = 0; j < full_d_vector.size(); ++j)
     cout << full_d_vector[j] << " ";
     cout << '\n';*/

    /*--------------PLOTTING AND SAVING THEM------------------*/

    plot_and_save(all_nodeTags,full_d_vector,names);

    /*-------------POST PROCESSING---------------*/
    /*-------------Computing and saving strains/stresses-----------------*/
    // ELEMENTAL STRAIN/STRESS COMPUTATION for now (mean value of values computed at GP's)

    post_processing(integration_rule, full_d_vector, H_matrix, nu, E, names);

    finalize();
    return 0;
}

void initialisation(double &E,double &nu,double &rho,double &bx,double &by,map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys)
{
    /*--------get physical groups--------------*/
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    /*--------get number of nodes-----------------*/

    model::mesh::getNodes(all_nodeTags, all_nodecoord, all_nodeparametricCoord);

    /*-----------get physical properties of domain----------------------*/
    onelab::getNames(keys, "(Volumic Forces|Materials).+");
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
        while(getline(ss, word, '/'))// read string until '/'
        {
            words.push_back(word);
        }
        if(words.size()==3)
        {
            if(words[2].compare("Young") == 0)
            {
                E = value[0];
            }
            if(words[2].compare("Poisson") == 0)
            {
                nu = value[0];
            }
            if(words[2].compare("rho") == 0)
            {
                rho = value[0];
            }
            if(words[2].compare("bx") == 0)
            {
                bx = value[0];
            }
            if(words[2].compare("by") == 0)
            {
                by = value[0];
            }
        }
    }
}

void initialise_H_matrix(Matrix<double, 3, 3> &H_matrix, double nu, double E)
{
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

}

int compute_K_matrix(vector<size_t> &all_nodeTags, map<string, pair<int, int>> &groups, string &integration_rule, double bx, double rho, double by, Matrix<double, 3, 3> &H_matrix, SparseMatrix<double> &full_K_matrix,vector<double> &full_f_vector, string &groupname)
{
    // we first consider the K matrix and the volumic part of the f vector

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return 1;
    }

    // getting entities of the domain
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // looping over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                 nodeTags, dim, tags[k]);

        // looping over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            /*cout << "\telement type " << elementTypes[i] << ":\n";
             cout << "\t\telements: ";
             for (int j = 0; j < elementTags[i].size(); ++j)
             cout << elementTags[i][j] << " ";
             cout << '\n';
             cout << "\t\tnodes: ";
             for (int j = 0; j < nodeTags[i].size(); ++j)
             cout << nodeTags[i][j] << " ";
             cout << '\n';*/

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
            cout << "\t" << weights.size() << " integration points\n";

            cout << "\tweights: ";
            for (int j = 0; j < weights.size(); ++j)
                cout << weights[j] << " ";
            cout << '\n';

            cout << "\tlocalCoords (gmsh): ";
            for (int j = 0; j < localCoords.size(); ++j)
                cout << localCoords[j] << " ";
            cout << '\n';

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
            for (int el = 0; el < elementTags[i].size(); el++){
                //for (int el = 0; el < 1; el++){
                // elemental K matrix and f vector initialization
                Matrix<double, Dynamic, Dynamic> K_matrix(2*numNodes,2*numNodes);
                vector<double> f_vector(2*numNodes);
                for (int j = 0; j < 2*numNodes; ++j)
                {
                    for (int l = 0; l < 2*numNodes; ++l)
                    {
                        K_matrix(j,l)=0.;
                    }
                    f_vector[j] = 0.;
                }
                // looping over the Gauss points (formula of slide 17)
                for (int j = 0; j < weights.size(); ++j)
                {
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
                    for (int l = 0; l < numNodes; ++l)
                    { // looping over the number of shape functions
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
                /*cout << "\telemental f vector of element " << elementTags[i][el] << ": ";
                 for (int j = 0; j < f_vector.size(); ++j)
                 cout << f_vector[j] << " ";
                 cout << '\n';*/
                //cout << "\t elemental K matrix of element " << elementTags[i][el] << ":\n" << K_matrix << "\n";
                /*------------- ASSEMBLY PROCESS ------------*/
                // USE LIST OF TRIPLETS TO PUSH ONLY ONCE -> more rapid and easier to parallelize
                for (int j = 0; j < 2*numNodes; ++j)
                {
                    int first_node_index = (int) j/2;
                    int first_node = nodeTags[i][numNodes*el+first_node_index];
                    int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
                    for(int l = 0; l <= j; ++l)
                    { // to exploit symmetry
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
    /*cout << "\t full K matrix:\n" << full_K_matrix << "\n";
     cout << "\tfull f vector : ";
     for (int j = 0; j < full_f_vector.size(); ++j)
     cout << full_f_vector[j] << " ";
     cout << '\n';*/
    return 0;
}

int boundary_conditions(vector<string> &keys, string &groupname, double &tx, double &ty, map<string, pair<int, int>> &groups, vector<int> &tags, string &integration_rule, vector<double> &full_f_vector)
{
    /*-----------get group_names related to BC's----------------------*/
    onelab::getNames(keys, "(Boundary Conditions).+");

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
        while(getline(ss, word, '/'))// read string until '/'
        {
            words.push_back(word);
        }
        if(words.size()==3)
        {
            groupname = ""; // will remain empty only if tx == ty == 0
            if(words[2].compare("tx") == 0)
            {
                tx = value[0];
                string new_key = "Boundary Conditions/" + words[1] + "/ty";
                onelab::getNumber(new_key, value);
                ty = value[0];
                if(tx || ty) // there is no need the consider Neumann BC if it is homogeneous
                    groupname = words[1];
            }
            //here, if we are dealing with neumann BC's, we should have retrieved a non-empty group_name
            if(groupname.compare("") != 0)
            {
                cout << "------------------------------------" << "\n";
                cout << groupname << ", tx: " << tx << ", ty: " << ty << "\n";
                //loop over elements of group_name
                int dim = groups[groupname].first;
                int tag = groups[groupname].second;
                if (tag == 0)
                {
                    cerr << "Group '" << groupname << "' does not exist!\n";
                    return 1;
                }
                model::getEntitiesForPhysicalGroup(dim, tag, tags);

                cout << "Entities in group named '" << groupname << "': ";
                for (int i = 0; i < tags.size(); ++i)
                    cout << tags[i] << " ";
                cout << '\n';
                // loop over entities
                for (int k = 0; k < tags.size(); ++k)
                {
                    cout << "Entity (" << dim << "," << tags[k] << "):\n";
                    vector<int> elementTypes;
                    vector<vector<size_t>> elementTags;
                    vector<vector<size_t>> nodeTags;
                    model::mesh::getElements(elementTypes, elementTags,
                                             nodeTags, dim, tags[k]);

                    for (int i = 0; i < elementTypes.size(); ++i)
                    {
                        cout << "\telement type " << elementTypes[i] << ":\n";
                        cout << "\t\telements: ";
                        for (int j = 0; j < elementTags[i].size(); ++j)
                            cout << elementTags[i][j] << " ";
                        cout << '\n';
                        cout << "\t\tnodes: ";
                        for (int j = 0; j < nodeTags[i].size(); ++j)
                            cout << nodeTags[i][j] << " ";
                        cout << '\n';

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
                        cout << "\t" << weights.size() << " integration points\n";

                        cout << "\tweights: ";
                        for (int j = 0; j < weights.size(); ++j)
                            cout << weights[j] << " ";
                        cout << '\n';

                        cout << "\tlocalCoords (gmsh): ";
                        for (int j = 0; j < localCoords.size(); ++j)
                            cout << localCoords[j] << " ";
                        cout << '\n';

                        // get determinants and Jacobians (only determinants will be useful)
                        vector<double> jacobians, determinants, coords;
                        model::mesh::getJacobians(elementTypes[i],
                                                  localCoords, jacobians,
                                                  determinants, coords,tags[k]); // ai dû rajouter mon tags[k] pr focus on this particular entity
                        cout << "\tGot " << determinants.size()
                        << " Jacobians for type " << elementTypes[i] << "\n";

                        cout << "\tdeterminants (gmsh): ";
                        for (int j = 0; j < determinants.size(); ++j)
                            cout << determinants[j] << " ";
                        cout << '\n';
                        int numComponents, numOrientations;
                        vector<double> basisFunctions;
                        model::mesh::getBasisFunctions(elementTypes[i], localCoords,
                                                       "Lagrange", numComponents,
                                                       basisFunctions,
                                                       numOrientations);
                        cout << "\tGot " << basisFunctions.size()
                        << " basis functions for type "
                        << elementTypes[i] << '\n';
                        for (int j = 0; j < basisFunctions.size(); ++j)
                            cout << basisFunctions[j] << " ";
                        cout << '\n';
                        /*-----------computing the f vector components related to this Neumann BC---------*/
                        //looping over the elements
                        vector<double> f_vector(2*numNodes);
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
                            cout << "\telemental f vector of element " << elementTags[i][el] << ": ";
                            for (int j = 0; j < f_vector.size(); ++j)
                                cout << f_vector[j] << " ";
                            cout << '\n';
                            /*------------- ASSEMBLY PROCESS ------------*/
                            for (int j = 0; j < 2*numNodes; ++j){
                                int first_node_index = (int) j/2;
                                int first_node = nodeTags[i][numNodes*el+first_node_index];
                                int my_row = 2*(first_node-1) + j % 2; // node list begins at 1, '-1' correction to begin at 0
                                full_f_vector[my_row] += f_vector[j];
                                //cout << my_row << "\n";
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int dirichlet(vector<int> &new_indices, vector<double> &dir_BC, vector<string> &keys, double &ux, double &uy, string &groupname, map<string, pair<int, int>> &groups, vector<int> &tags, int &number_removed_lines)
{
    // new_indices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not
    // if the given index has been removed, it will set to -1 in new_indices
    
    for (int j = 0; j < new_indices.size(); ++j)
    {
        new_indices[j] = j;
    }

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    for (int j = 0; j < dir_BC.size(); ++j)
    {
        dir_BC[j] = 0;
    }

    // looping over the boundary conditions imposed in the .geo file, the purpose of the whole loop is
    // to fill the two vectors defined just above
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
                groupname = words[1];

                //loop over elements of group_name
                int dim = groups[groupname].first;
                int tag = groups[groupname].second;
                if (tag == 0)
                {
                    cerr << "Group '" << groupname << "' does not exist!\n";
                    return 1;
                }
                model::getEntitiesForPhysicalGroup(dim, tag, tags);

                cout << "Entities in group named '" << groupname << "': ";
                for (int i = 0; i < tags.size(); ++i)
                    cout << tags[i] << " ";
                cout << '\n';
                // loop over entities
                for (int k = 0; k < tags.size(); ++k)
                {
                    cout << "Entity (" << dim << "," << tags[k] << "):\n";
                    vector<int> elementTypes;
                    vector<vector<size_t>> elementTags;
                    vector<vector<size_t>> nodeTags;
                    model::mesh::getElements(elementTypes, elementTags,
                                             nodeTags, dim, tags[k]);

                    for (int i = 0; i < elementTypes.size(); ++i)
                    {
                        cout << "\telement type " << elementTypes[i] << ":\n";
                        cout << "\t\telements: ";
                        for (int j = 0; j < elementTags[i].size(); ++j)
                            cout << elementTags[i][j] << " ";
                        cout << '\n';
                        cout << "\t\tnodes: ";
                        for (int j = 0; j < nodeTags[i].size(); ++j)
                            cout << nodeTags[i][j] << " ";
                        cout << '\n';

                        // getting element properties
                        string elementName;
                        int dim, order, numNodes, numPrimaryNodes;
                        vector<double> localNodeCoord;
                        model::mesh::getElementProperties(elementTypes[i],
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

    /*cout << "\t\tdir_BC: ";
     for (int j = 0; j < dir_BC.size(); ++j)
     cout << dir_BC[j] << " ";
     cout << '\n';*/

    for (int j = 0; j < new_indices.size(); ++j){
        if(new_indices[j] < 0)
        { // line has to be removed
            number_removed_lines++;
        }
        else{
            new_indices[j] -= number_removed_lines;
        }
    }

    return 0;

    /*cout << "\t\tnew indices: ";
     for (int j = 0; j < new_indices.size(); ++j)
     cout << new_indices[j] << " ";
     cout << '\n';*/
}

void filling(vector<size_t> &all_nodeTags, int &number_removed_lines, vector<int> &new_indices, SparseMatrix<double> &full_K_matrix, vector<double> &full_f_vector, vector<double> &dir_BC,SparseMatrix<double> &reduced_K_matrix,vector<double> &reduced_f_vector)
{
    //the correction term on the f vector comes from the "FEM" theoretical course (JP Ponthot), Chap2 Slide 38/50

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
        for (SparseMatrix<double>::InnerIterator it(full_K_matrix,col); it; ++it){
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
    //cout << "\t full K matrix:\n" << full_K_matrix << "\n";
    //cout << "\t reduced K matrix:\n" << reduced_K_matrix << "\n";
    cout << "\treduced f vector : ";
    for (int j = 0; j < reduced_f_vector.size(); ++j)
        cout << reduced_f_vector[j] << " ";
    cout << '\n';
}

void reduce_system_solver(vector<size_t> &all_nodeTags, int &number_removed_lines, vector<double> &reduced_f_vector, SparseMatrix<double> &reduced_K_matrix, VectorXd &final_reduced_d, VectorXd &final_reduced_f)
{
    for (int j = 0; j < reduced_f_vector.size(); ++j)
    {
        final_reduced_f(j) = reduced_f_vector[j];
    }
    //cout << "final f vector: " << final_reduced_f << "\n";


    // solve Kd=f
    SimplicialLDLT<SparseMatrix<double> > solver;
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

    //cout << "final d vector: " << final_reduced_d << "\n";
}

void nodal_values(vector<double> &full_d_vector, vector<int> &new_indices, vector<double> &dir_BC, VectorXd &final_reduced_d)
{
    for (int j = 0; j < full_d_vector.size(); ++j){
        if(new_indices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = final_reduced_d(new_indices[j]);
        }
    }
}

void plot_and_save(vector<size_t> &all_nodeTags, vector<double> &full_d_vector, vector<string> &names)
{
    // loop over the nodes
    int viewtag1 = view::add("ux");
    int viewtag2 = view::add("uy");
    vector<vector<double>> data1(all_nodeTags.size());
    vector<vector<double>> data2(all_nodeTags.size());
    for (size_t i = 0; i < all_nodeTags.size(); ++i)
    {
        data1[i].resize(1);
        data1[i][0] = full_d_vector[2*i]; //ux
        data2[i].resize(1);
        data2[i][0] = full_d_vector[2*i+1]; //uy
    }

    model::list(names);
    view::addModelData(viewtag1, 0, names[0], "NodeData",
                       all_nodeTags, data1); //independent of time here
    view::addModelData(viewtag2, 0, names[0], "NodeData",
                       all_nodeTags, data2); //independent of time here

    // save the results to disk
    view::write(viewtag1, "ux.msh");
    view::write(viewtag2, "uy.msh");

    // set display options & view results
    option::setNumber("View[0].IntervalsType", 3);
    option::setNumber("View[0].AdaptVisualizationGrid", 1);
    option::setNumber("View[0].MaxRecursionLevel", 3);
    option::setNumber("View[0].TargetError", -0.0001);
    option::setNumber("View[1].IntervalsType", 3);
    option::setNumber("View[1].AdaptVisualizationGrid", 1);
    option::setNumber("View[1].MaxRecursionLevel", 3);
    option::setNumber("View[1].TargetError", -0.0001);

}

void post_processing(string &integration_rule, vector<double> &full_d_vector, Matrix<double, 3, 3> &H_matrix, double &nu, double &E, vector<string> &names)
{
    // get all the elements for dim=2
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elnodeTags;
    model::mesh::getElements(elementTypes, elementTags, elnodeTags, 2);

    // compute total number of 2D elements
    int nbelems=0;
    for(size_t i=0; i< elementTags.size(); ++i)
    {
        nbelems+=elementTags[i].size();
    }
    cout << "the mesh has " << nbelems << " elements\n";

    int viewtag3 = view::add("eps_x");
    int viewtag4 = view::add("eps_y");
    int viewtag5 = view::add("2eps_xy");
    int viewtag6 = view::add("eps_z");
    int viewtag7 = view::add("sigma_x");
    int viewtag8 = view::add("sigma_y");
    int viewtag9 = view::add("sigma_xy");
    int viewtag10 = view::add("sigma_VM");


    vector<vector<double>> data3(nbelems);
    vector<vector<double>> data4(nbelems);
    vector<vector<double>> data5(nbelems);
    vector<vector<double>> data6(nbelems);
    vector<vector<double>> data7(nbelems);
    vector<vector<double>> data8(nbelems);
    vector<vector<double>> data9(nbelems);
    vector<vector<double>> data10(nbelems);
    vector<size_t> elems2D(nbelems);

    // for each element type
    int k=0;
    for (size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
    {
        // get info about this type of element
        string name;
        int dim, order, numNodes, numPrimaryNodes;
        vector<double> paramCoord;
        model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                          numNodes, paramCoord, numPrimaryNodes);

        cout << "paramCoord: " << "\n";
        for(int j = 0; j < paramCoord.size(); ++j)
        {
            cout << paramCoord[j] << " ";
        }
        cout << "\n";

        // select element/node tags for this type of element
        auto &etags = elementTags[ityp];
        auto &enods = elnodeTags[ityp];
        cout << numNodes << "\n";

        // will contain the d vector inside each element
        VectorXd nodal_displacement(2*numNodes); //use eigen vector to use Eigen matric product later

        // get Gauss points coordinates and weights (only coordinates useful here)
        vector<double> localCoords, weights;
        model::mesh::getIntegrationPoints(elementTypes[ityp], integration_rule, localCoords, weights);// "Gauss2", en vrai ne pas le hardcoder (on l'utilise plusieurs fois dans le code)

        cout << "localCoords: " << "\n";
        for(int j = 0; j < localCoords.size(); ++j)
            cout << localCoords[j] << " ";
        cout << "\n";


        // get determinants and Jacobians (only jacob useful here)
        vector<double> jacobians, determinants, coords;
        model::mesh::getJacobians(elementTypes[ityp],
                                  localCoords, jacobians,
                                  determinants, coords);

        // derivatives of the shape functions at the Gauss points of the reference element
        // we are here working with derivatives in the reference space, same for all elements
        int numComponents, numOrientations;
        vector<double> basisFunctionsGradient;

        model::mesh::getBasisFunctions(elementTypes[ityp], localCoords,
                                       "GradLagrange", numComponents,
                                       basisFunctionsGradient,
                                       numOrientations);

        Vector3d epsilon; // will contain mean elemental values of the strains
        Vector3d sigma; // will contain mean elemental values of the stresses
        Vector3d tmp_epsilon; // will store temporary values

        // loop over elements of type "ityp"
        for (size_t i = 0; i < etags.size(); ++i)
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
            for (size_t j = 0; j < numNodes; ++j)
            {
                int nodeTag = enods[numNodes * i + j];
                nodal_displacement(2*j) = full_d_vector[2*(nodeTag-1)];
                nodal_displacement(2*j+1) = full_d_vector[2*(nodeTag-1)+1];
            }
            // looping over the Gauss points (formula of slide 17)
            for (int j = 0; j < weights.size(); ++j)
            {
                // converting jacobian to Eigen format
                Map<Matrix<double, Dynamic, Dynamic> >
                jacob(&jacobians[9*j+9*weights.size()*i], 3, 3);
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
                for (int l = 0; l < numNodes; ++l)
                { // looping over the number of shape functions
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
    view::addModelData(viewtag3, 0, names[0], "ElementData", elems2D, data3);
    view::addModelData(viewtag4, 0, names[0], "ElementData", elems2D, data4);
    view::addModelData(viewtag5, 0, names[0], "ElementData", elems2D, data5);
    view::addModelData(viewtag6, 0, names[0], "ElementData", elems2D, data6);
    view::addModelData(viewtag7, 0, names[0], "ElementData", elems2D, data7);
    view::addModelData(viewtag8, 0, names[0], "ElementData", elems2D, data8);
    view::addModelData(viewtag9, 0, names[0], "ElementData", elems2D, data9);
    view::addModelData(viewtag10, 0, names[0], "ElementData", elems2D, data10);

    /*-------------------ELEMENTAL-NODAL VALUES-----------------------*/
        int viewtag11 = view::add("nodal_eps_x");
        int viewtag12 = view::add("nodal_eps_y");
        int viewtag13 = view::add("nodal_2eps_xy");
        int viewtag14 = view::add("nodal_eps_z");
        int viewtag15 = view::add("nodal_sig_x");
        int viewtag16 = view::add("nodal_sig_y");
        int viewtag17 = view::add("nodal_sig_xy");
        int viewtag18 = view::add("nodal_sig_VM");
        vector<vector<double>> data11(nbelems);
        vector<vector<double>> data12(nbelems);
        vector<vector<double>> data13(nbelems);
        vector<vector<double>> data14(nbelems);
        vector<vector<double>> data15(nbelems);
        vector<vector<double>> data16(nbelems);
        vector<vector<double>> data17(nbelems);
        vector<vector<double>> data18(nbelems);
        vector<size_t> elems2D_bis(nbelems);

        // for each element type
        k=0;
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
            for( int j = 0; j < paramCoord2D.size()/2; ++j){
                paramCoord[3*j] = paramCoord2D[2*j];
                paramCoord[3*j+1] = paramCoord2D[2*j+1];
                paramCoord[3*j+2] = 0.;
            }

            /*cout << "paramCoord: " << "\n";
            for(int j = 0; j < paramCoord.size(); ++j)
                cout << paramCoord[j] << " ";
            cout << "\n";*/


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

            /*cout << "jacobians: " << "\n";
            for(int j = 0; j < jacobians.size(); ++j)
                cout << jacobians[j] << " ";
            cout << "\n";*/

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
            Vector3d tmp_epsilon; // will store temporary values

            // loop over elements of type "ityp"
            for (size_t i = 0; i < etags.size(); ++i)
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
                for (size_t j = 0; j < numNodes; ++j)
                {
                    int nodeTag = enods[numNodes * i + j];
                    nodal_displacement(2*j) = full_d_vector[2*(nodeTag-1)];
                    nodal_displacement(2*j+1) = full_d_vector[2*(nodeTag-1)+1];
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
                }

                k++;
            }
        }
        cout << "hello";
        view::addModelData(viewtag11, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data11, 1);  // the last ,1 is important!
        view::addModelData(viewtag12, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data12, 1);  // the last ,1 is important!
        view::addModelData(viewtag13, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data13, 1);  // the last ,1 is important!
        view::addModelData(viewtag14, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data14, 1);  // the last ,1 is important!
        view::addModelData(viewtag15, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data15, 1);  // the last ,1 is important!
        view::addModelData(viewtag16, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data16, 1);  // the last ,1 is important!
        view::addModelData(viewtag17, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data17, 1);  // the last ,1 is important!
        view::addModelData(viewtag18, 0, names[0], "ElementNodeData",
                                    elems2D_bis, data18, 1);  // the last ,1 is important!


    fltk::run(); // opens the gmsh window
}

//int code_simple(int argc, char **argv)
//{
//    if (argc < 2)
//    {
//        cout << "Usage: " << argv[0] << " <geo_file>\n";
//        return 0;
//    }
//    initialize();
//
//    open(argv[1]);
//    model::mesh::generate(2);
//
//    vectorpair dim_tags;
//    model::getPhysicalGroups(dim_tags);
//    //model::getEntities(dim_tags);
//    cout << dim_tags.size() << endl;
//
//    for(int i = 0; i < dim_tags.size(); i++)
//    {
//        int dim = dim_tags[i].first;
//        int tag = dim_tags[i].second;
//        string name;
//        model::getPhysicalName(dim, tag, name);
//        cout << "Physical group (" << dim << "D) " << "named '" << name << "' tag = " << tag << endl;
//    }
//
//    fltk::run(); // opens the gmsh window
//
//    finalize();
//    return 0;
//}
