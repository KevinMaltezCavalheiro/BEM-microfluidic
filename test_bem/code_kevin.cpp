#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include "code_kevin.hpp"

using namespace gmsh;
using namespace std;
using namespace Eigen;

void show_mat(MatrixXd &A);

void show_vec(VectorXd &A);

int code_kevin(int argc, char **argv);

//todo: séparer en différentes fonctions
//todo: comparer au livre BEM
//todo: faire fonctionner le code :p
//todo: voir pour gérer des trous dans les surfaces
//todo: combiner avec la FEM

int code_kevin(int argc, char **argv)
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
    double neumann(0); // value of one neumann boundary condition
    double dirichlet(0); // value of one dirichlet boundary condition
    vector<int> index; // contain the value 0 if the element has a neumann boundary condition and 1 if it is a dirichlet
    VectorXd boundary_conditions_values(0); // contain the values relative the boundary conditions. The index vector tells if it is a neumann or a dirichlet
    VectorXd sol_neumann(0);// contain the final solution for the boundary values of du/dn.
    VectorXd sol_dirichlet(0);// contain the final solution for the boundary values of u.
    int N = 0;//total number of elements


    //vectors to contain all the nodes information
    vector<size_t> all_nodeTags;
    vector<double> all_nodecoord;
    vector<double> all_nodeparametricCoord;

    //physical group
    map<string, pair<int, int>> groups;
    vector<string> names;
    vectorpair dimTags;

    //physical properties of the boundary
    vector<string> keys;
    string groupname = "boundary"; // on considère la frontière du domaine (ce sera principalement pour plus tard quand on combinera avec la FEM)

    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    model::getPhysicalGroups(dimTags,1);// 1 car dans notre cas la frontière est de dimension 1 (on n'a pas besoin de la surface)
    model::mesh::getNodes(all_nodeTags, all_nodecoord, all_nodeparametricCoord);
    N = all_nodeTags.size();//number of nodes = number of elements
    index.resize(N,-1);//set all the values at -1 to be sure it is fill with 0 or 1.
    boundary_conditions_values.resize(N);

    //fill the map groups to gather all the information about the geometry
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    //read the .geo file to obtain all information about the boundary conditions.
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
            groupname = "";
            if(words[2].compare("neumann") == 0)
            {
                neumann = value[0];
                groupname = words[1];
            }
            else if(words[2].compare("dirichlet") == 0)
            {
                dirichlet = value[0];
                groupname = words[1];
            }
            if(groupname.compare("") != 0)
            {
                cout << "------------------------------------" << "\n";
                if(words[2].compare("dirichlet") == 0)
                {
                    cout << groupname << ", dirichlet: " << dirichlet << endl;
                }
                if(words[2].compare("neumann") == 0)
                {
                    cout << groupname << ", neumann: " << neumann << endl;
                }
                //loop over elements of group_name
                int dim = groups[groupname].first;
                int tag = groups[groupname].second;
                if (tag == 0)
                {
                    cerr << "Group '" << groupname << "' does not exist!\n";
                    return 1;
                }
                vector<int> tags;
                model::getEntitiesForPhysicalGroup(dim, tag, tags);

                cout << "Entities in group named '" << groupname << "': ";
                for (int i = 0; i < tags.size(); ++i)
                {
                    cout << tags[i] << " ";
                }
                cout << endl;
                // loop over entities
                for (int k = 0; k < tags.size(); ++k)
                {
                    //cout << "Entity (" << dim << "," << tags[k] << "):\n";
                    vector<int> elementTypes;
                    vector<size_t> elementTags;
                    vector<size_t> nodeTags;
                    model::mesh::getElementsByType(1,elementTags,nodeTags,tags[k]);
                    for (int l = 0; l < elementTags.size(); l++)
                    {
                        if(words[2].compare("neumann") == 0)
                        {
                            //if the element has a neumann boundary condition, we store it this information by puting a 0 in the index vector and we store the value of this BC in the vector boundary_conditions_values
                            index[elementTags[l]-1] = 0;
                            boundary_conditions_values(elementTags[l]-1) = neumann;
                        }
                        else if(words[2].compare("dirichlet") == 0)
                        {
                            //if the element has a dirichlet boundary condition, we store it this information by puting a 1 in the index vector and we store the value of this BC in the vector boundary_conditions_values
                            index[elementTags[l]-1] = 1;
                            boundary_conditions_values(elementTags[l]-1) = dirichlet;
                        }
                    }
                }
            }
        }
    }

    /*--------Integration rule used in the code----------*/
    string integration_rule = "CompositeGauss4";

    //defines all matrices related to the BEM method
    MatrixXd H_matrix;
    MatrixXd G_matrix;
    MatrixXd A_matrix;
    MatrixXd B_matrix;
    VectorXd b_vector;

    /*-----------------compute the matrix G and H------------------*/

    //resize all the matrices and initilize them to 0.
    G_matrix.resize(N,N);
    H_matrix.resize(N,N);
    A_matrix.resize(N,N);
    B_matrix.resize(N,N);
    b_vector.resize(N);
    G_matrix.array() = 0.;
    H_matrix.array() = 0.;
    A_matrix.array() = 0.;
    B_matrix.array() = 0.;

    // first loop over the elements: indices i -> k -> l
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim_1 = dimTags[i].first;
        int tag_1 = dimTags[i].second;
        string name_1;
        model::getPhysicalName(dim_1, tag_1, name_1);
        //si mon name c'est dirichlet, on peut trier les conditions de bords et les mettre dans dirichlet par après on pourra parcourir la matrice et trier G et H en fonctions de ces dernières
        groups[name_1] = {dim_1, tag_1};
        vector<int> tags_1;
        model::getEntitiesForPhysicalGroup(dim_1, tag_1, tags_1);

        //check qu'il n'y a pas de trous dans la numérotation des noeuds
        //on suppose que all_nodeTags n'a pas de trou, est continu et commence à 1
        //check via une petite routine
        //model::mesh::renumberNodes();
        for(int k = 0; k < tags_1.size(); k++)
        {
            vector<size_t> elementTags_1, nodeTags_1;
            model::mesh::getElementsByType(1,elementTags_1,nodeTags_1,tags_1[k]);// on a que des segments droits
            //cout << elementTags.size() << " " << nodeTags.size() << endl;
            for(int l = 0; l < elementTags_1.size(); l++)
            {
                unsigned long element_l = elementTags_1[l]-1;//to avoid the -1 between the tag notation and the c++
                unsigned long node_1_tag_l = nodeTags_1[2*l];//first node of element l
                unsigned long node_2_tag_l = nodeTags_1[2*l+1];//second node of element l
                double node_1_coord_x_l = all_nodecoord[3*(node_1_tag_l-1)];//x coordinate of the first node of element l
                double node_1_coord_y_l = all_nodecoord[3*(node_1_tag_l-1)+1];//y coordinate of the first node of element l
                double node_2_coord_x_l = all_nodecoord[3*(node_2_tag_l-1)];//x coordinate of the second node of element l
                double node_2_coord_y_l = all_nodecoord[3*(node_2_tag_l-1)+1];//y coordinate of the second node of element l

                /*-------affiche les différents tags et les coordonnées de noeuds-------*/
                //                cout << "element " << elementTags[l] << " nodes = " << nodeTags[2*l] << "," << nodeTags[2*l+1] <<
                //                 " x1= " << all_nodecoord[3*(nodeTags[2*l]-1)] << " ,y1= " << all_nodecoord[3*(nodeTags[2*l]-1)+1] <<
                //                " x2= " << all_nodecoord[3*(nodeTags[2*l+1]-1)] << " ,y2= " << all_nodecoord[3*(nodeTags[2*l+1]-1)+1]<<               endl;
                //cout << "element " << element << " nodes = " << node_1_tag << "," << node_2_tag << " x1= " << node_1_coord_x << " ,y1= " << node_1_coord_y << " x2= " << node_2_coord_x << " ,y2= " << node_2_coord_y << endl;

                /*-------compute the node at the middle of the element-------*/
                double middle_node_element_x_l = (node_1_coord_x_l+node_2_coord_x_l)/2.0;
                double middle_node_element_y_l = (node_1_coord_y_l+node_2_coord_y_l)/2.0;
                //cout << "middle node element: (" << middle_node_element_x << "," << middle_node_element_y << ")" << endl;

                // second loop over the elements: indices j -> n -> m
                for(int j = 0; j < dimTags.size(); j++)
                {
                    int dim_2 = dimTags[j].first;
                    int tag_2 = dimTags[j].second;
                    string name_2;
                    model::getPhysicalName(dim_2, tag_2, name_2);
                    groups[name_2] = {dim_2, tag_2};
                    vector<int> tags_2;
                    model::getEntitiesForPhysicalGroup(dim_2, tag_2, tags_2);
                    for(int n = 0; n < tags_2.size(); n++)
                    {
                        vector<size_t> elementTags_2, nodeTags_2;
                        model::mesh::getElementsByType(1,elementTags_2,nodeTags_2,tags_2[n]);
                        for(int m = 0; m < elementTags_2.size(); m++)
                        {
                            unsigned long element_m = elementTags_2[m]-1;//to avoid the -1 between the tag notation and the c++

                            //just to display the elements for debugging
                            //cout << element_l << " " << element_m << endl;

                            //condition i == j in the slides, to compute Hii and Gii
                            if(element_l == element_m)
                            {
                                //computation of Gii
                                double length_element_l = sqrt((node_1_coord_x_l-node_2_coord_x_l)*(node_1_coord_x_l-node_2_coord_x_l)+(node_1_coord_y_l-node_2_coord_y_l)*(node_1_coord_y_l-node_2_coord_y_l)); //equivalent to li in the slides
                                double g_ii = length_element_l/(2.0*M_PI)*(log(length_element_l/2.0)-1);
                                G_matrix(element_l,element_m) = g_ii;

                                //computation of Hii
                                double h_ii = -1.0/2.0;
                                H_matrix(element_l,element_m) = h_ii;
                            }
                            // condition i != j in the slides
                            else if(element_l != element_m)
                            {
                                //I have defined temporary variables to help debug the code.
                                unsigned long node_1_tag_m = nodeTags_2[2*m];
                                unsigned long node_2_tag_m = nodeTags_2[2*m+1];
                                double node_1_coord_x_m = all_nodecoord[3*(node_1_tag_m-1)];
                                double node_1_coord_y_m = all_nodecoord[3*(node_1_tag_m-1)+1];
                                double node_2_coord_x_m = all_nodecoord[3*(node_2_tag_m-1)];
                                double node_2_coord_y_m = all_nodecoord[3*(node_2_tag_m-1)+1];
                                double temp1 = (node_2_coord_y_m - middle_node_element_y_l)/(node_2_coord_x_m - middle_node_element_x_l);
                                double temp2 = (node_1_coord_y_m - middle_node_element_y_l)/(node_1_coord_x_m - middle_node_element_x_l);
                                double a_jp1 = atan(temp1);
                                double a_j = atan(temp2);
                                //to check the angles
//                                cout << a_jp1*(180.0/M_PI) << endl;
//                                cout << a_j*(180.0/M_PI) << endl;
//                                cout << endl;
                                double temp3 = (a_jp1 - a_j)/(2*M_PI);
                                H_matrix(element_l,element_m) = temp3;

                                // pi = position of the middle node of element i, pj and pj+1 for the position of the nodes

                                //                                cout << "pi: " << "(" << middle_node_element_x_l << "," << middle_node_element_y_l << ")" << endl;
                                //                                cout << "pj: " << "(" << node_1_coord_x_m << "," << node_1_coord_y_m << ")" << endl;
                                //                                cout << "pj+1: " << "(" << node_2_coord_x_m << "," << node_2_coord_y_m << ")" << endl;
                                //                                cout << temp1 << endl;
                                //                                cout << temp2 << endl;
                                //                                cout << temp3 << endl;
                                //                                cout << endl;

                                //equivalent to lj in the slides
                                double length_element_m = sqrt((node_1_coord_x_m-node_2_coord_x_m)*(node_1_coord_x_m-node_2_coord_x_m)+(node_1_coord_y_m-node_2_coord_y_m)*(node_1_coord_y_m-node_2_coord_y_m));

                                string elementName;
                                int dim, order, numNodes, numPrimaryNodes;
                                vector<double> localNodeCoord;
                                model::mesh::getElementProperties(1,elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);

                                // get Gauss points coordinates and weights
                                vector<double> localCoords, weights;
                                model::mesh::getIntegrationPoints(1, integration_rule, localCoords, weights);

                                //To chech the weights and local coordinates given by gmsh

                                //                                cout << "\t" << weights.size() << " integration points\n";
                                //
                                //                                cout << "\tweights: ";
                                //                                for (int j = 0; j < weights.size(); ++j)
                                //                                    cout << weights[j] << " ";
                                //                                cout << '\n';
                                //
                                //                                cout << "\tlocalCoords (gmsh): ";
                                //                                for (int j = 0; j < localCoords.size(); ++j)
                                //                                    cout << localCoords[j] << " ";
                                //                                cout << '\n';

                                //compute the sum over the gauss points
                                double sum_gauss = 0;
                                vector<double> distance(weights.size(),0);
                                for(int g = 0; g < weights.size(); g++)
                                {
                                    //the variable xi is for the greek letter xi in the slides
                                    double xi = localCoords[3*g];
                                    //local coordinate x and y
                                    double x_local = (node_2_coord_x_m + node_1_coord_x_m)/2.0 + ((node_2_coord_x_m - node_1_coord_x_m)/2.0)*xi;
                                    double y_local = (node_2_coord_y_m + node_1_coord_y_m)/2.0 + ((node_2_coord_y_m - node_1_coord_y_m)/2.0)*xi;
                                    //compute the norm of the r vector
                                    double r = sqrt((x_local - middle_node_element_x_l)*(x_local - middle_node_element_x_l) + (y_local - middle_node_element_y_l)*(y_local - middle_node_element_y_l));
                                    distance[g] = r;
                                    //sum over the number of gauss points
                                    sum_gauss = sum_gauss + log(r)*weights[g];
                                }
                                distance;
                                G_matrix(element_l,element_m) = (length_element_m/(4.0*M_PI))*sum_gauss;
                            }
                            //create the matrix A and B(which need to be multiplied by the vector boundary conditions values to form vector b of the system Ax=b)
                            if(index[element_m] == 0)//0 = neumann
                            {
                                A_matrix(element_l,element_m) = H_matrix(element_l,element_m);
                                B_matrix(element_l,element_m) = G_matrix(element_l,element_m);
                            }
                            else if(index[element_m] == 1)//1 = dirichlet
                            {
                                A_matrix(element_l,element_m) = -G_matrix(element_l,element_m);
                                B_matrix(element_l,element_m) = -H_matrix(element_l,element_m);
                            }
                        }
                    }
                }
            }
        }
    }

    //Scans the index vector and displays whether it is a neumann or dirichlet condition
    cout << endl;
    cout << "index vector:" << endl;
    for(int i = 0; i < index.size(); i++)
    {
        if(index[i] == 0)
        {
            cout << "neumann" << endl;
        }
        if(index[i] == 1)
        {
            cout << "dirichlet" << endl;
        }
    }

    //display the differents matrices
    cout << endl;
    cout << "matrice G" << endl;
    show_mat(G_matrix);
    cout << endl;
    cout << "matrice H" << endl;
    show_mat(H_matrix);
    cout << endl;
    cout << "matrice A" << endl;
    show_mat(A_matrix);
    cout << endl;
    cout << "matrice B" << endl;
    show_mat(B_matrix);

    /*--------solve the linear system-----------------*/
    b_vector = B_matrix*boundary_conditions_values;

    cout << endl;
    cout << "vector boundary conditions values" << endl;
    show_vec(boundary_conditions_values);
    cout << endl;
    cout << "vector b" << endl;
    show_vec(b_vector);

    VectorXd x_vector = A_matrix.colPivHouseholderQr().solve(b_vector);

    cout << endl;
    cout << "vector x" << endl;
    show_vec(x_vector);

    //sort x_vector and boundary_conditions_values to gather the solution for neumann and for dirichlet
    sol_neumann.resize(N);
    sol_dirichlet.resize(N);

    for(int i = 0; i < N; i++)
    {
        if(index[i] == 0)
        {
            sol_neumann[i] = boundary_conditions_values[i];
            sol_dirichlet[i] = x_vector[i];
        }
        else if(index[i] == 1)
        {
            sol_dirichlet[i] = boundary_conditions_values[i];
            sol_neumann[i] = x_vector[i];
        }
    }

    cout << endl;
    cout << "solution neumann" << endl;
    show_vec(sol_neumann);

    cout << endl;
    cout << "solution dirichlet" << endl;
    show_vec(sol_dirichlet);

    //fltk::run(); // opens the gmsh window
    finalize();
    return 0;
}

//function to display the matrices
void show_mat(MatrixXd &A)
{
    for(int i = 0; i < A.rows(); i++)
    {
        cout << "\t";
        for(int j = 0; j < A.cols(); j++)
        {
            cout << "\t" << A(i,j);
        }
        cout << endl;
        cout << endl;
    }
}

//function to display the vectors
void show_vec(VectorXd &A)
{
    for(int i = 0; i < A.rows(); i++)
    {
        cout << "\t" << A(i);
        cout << endl;
        cout << endl;
    }
}

