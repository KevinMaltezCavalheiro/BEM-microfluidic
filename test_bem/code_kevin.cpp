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

int code_py(int argc, char **argv);

void range_angle(double &angle);

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
    VectorXd sol_internal(0);
    int N_boundary = 0;//total number of elements
    int N_domain = 0;
    int N_tot = 0;


    //vectors to contain all the nodes information
    vector<size_t> all_nodeTags;
    vector<double> all_nodecoord;
    vector<double> all_nodeparametricCoord;
    vector<size_t> nodeTags_domain;
    vector<double> nodecoord_domain;
    vector<double> nodeparametricCoord_domain;
    vector<int> elementTypes_boundary;
    vector<vector<size_t> > elementTags_boundary;
    vector<vector<size_t> > nodeTags_boundary;
    vector<size_t> sort_nodes_1;
    vector<size_t> sort_nodes_2;

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

    model::mesh::getNodes(nodeTags_domain, nodecoord_domain, nodeparametricCoord_domain,2);

    model::mesh::getElements(elementTypes_boundary,elementTags_boundary,nodeTags_boundary,1,-1);

    //fill the map groups to gather all the information about the geometry
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    N_boundary = elementTags_boundary[0].size();//number of nodes = number of elements
    index.resize(N_boundary,-1);//set all the values at -1 to be sure it is fill with 0 or 1.
    boundary_conditions_values.resize(N_boundary);
    N_domain = nodeTags_domain.size();
    sol_internal.resize(N_domain);
    N_tot = N_domain+N_boundary;
    sort_nodes_1.resize(N_boundary);
    sort_nodes_2.resize(N_boundary);

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
    string integration_rule = "Gauss4";

    //defines all matrices related to the BEM method
    MatrixXd A_matrix;
    MatrixXd B_matrix;
    VectorXd b_vector;

    /*-----------------compute the matrix G and H------------------*/

    //resize all the matrices and initilize them to 0.
    A_matrix.resize(N_boundary,N_boundary);
    B_matrix.resize(N_boundary,N_boundary);
    b_vector.resize(N_boundary);
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
                sort_nodes_1[element_l] = node_1_tag_l;
                sort_nodes_2[element_l] = node_2_tag_l;
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
                            double g = 0;
                            double h = 0;
                            if(element_l == element_m)
                            {
                                //computation of Gii
                                double length_element_l = sqrt((node_1_coord_x_l-node_2_coord_x_l)*(node_1_coord_x_l-node_2_coord_x_l)+(node_1_coord_y_l-node_2_coord_y_l)*(node_1_coord_y_l-node_2_coord_y_l)); //equivalent to li in the slides
                                g = length_element_l/(2.0*M_PI)*(log(length_element_l/2.0)-1);

                                //computation of Hii
                                h = -1.0/2.0;
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
                                double num_ajp1 = (node_2_coord_y_m - middle_node_element_y_l);
                                double den_ajp1 = (node_2_coord_x_m - middle_node_element_x_l);
                                double num_aj = (node_1_coord_y_m - middle_node_element_y_l);
                                double den_aj = (node_1_coord_x_m - middle_node_element_x_l);
                                double a_jp1 = atan2(num_ajp1,den_ajp1);
                                double a_j = atan2(num_aj,den_aj);
                                //to check the angles
                                //                                cout << a_jp1*(180.0/M_PI) << endl;
                                //                                cout << a_j*(180.0/M_PI) << endl;
                                //                                cout << endl;

                                range_angle(a_j);
                                range_angle(a_jp1);

                                if(a_j > a_jp1)
                                {
                                    a_jp1 = a_jp1 + 2*M_PI;
                                }

                                h = (a_jp1 - a_j)/(2*M_PI);

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
                                for(int g = 0; g < weights.size(); g++)
                                {
                                    //the variable xi is for the greek letter xi in the slides
                                    double xi = localCoords[3*g];
                                    //local coordinate x and y
                                    double x_local = (node_2_coord_x_m + node_1_coord_x_m)/2.0 + ((node_2_coord_x_m - node_1_coord_x_m)/2.0)*xi;
                                    double y_local = (node_2_coord_y_m + node_1_coord_y_m)/2.0 + ((node_2_coord_y_m - node_1_coord_y_m)/2.0)*xi;
                                    //compute the norm of the r vector
                                    double r = sqrt((x_local - middle_node_element_x_l)*(x_local - middle_node_element_x_l) + (y_local - middle_node_element_y_l)*(y_local - middle_node_element_y_l));
                                    //sum over the number of gauss points
                                    sum_gauss = sum_gauss + log(r)*weights[g];
                                }
                                g = (length_element_m/(4.0*M_PI))*sum_gauss;
                            }
                            //create the matrix A and B(which need to be multiplied by the vector boundary conditions values to form vector b of the system Ax=b)
                            if(index[element_m] == 0)//0 = neumann
                            {
                                A_matrix(element_l,element_m) = h;
                                B_matrix(element_l,element_m) = g;
                            }
                            else if(index[element_m] == 1)//1 = dirichlet
                            {
                                A_matrix(element_l,element_m) = -g;
                                B_matrix(element_l,element_m) = -h;
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
    sol_neumann.resize(N_boundary);
    sol_dirichlet.resize(N_boundary);

    for(int i = 0; i < N_boundary; i++)
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

    for(int i = 0; i < N_domain; i++)
    {
        unsigned long node_i_tag = nodeTags_domain[i];
        double node_i_coord_x = all_nodecoord[3*(node_i_tag-1)];
        double node_i_coord_y = all_nodecoord[3*(node_i_tag-1)+1];
        double potential = 0;

        // second loop over the elements: indices j -> n -> m
        for(int j = 0; j < dimTags.size(); j++)
        {
            int dim = dimTags[j].first;
            int tag = dimTags[j].second;
            string name;
            model::getPhysicalName(dim, tag, name);
            groups[name] = {dim, tag};
            vector<int> tags;
            model::getEntitiesForPhysicalGroup(dim, tag, tags);
            for(int n = 0; n < tags.size(); n++)
            {
                vector<size_t> elementTags, nodeTags;
                model::mesh::getElementsByType(1,elementTags,nodeTags,tags[n]);
                for(int m = 0; m < elementTags.size(); m++)
                {
                    double g = 0;
                    double h = 0;
                    unsigned long element = elementTags[m]-1;
                    unsigned long node_1_tag = nodeTags[2*m];
                    unsigned long node_2_tag = nodeTags[2*m+1];
                    double node_1_coord_x = all_nodecoord[3*(node_1_tag-1)];
                    double node_1_coord_y = all_nodecoord[3*(node_1_tag-1)+1];
                    double node_2_coord_x = all_nodecoord[3*(node_2_tag-1)];
                    double node_2_coord_y = all_nodecoord[3*(node_2_tag-1)+1];
                    double length_element = sqrt((node_1_coord_x-node_2_coord_x)*(node_1_coord_x-node_2_coord_x)+(node_1_coord_y-node_2_coord_y)*(node_1_coord_y-node_2_coord_y));

                    double middle_node_element_x = (node_1_coord_x+node_2_coord_x)/2.0;
                    double middle_node_element_y = (node_1_coord_y+node_2_coord_y)/2.0;
                    double num_ajp1 = (node_2_coord_y - node_i_coord_y);
                    double den_ajp1 = (node_2_coord_x - node_i_coord_x);
                    double num_aj = (node_1_coord_y - node_i_coord_y);
                    double den_aj = (node_1_coord_x - node_i_coord_x);
                    double a_jp1 = atan2(num_ajp1,den_ajp1);
                    double a_j = atan2(num_aj,den_aj);

                    range_angle(a_j);
                    range_angle(a_jp1);

                    if(a_j > a_jp1)
                    {
                        a_jp1 = a_jp1 + 2*M_PI;
                    }

                    h = (a_jp1 - a_j)/(2*M_PI);

                    vector<double> localCoords, weights;
                    model::mesh::getIntegrationPoints(1, integration_rule, localCoords, weights);

                    double sum_gauss = 0;
                    for(int k = 0; k < weights.size(); k++)
                    {
                        double xi = localCoords[3*k];
                        double x_local = (node_2_coord_x + node_1_coord_x)/2.0 + ((node_2_coord_x - node_1_coord_x)/2.0)*xi;
                        double y_local = (node_2_coord_y + node_1_coord_y)/2.0 + ((node_2_coord_y - node_1_coord_y)/2.0)*xi;
                        double r = sqrt((x_local - node_i_coord_x)*(x_local - node_i_coord_x) + (y_local - node_i_coord_y)*(y_local - node_i_coord_y));
                        sum_gauss = sum_gauss + log(r)*weights[k];
                    }
                    g = (length_element/(4.0*M_PI))*sum_gauss;
                    potential = potential + h*sol_dirichlet[element] - g*sol_neumann[element];
                }
            }
        }
        sol_internal[i] = potential;
    }

    cout << endl;
    cout << "solution interne" << endl;
    show_vec(sol_internal);

    /*--------------PLOTTING AND SAVING THEM------------------*/

    vector<vector<double>> data_boundary(N_boundary, vector<double> (1, 0));
    vector<vector<double>> data_domain(N_domain, vector<double> (1, 0));

    for(int i = 0; i < N_boundary; i++)
    {
        data_boundary[i][0] = sol_dirichlet[i];
    }
    for(int i = 0; i < N_domain; i++)
    {
        data_domain[i][0] = sol_internal[i];
    }

    vector<string> names_boundary;
    vector<string> names_domain;
    model::list(names_boundary);
    model::list(names_domain);

    int viewNodeTag = view::add("NodeData");
    //view::addModelData(viewNodeTag, 0, names_boundary[0], "NodeData", /*elementTags_boundary[0]*/ sort_nodes_1, data_boundary, 0.0, 1);
    view::addModelData(viewNodeTag, 0, names_boundary[0], "NodeData", /*elementTags_boundary[0]*/ sort_nodes_2, data_boundary, 0.0, 1);
    view::addModelData(viewNodeTag, 0, names_domain[0], "NodeData", nodeTags_domain, data_domain, 0.0, 1);
    view::write(viewNodeTag, "results_boundary.msh");
    view::write(viewNodeTag, "results_domain.msh");

    cout << "the file contains " << names_boundary.size() << " model names\n";
    cout << "the file contains " << names_domain.size() << " model names\n";

    option::setNumber("View[0].IntervalsType", 3);
    option::setNumber("View[0].AdaptVisualizationGrid", 1);
    option::setNumber("View[0].MaxRecursionLevel", 3);
    option::setNumber("View[0].TargetError", -0.0001);

    fltk::run(); // opens the gmsh window
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

void range_angle(double &angle)
{
    angle = fmod(angle, (2.0 * M_PI));
    if(angle < 0.0)
    {
        angle = angle + 2.0 * M_PI;
    }
}

void show_vec_double(vector<double> &A)
{
    for(int i = 0; i < A.size(); i++)
    {
        cout << "\t" << A[i];
        cout << endl;
        cout << endl;
    }
}

int code_py(int argc, char **argv)
{
    // Each element has all its information gathered as in this structure.
    struct elementStruct{
        int tag;            // tag of the element
        std::string bcType; // "dirichlet" or "neumann" (at the end, should replace this by an int or boolean)
        double bcValue [2]; // values of [u, u_n]
        int nodeTags [2];   // tags of the [first, second] nodes of the element
    };

    if (argc < 2)
        {
            std::cout << "Usage: " << argv[0] << " <geo_file>\n";
            return 0;
        }

        // Parameters for displaying variables (for debugging).
        int dispNodeCoord = 0;  // Coordinates of the nodes.
        int dispHierarchy = 0;  // Hierarchy (physical groups -> entities -> elements -> nodes).
        int dispAngles = 0;     // Angles used for the computation of h. /!\ Displays ~n^2 angles.
        int dispMatrices = 0;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
        int dispValues = 0;     // Elements and their values of u and u_n.
        int openGUI = 1;        // Graphical interface of GMSH.

        gmsh::initialize();

        gmsh::open(argv[1]);
        gmsh::model::mesh::generate(2);

        // Gets the dim and the tags of the physical groups
        gmsh::vectorpair physicalGroupsDimTags;
        gmsh::model::getPhysicalGroups(physicalGroupsDimTags, 1);

        // Gets the boundary conditions from onelab and create a map: "name of the physical group" -> ("dirichlet"/"neumann", value)
        std::map<std::string, std::pair<std::string, double>> bcMap;

        std::vector<std::string> bcKeys;
        gmsh::onelab::getNames(bcKeys, "Boundary Conditions.+");
        for(auto &key : bcKeys)
        {
            std::vector<double> values;
            gmsh::onelab::getNumber(key, values);

            std::stringstream ss(key);
            std::vector<std::string> words;
            std::string word;
            while(std::getline(ss, word, '/'))
                words.push_back(word);

            if(words.size() == 3)
                bcMap[words[1]] = {words[2], values[0]};
        }

        // Initialises a vector of elements (that will contain all the elements of interrest)
        std::vector<elementStruct> elementVector;

        // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
        int dirichletCount = 0;

        // Gets the tags, the coordinates and the parametric coordinates of the nodes.
        std::vector<std::size_t> allNodeTags;
        std::vector<double> allCoord;
        std::vector<double> allParametricCoord;
        gmsh::model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);

        // Print the coordinates of all nodes if dispNodeCoord != 0
        if(dispNodeCoord)
        {
            for(int i = 0; i < allNodeTags.size(); i++)
                std::cout << "Node '" << allNodeTags[i] << "': (" << allCoord[3*i] << "; " << allCoord[3*i + 1] << "; " << allCoord[3*i + 2] << ")\n";
            std::cout << "\n";
        }

        for (int i = 0; i < physicalGroupsDimTags.size(); ++i)
        {
            int physicalGroupsDim = physicalGroupsDimTags[i].first;
            int physicalGroupsTag = physicalGroupsDimTags[i].second;
            std::string name;
            gmsh::model::getPhysicalName(physicalGroupsDim, physicalGroupsTag, name);

            // Retrieves the type and the value of the boundary condition of this physical group.
            // /!\ Should check that this physical groups was well defined in the list of "Boundary Condition" !
            std::string bcType = bcMap[name].first;
            double bcValue = bcMap[name].second;

            // Prints the physical groups name, dim, tag, bcType and bcValue (if dispHierarchy != 0).
            if(dispHierarchy)
                std::cout << "Physical group '" << name << "' (" << physicalGroupsDim << "D) with tag '" << physicalGroupsTag << "' has a '" << bcType << "' type and the value " << bcValue << "\n";

            std::vector<int> entitiesTags;
            gmsh::model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entitiesTags);

            for (int k = 0; k < entitiesTags.size(); k++)
            {
                // Print the tag of the entities of each physical group (if dispHierarchy != 0).
                if(dispHierarchy)
                    std::cout << "\t contains entity '" << entitiesTags[k] << "' ";

                std::vector<std::size_t> elementTags;
                std::vector<std::size_t> nodeTags;
                gmsh::model::mesh::getElementsByType(1, elementTags, nodeTags, entitiesTags[k]);

                // Print the number of element of the given entity (if dispHierarchy != 0).
                if(dispHierarchy)
                    std::cout << "which has " << elementTags.size() << " elements:\n";

                for(int el = 0; el < elementTags.size(); el++)
                {
                    // Stores the tag, the boundary type, the boundary value and the node tags of the element.
                    int elTag = elementTags[el];

                    elementStruct element;
                    element.tag = elTag;
                    element.bcType = bcType;
                    element.nodeTags[0] = nodeTags[2*el];
                    element.nodeTags[1] = nodeTags[2*el + 1];
                    // The elementVector should be "sorted" => push from front or from back depending on element type.
                    if(bcType == "dirichlet")
                    {
                        element.bcValue[0] = bcValue;
                        elementVector.insert(elementVector.begin(), element);
                        dirichletCount++;
                    }
                    else if(bcType == "neumann")
                    {
                        element.bcValue[1] = bcValue;
                        elementVector.push_back(element);
                    }

                    // Print the tag of a given element and the tags of the nodes that composes it as well as their coordinates (if dispHierarchy != 0).
                    if(dispHierarchy)
                    {
                        std::cout << "\t\t-Element '" << elTag << "' (" << element.bcType << "):\n";
                        std::cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*(element.nodeTags[0]-1)] << "; " << allCoord[3*(element.nodeTags[0]-1) + 1] << ")\n";
                        std::cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*(element.nodeTags[1]-1)] << "; " << allCoord[3*(element.nodeTags[1]-1) + 1] << ")\n";
                    }
                }
                if(dispHierarchy)
                    std::cout << "\n";
            }
        }

        // Declaration of matrices A, B and vector c (the matrix type is perhaps not the most efficient).
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
        Eigen::MatrixXd c = Eigen::MatrixXd::Zero(elementVector.size(), 1);

        //  Double loop over each element in order to compute the element h and g that will be given to the matrices A and B.
        for(int i = 0; i < elementVector.size(); i++)
        {
            // All these usefull i/j values could be stored in the elementStruct type, so that we only compute it once.
            // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element 'i' here).
            int i1Tag = elementVector[i].nodeTags[0];
            int i2Tag = elementVector[i].nodeTags[1];

            // Get the coordinates of the nodes from their tag.
            double i1X = allCoord[3*(i1Tag-1)];
            double i1Y = allCoord[3*(i1Tag-1) + 1];
            double i2X = allCoord[3*(i2Tag-1)];
            double i2Y = allCoord[3*(i2Tag-1) + 1];

            // Compute the coordinates of the mid point of element 'i' and its length.
            double iMidX = (i1X + i2X) / 2;
            double iMidY = (i1Y + i2Y) / 2;
            double iLength = sqrt(pow(i1X - i2X, 2.0) + pow(i1Y - i2Y, 2.0));

            for(int j = 0; j < elementVector.size(); j++)
            {
                // Will be computed and stored in A and B.
                double h;
                double g;

                if(i == j)  // Formula for diagonal elements.
                {
                    h = -0.5;
                    g = (iLength/(2*M_PI)) * (log(iLength/2)-1);
                }
                else    // General case
                {
                    // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element 'j' here).
                    int j1Tag = elementVector[j].nodeTags[0];
                    int j2Tag = elementVector[j].nodeTags[1];

                    // Get the coordinates of the nodes from their tag.
                    double j1X = allCoord[3*(j1Tag-1)];
                    double j1Y = allCoord[3*(j1Tag-1) + 1];
                    double j2X = allCoord[3*(j2Tag-1)];
                    double j2Y = allCoord[3*(j2Tag-1) + 1];

                    // Compute the coordinates of the mid point of element 'j', deltaX and deltaY and its length.
                    double jMidX = (j1X + j2X) / 2;
                    double jMidY = (j1Y + j2Y) / 2;
                    double jDeltaX = j2X - j1X;
                    double jDeltaY = j2Y - j1Y;
                    double jLength = sqrt(pow(jDeltaX, 2.0) + pow(jDeltaY, 2.0));

                    // Compute the angles between mid point of element 'i' and nodes 1/2 of element 'j' (in radian).
                    double a1 = atan2(j1Y - iMidY, j1X - iMidX);
                    double a2 = atan2(j2Y - iMidY, j2X - iMidX);
                    while(a2 < a1)  // This is not general and should be reviewed (concave surfaces)!
                        a2 += 2*M_PI;

                    // Angle verification
                    if(dispAngles)
                    {
                        double angleConst = 180 / M_PI;
                        double a1_deg = angleConst * a1;
                        double a2_deg = angleConst * a2;
                        std::cout << "Angles from element with tag '" << elementVector[i].tag << "' to element with tag '" << elementVector[j].tag << "':\n";
                        std::cout << "\tNode: tag '" << j1Tag << "': " << a1_deg << "\n";
                        std::cout << "\tNode: tag '" << j2Tag << "': " << a2_deg << "\n";
                    }

                    h = (a2 - a1) / (2*M_PI);

                    // Computation of the integration points (parametric) and weights.
                    std::string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.
                    std::vector<double> localCoord;
                    std::vector<double> weights;
                    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

                    // Computation of the sum (of integrand evaluated at integration point * weight)
                    double gSum = 0;
                    for(int k = 0; k < weights.size(); k++)
                    {
                        // Parametric coordinate
                        double xi = localCoord[3*k];
                        // Global coordinates
                        double gX = jMidX + jDeltaX*xi/2;
                        double gY = jMidY + jDeltaY*xi/2;

                        // Distance between mid point of element 'i' and integration point.
                        double r = sqrt(pow(iMidX - gX, 2.0) + pow(iMidY - gY, 2.0));
                        gSum += log(r) * weights[k];
                    }

                    g = gSum * jLength / (4*M_PI);
                }

                // Storage of h and g at the right place (in matrices A and B).
                if(j < dirichletCount)
                {
                    A(i,j) = g;
                    B(i,j) = h;
                    c(j) = elementVector[j].bcValue[0];
                }
                else
                {
                    A(i,j) = -h;
                    B(i,j) = -g;
                    c(j) = elementVector[j].bcValue[1];
                }
            }
        }

        // Compute b then x.
        Eigen::MatrixXd b = B*c;
        //Eigen::MatrixXd x = A.householderQr().solve(b);
        Eigen::MatrixXd x = A.fullPivLu().solve(b);

        // Matrices verification
        if(dispMatrices)
        {
            std::cout << "A matrix:\n";
            std::cout << A << std::endl;
            std::cout << "\nB matrix:\n";
            std::cout << B << std::endl;
            std::cout << "\nc vector:\n";
            std::cout << c << std::endl;
            std::cout << "\nb vector:\n";
            std::cout << b << std::endl;
            std::cout << "\nx vector:\n";
            std::cout << x << std::endl;
        }

        // Assign each value of x to the right element.
        for(int i = 0; i < elementVector.size(); i++)
        {
            if(i < dirichletCount)
                elementVector[i].bcValue[1] = x(i);
            else
                elementVector[i].bcValue[0] = x(i);

            if(dispValues)
                std::cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
        }

        /*
         *  Visualisation
         */

        std::vector<std::string> names;
        gmsh::model::list(names);
        std::cout << "the file contains " << names.size() << " model names\n";

        int viewNodeTag = gmsh::view::add("NodeData");
        std::vector<std::vector<double>> data(allNodeTags.size());

        for(int i = 0; i < allNodeTags.size(); i++)
        {
            // Coordinates of the node
            double x = allCoord[3*i];
            double y = allCoord[3*i + 1];

            double phi = 0;
            for(int j = 0; j < elementVector.size(); j++)
            {
                double h;
                double g;

                // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element i here).
                int j1Tag = elementVector[j].nodeTags[0];
                int j2Tag = elementVector[j].nodeTags[1];

                // This is certainly not the most efficient way to do.
                // Should also find an alternative for the ambiguity at the corner nodes.
                if(j1Tag == i+1 || j2Tag == i+1)
                {
                    phi = elementVector[j].bcValue[0];
                    break;
                }

                // Get the coordinates of the nodes from their tag.
                double j1X = allCoord[3*(j1Tag-1)];
                double j1Y = allCoord[3*(j1Tag-1) + 1];
                double j2X = allCoord[3*(j2Tag-1)];
                double j2Y = allCoord[3*(j2Tag-1) + 1];

                // Compute the coordinates of the mid point of element i and its length.
                double jMidX = (j1X + j2X) / 2;
                double jMidY = (j1Y + j2Y) / 2;
                double jDeltaX = j2X - j1X;
                double jDeltaY = j2Y - j1Y;
                double jLength = sqrt(pow(jDeltaX, 2.0) + pow(jDeltaY, 2.0));

                double a1 = atan2((j1Y - y),(j1X - x));
                double a2 = atan2((j2Y - y),(j2X - x));
                while(a2 < a1)
                    a2 += 2*M_PI;
                h = (a2 - a1) / (2*M_PI);

                // Computation of the integration points (parametric) and weights.
                std::string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.
                std::vector<double> localCoord;
                std::vector<double> weights;
                gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

                // Computation of the sum (of integrand evaluated at integration point * weight)
                double gSum = 0;
                for(int k = 0; k < weights.size(); k++)
                {
                    // Parametric coordinate
                    double xi = localCoord[3*k];
                    // Global coordinates
                    double gX = jMidX + jDeltaX*xi/2;
                    double gY = jMidY + jDeltaY*xi/2;

                    // Distance between mid point of element i and integration point.
                    double r = sqrt(pow(x - gX, 2.0) + pow(y - gY, 2.0));
                    gSum += log(r) * weights[k];
                }

                g = gSum * jLength / (4*M_PI);

                phi += h*elementVector[j].bcValue[0] - g*elementVector[j].bcValue[1];
            }
            // Potential at the node is phi.
            data[i].resize(1);
            data[i][0] = phi;
        }

    for(int i = 0; i < data.size(); i++)
    {
        cout << data[i][0] << endl;
    }

        gmsh::view::addModelData(viewNodeTag, 0, names[0], "NodeData", allNodeTags, data, 0.0, 1);
        gmsh::view::write(viewNodeTag, "results.msh");

        gmsh::option::setNumber("View[0].IntervalsType", 3);
        gmsh::option::setNumber("View[0].AdaptVisualizationGrid", 1);
        gmsh::option::setNumber("View[0].MaxRecursionLevel", 3);
        gmsh::option::setNumber("View[0].TargetError", -0.0001);

        if(openGUI)
            gmsh::fltk::run();

        gmsh::finalize();
        return 0;
}


