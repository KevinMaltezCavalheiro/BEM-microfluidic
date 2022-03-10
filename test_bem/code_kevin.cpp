#include <gmsh.h>       // gmsh main header
#include <iostream>     // for cout
#include <map>
#include <sstream>  // for stringstream
#include <Eigen/Dense> // Eigen library
#include <Eigen/Sparse> // Eigen library for sparse matrices (QUESTION: OK to use it ?)
#include <Eigen/SparseCholesky> // solving sparse linear systems (QUESTION: OK to use it ?)
#include <cmath>
#include "code_kevin.hpp"

using namespace gmsh;
using namespace std;
using namespace Eigen;

void show_mat(Matrix<double, Dynamic, Dynamic> &A);

int code_kevin(int argc, char **argv);

//todo:séparer en différentes fonctions
//todo:gérer les conditions de dirichlet et neumann
//todo: comparer au livre BEM
//todo: générer matrice H et G
// ordre (pour l'instant) initialisation, générer G, générer H

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
    int number_removed_lines = 0;
    vector<string> names;
    vectorpair dimTags;

    //physical properties of the domain
    vector<string> keys;
    string groupname = "boundary"; // on considère la frontière du domaine

    /*--------Integration rule used in the code----------*/
    // could be useful to pass it as an argument in the .geo file
    string integration_rule = "CompositeGauss4"; // tested with other methods, OK too.
    //H matrix
    Matrix<double, Dynamic, Dynamic> H_matrix; //projet d'avenir
    Matrix<double, Dynamic, Dynamic> G_matrix; //projet d'avenir

    /*-----------------compute the matrix G------------------*/

    //compute_G_matrix(groups, all_nodeTags, all_nodecoord, all_nodeparametricCoord, keys, groupname, tags, integration_rule, dimTags); //en travaux pour l'instant


    //initialise_G_matrix(all_nodeTags, groups, integration_rule,bx,rho,by,G_matrix,full_G_matrix,groupname);

    /*-----------------compute the matrix H------------------*/
    //initialise_H_matrix(H_matrix,nu,E);

    /*--------get physical groups--------------*/
    model::getPhysicalGroups(dimTags,1);// 1 car dans notre cas la frontière est de dimension 1.
    model::mesh::getNodes(all_nodeTags, all_nodecoord, all_nodeparametricCoord);
    G_matrix.resize(all_nodeTags.size(), all_nodeTags.size());
    H_matrix.resize(all_nodeTags.size(), all_nodeTags.size());
    G_matrix.array() = 0.;
    H_matrix.array() = 0.;
    show_mat(G_matrix);
    cout << endl;
    show_mat(H_matrix);

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
                unsigned long element_l = elementTags_1[l]-1;
                unsigned long node_1_tag_l = nodeTags_1[2*l];
                unsigned long node_2_tag_l = nodeTags_1[2*l+1];
                double node_1_coord_x_l = all_nodecoord[3*(node_1_tag_l-1)];
                double node_1_coord_y_l = all_nodecoord[3*(node_1_tag_l-1)+1];
                double node_2_coord_x_l = all_nodecoord[3*(node_2_tag_l-1)];
                double node_2_coord_y_l = all_nodecoord[3*(node_2_tag_l-1)+1];

                /*-------affiche les différents tags et les coordonnées de noeuds-------*/
                //                cout << "element " << elementTags[l] << " nodes = " << nodeTags[2*l] << "," << nodeTags[2*l+1] <<
                //                 " x1= " << all_nodecoord[3*(nodeTags[2*l]-1)] << " ,y1= " << all_nodecoord[3*(nodeTags[2*l]-1)+1] <<
                //                " x2= " << all_nodecoord[3*(nodeTags[2*l+1]-1)] << " ,y2= " << all_nodecoord[3*(nodeTags[2*l+1]-1)+1]<<               endl;
                //cout << "element " << element << " nodes = " << node_1_tag << "," << node_2_tag << " x1= " << node_1_coord_x << " ,y1= " << node_1_coord_y << " x2= " << node_2_coord_x << " ,y2= " << node_2_coord_y << endl;

                /*-------calcule ne noeud au milieu du segment-------*/
                double middle_node_element_x_l = (node_1_coord_x_l+node_2_coord_x_l)/2.0;
                double middle_node_element_y_l = (node_1_coord_y_l+node_2_coord_y_l)/2.0;
                //cout << "middle node element: (" << middle_node_element_x << "," << middle_node_element_y << ")" << endl;

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
                        for(int m = 0; m < elementTags_2.size(); m++)//chez nous element type = 1;
                        {
                            unsigned long element_m = elementTags_2[m]-1;
                            if(m == l)
                            {
                                double length_element_l = sqrt((node_1_coord_x_l-node_2_coord_x_l)*(node_1_coord_x_l-node_2_coord_x_l)+(node_1_coord_y_l-node_2_coord_y_l)*(node_1_coord_y_l-node_2_coord_y_l));
                                double g_ii = length_element_l/(2.0*M_PI)*(log(length_element_l/2.0)-1);
                                G_matrix(element_l,element_m) = g_ii;
                                double h_ii = -1.0/2.0;
                                H_matrix(element_l,element_m) = h_ii;
                            }
                            if(m != l)
                            {
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
                                double temp3 = (a_jp1 - a_j)/(2*M_PI);
                                H_matrix(element_l,element_m) = temp3;
//                                cout << "pi: " << "(" << middle_node_element_x_l << "," << middle_node_element_y_l << ")" << endl;
//                                cout << "pj: " << "(" << node_1_coord_x_m << "," << node_1_coord_y_m << ")" << endl;
//                                cout << "pj+1: " << "(" << node_2_coord_x_m << "," << node_2_coord_y_m << ")" << endl;
//                                cout << temp1 << endl;
//                                cout << temp2 << endl;
//                                cout << temp3 << endl;
//                                cout << endl;
                                double length_element_m = sqrt((node_1_coord_x_m-node_2_coord_x_m)*(node_1_coord_x_m-node_2_coord_x_m)+(node_1_coord_y_m-node_2_coord_y_m)*(node_1_coord_y_m-node_2_coord_y_m));

                                string elementName;
                                int dim, order, numNodes, numPrimaryNodes;
                                vector<double> localNodeCoord;
                                model::mesh::getElementProperties(1,elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);

                                // get Gauss points coordinates and weights
                                vector<double> localCoords, weights;
                                model::mesh::getIntegrationPoints(1, integration_rule, localCoords, weights);
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
                                double sum_gauss = 0;
                                for(int g = 0; g < weights.size(); g++)
                                {
                                    double xi_local_x = localCoords[3*g];
                                    double xi_local_y = localCoords[3*g+1];
                                    double x_local = (node_2_coord_x_m + node_1_coord_x_m)/2.0 + ((node_2_coord_x_m - node_1_coord_x_m)/2.0)*xi_local_x;
                                    double y_local = (node_2_coord_y_m + node_1_coord_y_m)/2.0 + ((node_2_coord_y_m - node_1_coord_y_m)/2.0)*xi_local_y;
                                    double r = sqrt((x_local - middle_node_element_x_l)*(x_local - middle_node_element_x_l) + (y_local - middle_node_element_y_l)*(y_local - middle_node_element_y_l));
                                    sum_gauss = sum_gauss + log(r)*weights[g];
                                }
                                G_matrix(element_l,element_m) = (length_element_m/(4.0*M_PI))*sum_gauss;
                            }
                        }
                    }
                }
            }
            //double boucle sur les éléments, faire G et H
            //première boucle
            //position du point i milieu
            //calculer aj et aj+1
            //calculer tous les élements G et H, vérifier s'ils sont bon
            //trier en fonction des conditions

            //attention il y a une somme sur les j pour H et G avec les u
            //parcourir la boucle 2 fois pour trier ce qui est neumann et ce qui est dirichlet
        }
        //cout << "Entities in group named: " << name << endl;
    }

    cout << endl;
    show_mat(G_matrix);
    cout << endl;
    show_mat(H_matrix);

    /*--------get number of nodes-----------------*/

    /*for(int i = 0; i < all_nodecoord.size(); i++)
     {
     cout << all_nodecoord[i] << endl;
     }

     for(int i = 0; i < all_nodeTags.size(); i++)
     {
     double x = all_nodecoord[i*3];
     double y = all_nodecoord[i*3+1];
     double z = all_nodecoord[i*3+2];
     cout << "x:" << x << " y:" << y << " z:" << z << endl;
     }*/


    //fltk::run(); // opens the gmsh window
    finalize();
    return 0;
}

void show_mat(Matrix<double, Dynamic, Dynamic> &A)
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

