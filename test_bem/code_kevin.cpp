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

int initialisation(map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys, string &groupname,vector<int> &tags,string &integration_rule, vectorpair &dimTags);

int compute_G_matrix(map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys, string &groupname,vector<int> &tags,string &integration_rule, vectorpair &dimTags);

void sorting_nodes(vector<size_t> &nodeTags,vector<double> &boundary_nodes_coordinates_x,vector<double> &boundary_nodes_coordinates_y,vector<double> &boundary_sort_nodes_tags, vector<double> &boundary_sort_length_elements);

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
    vector<int> tags;
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
    Matrix<double, 3, 3> H_matrix; //projet d'avenir
    Matrix<double, 3, 3> G_matrix; //projet d'avenir

    /*FONCTIONS*/
    //initialise toute la frontière physique
    initialisation(groups,all_nodeTags,all_nodecoord,all_nodeparametricCoord,keys,groupname,tags,integration_rule, dimTags);

    /*-----------------compute the matrix G------------------*/
    //Matrix<double,3,3> full_G_matrix;

    //compute_G_matrix(groups, all_nodeTags, all_nodecoord, all_nodeparametricCoord, keys, groupname, tags, integration_rule, dimTags); //en travaux pour l'instant


    //initialise_G_matrix(all_nodeTags, groups, integration_rule,bx,rho,by,G_matrix,full_G_matrix,groupname);

    /*-----------------compute the matrix H------------------*/
    //initialise_H_matrix(H_matrix,nu,E);

    fltk::run(); // opens the gmsh window
    finalize();
    return 0;
}

int initialisation(map<string, pair<int, int>> &groups,vector<size_t> &all_nodeTags,vector<double> &all_nodecoord,vector<double> &all_nodeparametricCoord,vector<string> &keys, string &groupname,vector<int> &tags,string &integration_rule,vectorpair &dimTags)
{
    /*--------get physical groups--------------*/
    model::getPhysicalGroups(dimTags,1);// 1 car dans notre cas la frontière est de dimension 1.
    model::mesh::getNodes(all_nodeTags, all_nodecoord, all_nodeparametricCoord);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        //si mon name c'est dirichlet, on peut trier les conditions de bords et les mettre dans dirichlet par après on pourra parcourir la matrice et trier G et H en fonctions de ces dernières
        groups[name] = {dim, tag};
        model::getEntitiesForPhysicalGroup(dim, tag, tags);

        //check qu'il n'y a pas de trous dans la numérotation des noeuds
        //on suppose que all_nodeTags n'a pas de trou, est continu et commence à 1
        //check via une petite routine
        //model::mesh::renumberNodes();
        for(int k = 0; k<tags.size(); k++)
        {
            vector<size_t> elementTags, nodeTags;
            model::mesh::getElementsByType(1,elementTags,nodeTags,tags[k]);// on a que des segments droits
            //cout << elementTags.size() << " " << nodeTags.size() << endl;
            for(int l = 0; l < elementTags.size(); l++)
            {
                int element = elementTags[l];
                int node_1_tag = nodeTags[2*l];
                int node_2_tag = nodeTags[2*l+1];
                double node_1_coord_x = all_nodecoord[3*(node_1_tag-1)];
                double node_1_coord_y = all_nodecoord[3*(node_1_tag-1)+1];
                double node_2_coord_x = all_nodecoord[3*(node_2_tag-1)];
                double node_2_coord_y = all_nodecoord[3*(node_2_tag-1)+1];
//                cout << "element " << elementTags[l] << " nodes = " << nodeTags[2*l] << "," << nodeTags[2*l+1] <<
//                 " x1= " << all_nodecoord[3*(nodeTags[2*l]-1)] << " ,y1= " << all_nodecoord[3*(nodeTags[2*l]-1)+1] <<
//                " x2= " << all_nodecoord[3*(nodeTags[2*l+1]-1)] << " ,y2= " << all_nodecoord[3*(nodeTags[2*l+1]-1)+1]<<               endl;
                cout << "element " << element << " nodes = " << node_1_tag << "," << node_2_tag << " x1= " << node_1_coord_x << " ,y1= " << node_1_coord_y << " x2= " << node_2_coord_x << " ,y2= " << node_2_coord_y << endl;
            }
            //double boucle sur les éléments, faire G et H
            //première boucle
            //position du point i milieu
            //calculer aj et aj+1
            //calculer tous les élements G et H, vérifier s'ils sont bon
            //trier en fonction des conditions
        }
        //cout << "Entities in group named: " << name << endl;
    }

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

    return 0;
}
