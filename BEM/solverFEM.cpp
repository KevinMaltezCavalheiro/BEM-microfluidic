#include "functionsFEM.hpp"

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

int solverFEM(int argc, char **argv, map<int, double> &electrostaticPressure, const int nbViews)
{
    /*--------------INIT----------------*/
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    bool openGUI = false;        // Graphical interface of GMSH.

    //double start = omp_get_wtime();
    int max_threads = omp_get_max_threads();


    /*--------Integration rule used in the code----------*/
    // could be useful to pass it as an argument in the .geo file
    string integration_rule = "CompositeGauss4"; // tested with other methods, OK too.

    /*--------get physical groups--------------*/
    map<string, pair<int, int>> groups;
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (size_t i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    /*--------get nodes related to FEM-----------------*/
    string groupname = "FEM_domain"; // we consider the domain as containing the FEM 2d model
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;

    //loop over the entities of the domain to retrieve all nodeTags by taking care not to take some nodes twice
    //using sort and set_union methods
    vector<size_t> FEM_nodeTags;
    vector<double> FEM_nodecoord;
    vector<double> FEM_nodeparametricCoord;
    // Getting entities of the domain to fill FEM_nodeTags, en faire une fonction ?
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);
    vector<vector<size_t>> tmp_nodeTags(tags.size()); // will store the nodeTags associated to one given entity of the domain
    model::mesh::getNodes(tmp_nodeTags[0], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[0], true);
    sort(tmp_nodeTags[0].begin(), tmp_nodeTags[0].end()); // need to sort the vectors in order to merge them afterwards
    for (size_t i = 1; i < tags.size(); i++)
    { // in the case the domain contains multiple entities
        model::mesh::getNodes(tmp_nodeTags[i], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[i], true);
        sort(tmp_nodeTags[i].begin(), tmp_nodeTags[i].end());
        set_union(tmp_nodeTags[i-1].begin(), tmp_nodeTags[i-1].end(),
                       tmp_nodeTags[i].begin(), tmp_nodeTags[i].end(), back_inserter(FEM_nodeTags));
        
    }
    if(tags.size() == 1){ // if the domain contains one single entity
        FEM_nodeTags = tmp_nodeTags[0];
    }

    

    /*-------map between nodeTag and index of first nodal displacement associated to the node (second one is the same +1)------*/
    map<int, int> node_map;
    for (size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        int nodeTag = FEM_nodeTags[i];// use i à la place de nodeTag pour traiter cas ou le nodeTag n'est pas une suite continue
        int node_first_index = 2*i; 
        node_map[nodeTag] = node_first_index;
    }

    /*-----------get physical properties of domain----------------------*/
    vector<string> keys;
    onelab::getNames(keys, "(Volumic Forces|Materials).+");
    double E = 0; // Young modulus
    double nu = 0; // Poisson ratio
    double rho = 0; // volumic mass
    double bx = 0; // volumic force along x axis
    double by = 0; // volumic force along y axis

    Get_physical_Prop_Domain(keys, E, nu, rho, bx, by);

    /*-----------------H matrix involved in stress computation------------------*/
    Matrix<double, 3, 3> H_matrix = Compute_H_matrix(E, nu);
    //double startbis = omp_get_wtime();
    //cout << "preliminary: " << startbis-start << "\n";

    /*------------computing K matrix and f vector------------------------*/
    // we first consider the K matrix and the volumic part of the f vector
    SparseMatrix<double> full_K_matrix(2*FEM_nodeTags.size(),2*FEM_nodeTags.size());
    vector<double> full_f_vector(2*FEM_nodeTags.size());
    for (size_t i = 0; i < full_f_vector.size(); i++)
    {
        full_f_vector[i] = 0;
    }
    //vector<vector<double>> thread_f_vector(max_threads);

    vector<list<Triplet <double>>> my_k_list(max_threads);
    vector<list<pair<int,double>>> thread_f_vector(max_threads);

    // Looping over entities and compute K and F
    Loop_K_F(H_matrix, rho, bx, by, integration_rule, dim, tags, node_map, my_k_list, thread_f_vector);


    list<Triplet <double>> k_tripletList = my_k_list[0];
    //cout << "number of elements for thread 0: " << my_k_list[0].size() << "\n";

    if(max_threads > 1)
    {
        for (int i = 1; i < max_threads; i++)
        {
            //cout << "number of elements for thread " << i << ": " << my_k_list[i].size() << "\n";
            k_tripletList.insert(k_tripletList.end(), my_k_list[i].begin(), my_k_list[i].end());
        }
    }
    //cout << k_tripletList.size() << "\n";
    full_K_matrix.setFromTriplets(k_tripletList.begin(), k_tripletList.end());
    for(int i = 0; i < max_threads; i++){
        list<pair<int,double>> &tmp_vector = thread_f_vector[i];
        for(auto &p : tmp_vector) {
            full_f_vector[p.first] += p.second;
            //full_f_vector[tmp_vector[j].first] += tmp_vector[j].second;
        }
    }

    //double start2 = omp_get_wtime();
    //cout << "time for K assembly: " << start2-startbis << "\n";


    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    onelab::getNames(keys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    Neumann_BC(keys, integration_rule, node_map, groups, full_f_vector);

    //double start3 = omp_get_wtime();
    //cout << "time for neumann: " << start3 - start2 << "\n";

    /* ELECTROSTATIC PRESSURE */
    applyElecPressure(integration_rule, node_map, groups, electrostaticPressure, full_f_vector);
    //double start4 = omp_get_wtime();
    //cout << "time for elec pressure: " << start4 - start3 << "\n";

    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    // new_indices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not
    // if the given index has been removed, it will be set to -1 in new_indices
    vector<int> new_indices(2*FEM_nodeTags.size());

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    vector<double> dir_BC(2*FEM_nodeTags.size());

    // Fill the two vectors defined just above
    Fill_New_indices_Dir_BC(keys, node_map, groups, new_indices, dir_BC);

    int number_removed_lines = 0;
    for (size_t j = 0; j < new_indices.size(); ++j){
        if(new_indices[j] < 0){ // line has to be removed
            number_removed_lines++;
        }
        else{
            new_indices[j] -= number_removed_lines;
        }
    }

    //double start3 = omp_get_wtime();
    //cout << "time for boundary conditions: " << start3 - start2 << "\n";

    /*------------FILLING THE REDUCED K MATRIX AND F VECTOR--------------*/
    //the correction term on the f vector comes from the "FEM" theoretical course (JP Ponthot), Chap2 Slide 38/50
    SparseMatrix<double> reduced_K_matrix(2*FEM_nodeTags.size() - number_removed_lines,2*FEM_nodeTags.size() - number_removed_lines);
    vector<double> reduced_f_vector(2*FEM_nodeTags.size() - number_removed_lines);

    list<Triplet <double>> reduced_k_tripletList;

    Fill_K_F_Reduced(full_f_vector, full_K_matrix, new_indices, dir_BC, reduced_K_matrix, reduced_f_vector, reduced_k_tripletList);

    /*--------------SOLVING THE REDUCED SYSTEM--------------------*/
    VectorXd final_reduced_d(2*FEM_nodeTags.size() - number_removed_lines);

    VectorXd final_reduced_f(2*FEM_nodeTags.size() - number_removed_lines);
    for (size_t j = 0; j < reduced_f_vector.size(); ++j){
        final_reduced_f(j) = reduced_f_vector[j];
    }

    // solve Kd=f
    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(reduced_K_matrix);

    final_reduced_d = solver.solve(final_reduced_f);

    //double start4 = omp_get_wtime();
    //cout << "time to solve system: " << start4 -start3 << "\n";


    /*-----------RECONSTRUCTING FULL NODAL VALUES--------------*/
    vector<double> full_d_vector(2*FEM_nodeTags.size());
    for (size_t j = 0; j < full_d_vector.size(); ++j){
        if(new_indices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = final_reduced_d(new_indices[j]);
        }
    }

    /*--------------COMPUTING NODAL FORCES-------------*/
    VectorXd eigen_d(2*FEM_nodeTags.size());
    for(size_t i = 0 ; i < 2*FEM_nodeTags.size(); ++i){
        eigen_d(i) = full_d_vector[i];
    }
    VectorXd nodal_forces(2*FEM_nodeTags.size());
    nodal_forces = full_K_matrix*eigen_d;

    /*--------------REACTION FORCES COMPUTATION------------------*/

    double Rx = 0.;
    double Ry = 0.;
    Reaction_Forces(Rx, Ry, groups, node_map, nodal_forces, keys);

    //double start5 = omp_get_wtime();
    //cout << "time for rest: " << start5 - start4 << "\n";

    /*--------------PLOTTING AND SAVING NODAL VALUES------------------*/
    // loop over the nodes
    vector<string> names = Plot_Save_ux_uy(FEM_nodeTags, full_d_vector, nbViews);

    /*-------------POST PROCESSING---------------*/
    /*-------------Computing and saving strains/stresses-----------------*/

    // get all the elements in the domain

    groupname = "FEM_domain"; // we consider the domain as containing the FEM 2d model
    dim = groups[groupname].first;
    tag = groups[groupname].second;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elnodeTags;
    
    // compute total number of elements in the domain
    int nbelems=0;
    for(size_t i=0; i< tags.size(); ++i)
    {
        model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(size_t j=0; j< elementTags.size(); ++j)
        {
            nbelems+=elementTags[j].size();
        }
    }


    /*-------------------ELEMENTAL-NODAL VALUES-----------------------*/
    // éventuellement faire un tableau ici qui regroupe tout le bazar pour ne pas avoir la même ligne 7 fois ...
    int viewtag3 = view::add("nodal_eps_x");
    int viewtag4 = view::add("nodal_eps_y");
    int viewtag5 = view::add("nodal_2eps_xy");
    int viewtag6 = view::add("nodal_eps_z");
    int viewtag7 = view::add("nodal_sig_x");
    int viewtag8 = view::add("nodal_sig_y");
    int viewtag9 = view::add("nodal_sig_xy");
    int viewtag10 = view::add("nodal_sig_VM");
    vector<vector<double>> data3(nbelems);
    vector<vector<double>> data4(nbelems);
    vector<vector<double>> data5(nbelems);
    vector<vector<double>> data6(nbelems);
    vector<vector<double>> data7(nbelems);
    vector<vector<double>> data8(nbelems);
    vector<vector<double>> data9(nbelems);
    vector<vector<double>> data10(nbelems);
    vector<size_t> elems2D(nbelems);
    
    // computing the elemental nodal data for each physical group of the domain
    int k = 0; // current index of data vectors, k is incremented in Compute_Stress_Strain_Elemental_Nodal at each element
    for(size_t i=0; i< tags.size(); ++i)
    {
        model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        Compute_Stress_Strain_Elemental_Nodal(elementTypes, elementTags, elnodeTags, full_d_vector, elems2D, H_matrix, 
                                    nu, E, node_map, data3, data4, data5, data6, data7, data8, data9, data10, k);
    }

    view::addModelData(viewtag3, 0, names[0], "ElementNodeData",
                                elems2D, data3, 1);  // the last ,1 is important!
    view::addModelData(viewtag4, 0, names[0], "ElementNodeData",
                                elems2D, data4, 1);  // the last ,1 is important!
    view::addModelData(viewtag5, 0, names[0], "ElementNodeData",
                                elems2D, data5, 1);  // the last ,1 is important!
    view::addModelData(viewtag6, 0, names[0], "ElementNodeData",
                                elems2D, data6, 1);  // the last ,1 is important!
    view::addModelData(viewtag7, 0, names[0], "ElementNodeData",
                                elems2D, data7, 1);  // the last ,1 is important!
    view::addModelData(viewtag8, 0, names[0], "ElementNodeData",
                                elems2D, data8, 1);  // the last ,1 is important!
    view::addModelData(viewtag9, 0, names[0], "ElementNodeData",
                                elems2D, data9, 1);  // the last ,1 is important!
    view::addModelData(viewtag10, 0, names[0], "ElementNodeData",
                                elems2D, data10, 1);  // the last ,1 is important!

    // representing vector displacement field and nodal forces
    int viewtag11 = view::add("u");
    int viewtag12 = view::add("f");
    vector<vector<double>> data11(FEM_nodeTags.size());
    vector<vector<double>> data12(FEM_nodeTags.size());
    for (size_t i = 0; i < FEM_nodeTags.size(); ++i)
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

    view::addModelData(viewtag11, 0, names[0], "NodeData",
                             FEM_nodeTags, data11); //independent of time here
    view::addModelData(viewtag12, 0, names[0], "NodeData",
                             FEM_nodeTags, data12); //independent of time here



    // empeche que tout se plot l'un sur l'autre
    int nb_views = 12 + nbViews; //hardcoded
    for( int i = 0; i < nb_views; ++i){
        string my_string = "View[" + to_string(i) + "].Visible";
        option::setNumber(my_string, 0);
    }

    //double start6 = omp_get_wtime();
    //cout << "time for saving results: " << start6 - start5 << "\n";


    if(openGUI)
        fltk::run(); // opens the gmsh window
    
    //finalize();
    return 0;
}
