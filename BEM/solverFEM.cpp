#include "functionsFEM.hpp"

#ifdef _MSC_VER
#include <gmsh.h_cwrap>       // gmsh main header
#else
#include <gmsh.h>       // gmsh main header
#endif

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

void solverFEMnonLinear(std::map<int, double> &electrostaticPressure, std::map<int,std::pair<double,double>> & boundaryDisplacementMap, const int nbViews, bool postProcessing, const int iteration, const int viewTagU, const int viewTagF)
{
    /*--------------INITIALIZATION OF ITERATIVE ALGORITHM----------------*/
    /* FIRST STEP: 
    - construct and compute the vector of applied forces
    - set up the kinematics of each element and compute the local stiffness matrix
    - initialize the full d vector with Dirichlet boundary conditions and set pg = 0
    - important remark: in pg we focus only on the free DOFs (we do not focus on the DOFs fixed by dirichlet BCs)*/

    bool openGUI = false;        // Graphical interface of GMSH.

    //double start = omp_get_wtime();
    int max_threads = omp_get_max_threads();


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

    /*--------get nodes related to FEM-----------------*/
    std::string groupname = "FEM_domain"; // we consider the domain as containing the FEM 2d model
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;

    //loop over the entities of the domain to retrieve all nodeTags by taking care not to take some nodes twice
    //using sort and set_union methods + retrieving the node coordinates and storing them in a map
    std::vector<std::size_t> FEM_nodeTags;
    std::vector<double> FEM_nodecoord;
    std::vector<double> FEM_nodeparametricCoord;
    std::map<int,std::pair<double,double>> FEM_nodeCoordsMap; // nodeTag -> (x,y)
    // Getting entities of the domain to fill FEM_nodeTags, en faire une fonction ?
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
    std::vector<std::vector<std::size_t>> tmp_nodeTags(tags.size()); // will store the nodeTags associated to one given entity of the domain
    gmsh::model::mesh::getNodes(tmp_nodeTags[0], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[0], true); 
    for(std::size_t j = 0; j < tmp_nodeTags[0].size(); j++)
    {
        FEM_nodeCoordsMap[tmp_nodeTags[0][j]] = std::pair<double, double>(FEM_nodecoord[3*j], FEM_nodecoord[3*j+1]);
    }
    std::sort(tmp_nodeTags[0].begin(), tmp_nodeTags[0].end()); // need to sort the vectors in order to merge them afterwards
    for (std::size_t i = 1; i < tags.size(); i++)
    { // in the case the domain contains multiple entities
        gmsh::model::mesh::getNodes(tmp_nodeTags[i], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[i], true);
        for(std::size_t j = 0; j < tmp_nodeTags[i].size(); j++)
        {
            FEM_nodeCoordsMap[tmp_nodeTags[i][j]] = std::pair<double, double>(FEM_nodecoord[3*j], FEM_nodecoord[3*j+1]);
        }
        std::sort(tmp_nodeTags[i].begin(), tmp_nodeTags[i].end());
        std::set_union(tmp_nodeTags[i-1].begin(), tmp_nodeTags[i-1].end(),
                       tmp_nodeTags[i].begin(), tmp_nodeTags[i].end(), std::back_inserter(FEM_nodeTags));
        
    }
    if(tags.size() == 1){ // if the domain contains one single entity
        FEM_nodeTags = tmp_nodeTags[0];
    }

    /*-----get elements related to FEM--------------*/
    std::vector<std::size_t> FEM_elementTags;
    std::map<int,elementData> FEM_kinematicsMap; // map between elementTag and all the kinematics associated to this element
    // map is initialized, see functionsFEM.hpp for structure of elementData type
    int nbelems = initKinematics(FEM_elementTags, FEM_kinematicsMap, FEM_nodeCoordsMap, dim, tags);

    /*-------map between nodeTag and index of first nodal displacement associated to the node (second one is the same +1)------*/
    std::map<int, int> nodeIndexMap;
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        int nodeTag = FEM_nodeTags[i];// use i à la place de nodeTag pour traiter cas ou le nodeTag n'est pas une suite continue
        int node_first_index = 2*i; 
        nodeIndexMap[nodeTag] = node_first_index;
    }

    /*-----------get physical properties of domain----------------------*/
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Volumic Forces|Materials).+");
    double E = 0; // Young modulus
    double nu = 0; // Poisson ratio
    double rho = 0; // volumic mass
    double bx = 0; // volumic force along x axis
    double by = 0; // volumic force along y axis

    getPhysicalProperties(keys, E, nu, rho, bx, by);

    /*-----------------H matrix involved in stress computation------------------*/
    Eigen::Matrix<double, 3, 3> H_matrix = computeHmatrix(E, nu);
    //double startbis = omp_get_wtime();
    //std::cout << "preliminary: " << startbis-start << "\n";

    /*------------computing full f vector------------------------*/
    // we first consider volumic part of the f vector, computed simultaneously with local K matrices
    std::vector<double> full_f_vector(2*FEM_nodeTags.size());
    for (size_t i = 0; i < full_f_vector.size(); i++)
    {
        full_f_vector[i] = 0;
    }

    std::vector<std::list<Eigen::Triplet <double>>> my_k_list(max_threads);
    std::vector<std::list<std::pair<int,double>>> thread_f_vector(max_threads);

    // Looping over entities and compute local K matrices and volumic f vector
    fillElementalKandF(H_matrix, rho, bx, by, integration_rule, dim, tags, nodeIndexMap, my_k_list, thread_f_vector, FEM_kinematicsMap);

    for(int i = 0; i < max_threads; i++){
        std::list<std::pair<int,double>> &tmp_vector = thread_f_vector[i];
        for(auto &p : tmp_vector) {
            full_f_vector[p.first] += p.second;
        }
    }

    //double start2 = omp_get_wtime();
    //std::cout << "time for K assembly: " << start2-startbis << "\n";


    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    gmsh::onelab::getNames(keys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    includeNeumannBC(keys, integration_rule, nodeIndexMap, groups, full_f_vector);

    //double start3 = omp_get_wtime();
    //std::cout << "time for neumann: " << start3 - start2 << "\n";

    /* taking ELECTROSTATIC PRESSURE */
    if(electrostaticPressure.size() != 0)
        applyElecPressure(integration_rule, FEM_nodeCoordsMap, nodeIndexMap, groups, electrostaticPressure, full_f_vector);
    //double start4 = omp_get_wtime();
    //std::cout << "time for elec pressure: " << start4 - start3 << "\n";

    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    // new_DOFindices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not
    // if the given index has been removed, it will be set to -1 in new_DOFindices
    std::vector<int> new_DOFindices(2*FEM_nodeTags.size());

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    std::vector<double> dir_BC(2*FEM_nodeTags.size());

    // Fill the two vectors defined just above
    fillDOFindicesDirBC(keys, nodeIndexMap, groups, new_DOFindices, dir_BC);

    int number_removed_lines = 0;
    for (std::size_t j = 0; j < new_DOFindices.size(); ++j){
        if(new_DOFindices[j] < 0){ // line has to be removed
            number_removed_lines++;
        }
        else{
            new_DOFindices[j] -= number_removed_lines;
        }
    }

    std::vector<double> full_d_vector(2*FEM_nodeTags.size());
    for (std::size_t j = 0; j < full_d_vector.size(); ++j){
        if(new_DOFindices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = 0;
        }
    }
    // initializing the applied forces vector (will remain constant) and the global p (p_g) vector
    std::vector<double> f_app_vector(2*FEM_nodeTags.size() - number_removed_lines);
    std::vector<double> global_p_vector(2*FEM_nodeTags.size() - number_removed_lines);
    int nb_free_DOFs = 0;
    for(size_t i = 0; i < 2*FEM_nodeTags.size(); i++)
    {
        if(new_DOFindices[i] >= 0)
        {
            f_app_vector[nb_free_DOFs] = full_f_vector[i];
            global_p_vector[nb_free_DOFs] = full_d_vector[i];
            nb_free_DOFs++;
        }
    }

    // computing the norm of the applied forces vector
    double f_app_norm = 0;
    for(size_t i = 0; i < f_app_vector.size(); i++)
    {
        f_app_norm += f_app_vector[i]*f_app_vector[i];
    }
    f_app_norm = sqrt(f_app_norm);
    std::cout << "f app norm: " << f_app_norm << "\n";

    std::vector<double> global_f_vector(nb_free_DOFs);

    Eigen::Matrix<double, Eigen::Dynamic, 1> residual(nb_free_DOFs, 1);

    double residualNorm = 100000000; // hardcoded and arbitrary
    double residualTolerance = 2e-12; // empirical

    /*double max_applied_force = 0;
    for(int i = 0; i < nb_free_DOFs; i++)
    {
        if(abs(f_app_vector[i]) > max_applied_force)
            max_applied_force = abs(f_app_vector[i]);
    }
    std::cout << "max applied force: " << max_applied_force << "\n";*/


    /*---------ITERATIVE PART---------------*/
    int max_number_of_steps = 10; // purely arbitrary
    for(int step = 0; step < max_number_of_steps; step++)
    {
        /*------------STEP 2 of the algorithm----------------*/
        /* 
        - update the kinematics from last iteration (or from dirichlet boundary conditions)
        - compute the global nodal forces f_g based on previous results
        - compute the residual = norm(f_g - f_app)*/

        /*--------- UPDATING KINEMATICS---------------*/
        updateKinematics(FEM_kinematicsMap, full_d_vector, nodeIndexMap, FEM_elementTags);

        /*--------COMPUTING AND ASSEMBLING THE GLOBAL F VECTOR BASED ON PREVIOUS ITERATION-----------*/
        for(int i = 0; i < nb_free_DOFs; i++)
        {
            global_f_vector[i] = 0;
        }
        assembleGlobalF(global_f_vector, FEM_kinematicsMap, nodeIndexMap, new_DOFindices, FEM_elementTags);

        residualNorm = computeResidualNorm(global_f_vector, f_app_vector, residual);

        std::cout << "residual value: " << residualNorm/f_app_norm << "\n";

        /*------------STEP 3: EVALUATING THE RESIDUAL, have we converged ??----------*/
        if (residualNorm/f_app_norm < residualTolerance)
            break;

        /*------------STEP 4: ASSEMBLING THE GLOBAL TANGENT STIFFNESS MATRIX----------*/
        // it is assembled as a sparse matrix
        Eigen::SparseMatrix<double> global_tangent_matrix(nb_free_DOFs, nb_free_DOFs);
        assembleGlobalKg(global_tangent_matrix, FEM_kinematicsMap, nodeIndexMap, new_DOFindices, FEM_elementTags);

        //std::cout << "global tangent matrix: " << global_tangent_matrix << "\n";

        /*-----------STEP 5: SOLVING THE SYSTEM----------*/
        Eigen::Matrix<double, Eigen::Dynamic, 1> global_delta_p(nb_free_DOFs, 1);
        
        // solve Kg*Deltapg=-res
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(global_tangent_matrix);

        global_delta_p = solver.solve(-residual);

        /*----------STEP 6: updating the global p vector and the full d vector of nodal values-----------*/
        for(size_t i = 0; i < global_p_vector.size(); i++)
        {
            global_p_vector[i] += global_delta_p(i,0);
        }
        updateNodalDispVector(full_d_vector, global_p_vector, new_DOFindices);
    }
    
    //updating the kinematics according to the results of last iteration
    updateKinematics(FEM_kinematicsMap, full_d_vector, nodeIndexMap, FEM_elementTags);

    /*--------------COMPUTING NODAL FORCES-------------*/
    std::vector<double> final_nodal_forces(2*FEM_nodeTags.size());
    retrieveFinalNodalForces(final_nodal_forces, FEM_kinematicsMap, nodeIndexMap, FEM_elementTags);

    /*--------------RETRIEVING NODAL DISPLACEMENTS ON BEM-FEM BOUNDARY------------*/
    //those displacements will be sent to the BEM code in the non-linear iterative solver
    //std::map<int,std::pair<double,double>> boundaryDisplacementMap;
    // mapping nodeTag (of node on the boundary) into (u,v) displacements of the node
    retrieveBoundaryDisplacements(boundaryDisplacementMap, full_d_vector, nodeIndexMap, groups);

    std::vector<std::string> names;
    gmsh::model::list(names);

    if(postProcessing)
    {
        /*--------------REACTION FORCES COMPUTATION------------------*/
        computeReactionForces(groups, nodeIndexMap, final_nodal_forces, keys);

        //double start5 = omp_get_wtime();
        //std::cout << "time for rest: " << start5 - start4 << "\n";

        /*--------------PLOTTING NODAL VALUES------------------*/
        plotUxUy(FEM_nodeTags, full_d_vector, nbViews, names);

        /*-------------POST PROCESSING---------------*/
        /*-------------Computing and saving strains/stresses-----------------*/


        /*-------------------ELEMENTAL-NODAL VALUES-----------------------*/
        // éventuellement faire un tableau ici qui regroupe tout le bazar pour ne pas avoir la même ligne 7 fois ...
        int viewtag3 = gmsh::view::add("nodal_eps_x");
        int viewtag4 = gmsh::view::add("nodal_eps_y");
        int viewtag5 = gmsh::view::add("nodal_eps_xy");
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

        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> elnodeTags;
        
        // computing the elemental nodal data for each physical group of the domain
        int k = 0; // current index of data vectors, k is incremented in computeStressStrainElNodal at each element
        for(std::size_t i=0; i< tags.size(); ++i)
        {
            gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
            computeStressStrainElNodal(elementTypes, elementTags, elnodeTags, elems2D, H_matrix, 
                                        nu, E, FEM_kinematicsMap, data3, data4, data5, data6, data7, data8, data9, data10, k);
        }

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
    }
    // representing vector displacement field and nodal forces
    //int viewTagU = gmsh::view::add("u");
    //int viewTagF = gmsh::view::add("f");
    std::vector<std::vector<double>> data11(FEM_nodeTags.size());
    std::vector<std::vector<double>> data12(FEM_nodeTags.size());
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        data11[i].resize(3);
        data12[i].resize(3);
        data11[i][0] = full_d_vector[2*i]; //ux, here use NodeIndexMap
        data11[i][1] = full_d_vector[2*i+1]; //uy
        data11[i][2] = 0.; //uz
        data12[i][0] = final_nodal_forces[2*i]; //fx
        data12[i][1] = final_nodal_forces[2*i+1]; //fy
        data12[i][2] = 0.; //fz
    }

    gmsh::view::addModelData(viewTagU, iteration, names[0], "NodeData",
                             FEM_nodeTags, data11, iteration); //independent of time here
    gmsh::view::addModelData(viewTagF, iteration, names[0], "NodeData",
                             FEM_nodeTags, data12, iteration); //independent of time here


    if(postProcessing)
    {
        // empeche que tout se plot l'un sur l'autre
        int nb_views = 12 + nbViews; //hardcoded
        for( int i = 0; i < nb_views; ++i){
            std::string my_string = "View[" + std::to_string(i) + "].Visible";
            gmsh::option::setNumber(my_string, 0);
        }
        gmsh::option::setNumber("View[0].VectorType", 5); // plot automatiquement sous forme de déplacement
        gmsh::option::setNumber("View[0].RangeType", 3); // abscisses par step
    }


    if(openGUI && postProcessing)
        gmsh::fltk::run(); // opens the gmsh window
}


void solverFEM(std::map<int, double> &electrostaticPressure, const int nbViews)
{
    /*--------------INIT----------------*/

    bool openGUI = false;        // Graphical interface of GMSH.

    //double start = omp_get_wtime();
    int max_threads = omp_get_max_threads();


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

    /*--------get nodes related to FEM-----------------*/
    std::string groupname = "FEM_domain"; // we consider the domain as containing the FEM 2d model
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;

    //loop over the entities of the domain to retrieve all nodeTags by taking care not to take some nodes twice
    //using sort and set_union methods
    std::vector<std::size_t> FEM_nodeTags;
    std::vector<double> FEM_nodecoord;
    std::vector<double> FEM_nodeparametricCoord;
    // Getting entities of the domain to fill FEM_nodeTags, en faire une fonction ?
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
    std::vector<std::vector<std::size_t>> tmp_nodeTags(tags.size()); // will store the nodeTags associated to one given entity of the domain
    gmsh::model::mesh::getNodes(tmp_nodeTags[0], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[0], true); 
    std::sort(tmp_nodeTags[0].begin(), tmp_nodeTags[0].end()); // need to sort the vectors in order to merge them afterwards
    for (std::size_t i = 1; i < tags.size(); i++)
    { // in the case the domain contains multiple entities
        gmsh::model::mesh::getNodes(tmp_nodeTags[i], FEM_nodecoord, FEM_nodeparametricCoord, dim, tags[i], true);
        std::sort(tmp_nodeTags[i].begin(), tmp_nodeTags[i].end());
        std::set_union(tmp_nodeTags[i-1].begin(), tmp_nodeTags[i-1].end(),
                       tmp_nodeTags[i].begin(), tmp_nodeTags[i].end(), std::back_inserter(FEM_nodeTags));
        
    }
    if(tags.size() == 1){ // if the domain contains one single entity
        FEM_nodeTags = tmp_nodeTags[0];
    }

    

    /*-------map between nodeTag and index of first nodal displacement associated to the node (second one is the same +1)------*/
    std::map<int, int> nodeIndexMap;
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        int nodeTag = FEM_nodeTags[i];// use i à la place de nodeTag pour traiter cas ou le nodeTag n'est pas une suite continue
        int node_first_index = 2*i; 
        nodeIndexMap[nodeTag] = node_first_index;
    }

    /*-----------get physical properties of domain----------------------*/
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Volumic Forces|Materials).+");
    double E = 0; // Young modulus
    double nu = 0; // Poisson ratio
    double rho = 0; // volumic mass
    double bx = 0; // volumic force along x axis
    double by = 0; // volumic force along y axis

    getPhysicalProperties(keys, E, nu, rho, bx, by);

    /*-----------------H matrix involved in stress computation------------------*/
    Eigen::Matrix<double, 3, 3> H_matrix = computeHmatrix(E, nu);
    //double startbis = omp_get_wtime();
    //std::cout << "preliminary: " << startbis-start << "\n";

    /*------------computing K matrix and f vector------------------------*/
    // we first consider the K matrix and the volumic part of the f vector
    Eigen::SparseMatrix<double> full_K_matrix(2*FEM_nodeTags.size(),2*FEM_nodeTags.size());
    std::vector<double> full_f_vector(2*FEM_nodeTags.size());
    for (size_t i = 0; i < full_f_vector.size(); i++)
    {
        full_f_vector[i] = 0;
    }

    std::vector<std::list<Eigen::Triplet <double>>> my_k_list(max_threads);
    std::vector<std::list<std::pair<int,double>>> thread_f_vector(max_threads);

    // Looping over entities and compute full K (under tripletlist form) and F
    assembleKlistVolumicF(H_matrix, rho, bx, by, integration_rule, dim, tags, nodeIndexMap, my_k_list, thread_f_vector);


    std::list<Eigen::Triplet <double>> k_tripletList = my_k_list[0];

    if(max_threads > 1)
    {
        for (int i = 1; i < max_threads; i++)
        {
            k_tripletList.insert(k_tripletList.end(), my_k_list[i].begin(), my_k_list[i].end());
        }
    }
    full_K_matrix.setFromTriplets(k_tripletList.begin(), k_tripletList.end());
    for(int i = 0; i < max_threads; i++){
        std::list<std::pair<int,double>> &tmp_vector = thread_f_vector[i];
        for(auto &p : tmp_vector) {
            full_f_vector[p.first] += p.second;
        }
    }

    //double start2 = omp_get_wtime();
    //std::cout << "time for K assembly: " << start2-startbis << "\n";


    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    gmsh::onelab::getNames(keys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    includeNeumannBC(keys, integration_rule, nodeIndexMap, groups, full_f_vector);

    //double start3 = omp_get_wtime();
    //std::cout << "time for neumann: " << start3 - start2 << "\n";

    /* ELECTROSTATIC PRESSURE */
    if(electrostaticPressure.size() != 0)
    {
        std::string groupnameBoundary = "BEM_FEM_boundary";
        int dimBoundary = groups[groupnameBoundary].first;
        int tagBoundary = groups[groupnameBoundary].second;
        std::vector<int> tagsBoundary;
        gmsh::model::getEntitiesForPhysicalGroup(dimBoundary, tagBoundary, tagsBoundary);
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> elnodeTags;

        std::map<int,std::pair<double,double>> nodeCoordsMap;
        //filling the map with : nodeTag -> (x,y)
        for(std::size_t i=0; i< tags.size(); ++i)
        {
            std::vector<std::size_t> entity_nodeTags;
            std::vector<double> entity_nodecoord;
            std::vector<double> entity_nodeparametricCoord;
            gmsh::model::mesh::getNodes(entity_nodeTags, entity_nodecoord, entity_nodeparametricCoord, dim, tags[i], true); 
            for ( std::size_t j=0; j < entity_nodeTags.size(); j++){
                nodeCoordsMap[entity_nodeTags[j]] = std::pair<double, double>(entity_nodecoord[3*j], entity_nodecoord[3*j+1]);
            }
        }
        applyElecPressure(integration_rule, nodeCoordsMap, nodeIndexMap, groups, electrostaticPressure, full_f_vector);
    }
        
    //double start4 = omp_get_wtime();
    //std::cout << "time for elec pressure: " << start4 - start3 << "\n";

    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    // new_DOFindices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not
    // if the given index has been removed, it will be set to -1 in new_DOFindices
    std::vector<int> new_DOFindices(2*FEM_nodeTags.size());

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    std::vector<double> dir_BC(2*FEM_nodeTags.size());

    // Fill the two vectors defined just above
    fillDOFindicesDirBC(keys, nodeIndexMap, groups, new_DOFindices, dir_BC);

    int number_removed_lines = 0;
    for (std::size_t j = 0; j < new_DOFindices.size(); ++j){
        if(new_DOFindices[j] < 0){ // line has to be removed
            number_removed_lines++;
        }
        else{
            new_DOFindices[j] -= number_removed_lines;
        }
    }

    //double start3 = omp_get_wtime();
    //std::cout << "time for boundary conditions: " << start3 - start2 << "\n";

    /*------------FILLING THE REDUCED K MATRIX AND F VECTOR--------------*/
    //the correction term on the f vector comes from the "FEM" theoretical course (JP Ponthot), Chap2 Slide 38/50
    Eigen::SparseMatrix<double> reduced_K_matrix(2*FEM_nodeTags.size() - number_removed_lines,2*FEM_nodeTags.size() - number_removed_lines);
    std::vector<double> reduced_f_vector(2*FEM_nodeTags.size() - number_removed_lines);

    std::list<Eigen::Triplet <double>> reduced_k_tripletList;

    fillReducedKf(full_f_vector, full_K_matrix, new_DOFindices, dir_BC, reduced_K_matrix, reduced_f_vector, reduced_k_tripletList);

    /*--------------SOLVING THE REDUCED SYSTEM--------------------*/
    Eigen::VectorXd final_reduced_d(2*FEM_nodeTags.size() - number_removed_lines);

    Eigen::VectorXd final_reduced_f(2*FEM_nodeTags.size() - number_removed_lines);
    for (std::size_t j = 0; j < reduced_f_vector.size(); ++j){
        final_reduced_f(j) = reduced_f_vector[j];
    }

    // solve Kd=f
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(reduced_K_matrix);

    final_reduced_d = solver.solve(final_reduced_f);

    //double start4 = omp_get_wtime();
    //std::cout << "time to solve system: " << start4 -start3 << "\n";


    /*-----------RECONSTRUCTING FULL NODAL VALUES--------------*/
    std::vector<double> full_d_vector(2*FEM_nodeTags.size());
    for (std::size_t j = 0; j < full_d_vector.size(); ++j){
        if(new_DOFindices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = final_reduced_d(new_DOFindices[j]);
        }
    }

    /*--------------COMPUTING NODAL FORCES-------------*/
    Eigen::VectorXd eigen_d(2*FEM_nodeTags.size());
    for(size_t i = 0 ; i < 2*FEM_nodeTags.size(); ++i){
        eigen_d(i) = full_d_vector[i];
    }
    Eigen::VectorXd nodalForces(2*FEM_nodeTags.size());
    nodalForces = full_K_matrix*eigen_d;

    std::vector<double> nodal_forces(2*FEM_nodeTags.size());
    for(size_t i = 0; i < 2*FEM_nodeTags.size(); i++)
    {
        nodal_forces[i] = nodalForces(i);
    }

    /*--------------REACTION FORCES COMPUTATION------------------*/
    computeReactionForces(groups, nodeIndexMap, nodal_forces, keys);

    //double start5 = omp_get_wtime();
    //std::cout << "time for rest: " << start5 - start4 << "\n";

    /*--------------PLOTTING AND SAVING NODAL VALUES------------------*/
    std::vector<std::string> names;
    gmsh::model::list(names);
    plotUxUy(FEM_nodeTags, full_d_vector, nbViews, names);

    /*-------------POST PROCESSING---------------*/
    /*-------------Computing and saving strains/stresses-----------------*/

    groupname = "FEM_domain"; // we consider the domain as containing the FEM 2d model
    dim = groups[groupname].first;
    tag = groups[groupname].second;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elnodeTags;
    
    // compute total number of elements in the domain
    int nbelems=0;
    for(std::size_t i=0; i< tags.size(); ++i)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        for(std::size_t j=0; j< elementTags.size(); ++j)
        {
            nbelems+=elementTags[j].size();
        }
    }


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
    
    // computing the elemental nodal data for each physical group of the domain
    int k = 0; // current index of data vectors, k is incremented in Compute_Stress_Strain_Elemental_Nodal at each element
    for(std::size_t i=0; i< tags.size(); ++i)
    {
        gmsh::model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, tags[i]);
        computeStressStrainLinearWay(elementTypes, elementTags, elnodeTags, full_d_vector, elems2D, H_matrix, 
                                    nu, E, nodeIndexMap, data3, data4, data5, data6, data7, data8, data9, data10, k);
    }

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
    int nodeIndex;
    std::vector<std::vector<double>> data11(FEM_nodeTags.size());
    std::vector<std::vector<double>> data12(FEM_nodeTags.size());
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        nodeIndex = nodeIndexMap[FEM_nodeTags[i]];
        data11[i].resize(3);
        data12[i].resize(3);
        data11[i][0] = full_d_vector[nodeIndex]; //ux
        data11[i][1] = full_d_vector[nodeIndex+1]; //uy
        data11[i][2] = 0.; //uz
        data12[i][0] = nodal_forces[nodeIndex]; //fx
        data12[i][1] = nodal_forces[nodeIndex+1]; //fy
        data12[i][2] = 0.; //fz
    }

    gmsh::view::addModelData(viewtag11, 0, names[0], "NodeData",
                             FEM_nodeTags, data11); //independent of time here
    gmsh::view::addModelData(viewtag12, 0, names[0], "NodeData",
                             FEM_nodeTags, data12); //independent of time here



    // empeche que tout se plot l'un sur l'autre
    int nb_views = 12 + nbViews; //hardcoded
    for( int i = 0; i < nb_views; ++i){
        std::string my_string = "View[" + std::to_string(i) + "].Visible";
        gmsh::option::setNumber(my_string, 0);
    }

    //double start6 = omp_get_wtime();
    //std::cout << "time for saving results: " << start6 - start5 << "\n";


    if(openGUI)
        gmsh::fltk::run(); // opens the gmsh window
}
