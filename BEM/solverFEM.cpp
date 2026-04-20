#include "functionsFEM.hpp"

#ifdef _MSC_VER
#include <gmsh.h_cwrap>
#else
#include <gmsh.h>
#endif

#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Core>
#include <cmath>
#include <algorithm>
#include <omp.h>

// solves the non-linear (large displacements, no large deformations) elastic problem.
// if coupled to the BEM solver, "electrostaticPressure" is the pressure computed from the electric field
// on the BEM_FEM boundary (mapping element tag -> pressure). "boundaryDisplacementMap" is filled for all
// nodes on the BEM_FEM_boundary as (node tag -> (u,v) displacement)
// "nbViews" is the current number of views already allocated in the gmsh window.
// "postProcessing" boolean activates the post processing operations (must only be computed once at the end of the
// coupled solver). "iteration" denotes the current iteration of the two-way coupled solver (starting at 1).
// "viewTagU" and "viewTagF" correspond to the view tags of the nodal displacements and the nodal forces in the gmsh window.
// when "untangleMesh" is activated, the nodal coordinates of the nodes in the FEM_domain are updated according to
// their final nodal displacements.
void solverFEMnonLinear(std::map<int, double> &electrostaticPressure, 
                    std::map<int,std::pair<double,double>> & boundaryDisplacementMap, int &nbViews, bool postProcessing, 
                    const int iteration, const int viewTagU, const int viewTagF, bool untangleMesh)
{
    /*--------------INITIALIZATION OF ITERATIVE ALGORITHM----------------*/
    /* FIRST STEP: 
    - construct and compute the vector of applied forces
    - set up the kinematics of each element and compute the local stiffness matrix
    - initialize the full d vector with Dirichlet boundary conditions and set pg = 0
    - important remark: in pg we focus only on the free DOFs (we do not focus on the DOFs fixed by dirichlet BCs)*/

    bool display_time = false; //displays time of every different step of the algorithm
    bool validation = false; //retrieves displacement of point A for large_rotation_validation.geo

    double start = omp_get_wtime();

    /*--------Integration rule used in the code for surface integration (neumann BC)----------*/
    std::string integration_rule = "CompositeGauss4"; 

    /*--------get all physical groups--------------*/
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
    std::string groupname = "FEM_domain"; // we consider the FEM_domain as containing the FEM 2d model
    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    // Getting entities of the domain to fill FEM_nodeTags
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

    std::vector<std::size_t> FEM_nodeTags;
    std::map<int,std::pair<double,double>> FEM_nodeCoordsMap; // nodeTag -> (x,y)
    RetrieveFEMNodeTagsAndCoords(FEM_nodeTags, FEM_nodeCoordsMap, dim, tags);

    /*-----get elements related to FEM--------------*/
    std::vector<std::size_t> FEM_elementTags;
    std::map<int,elementData> FEM_kinematicsMap; // map between elementTag and all the kinematics associated to this element
    // map is initialized, see functionsFEM.hpp for structure of elementData type
    int nbelems = initKinematics(FEM_elementTags, FEM_kinematicsMap, FEM_nodeCoordsMap, dim, tags);

    if(iteration == 1)
    {
        std::cout << "The FEM domain contains " << nbelems << " element(s) and " << FEM_nodeTags.size() << " nodes\n";
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
    double E = 0; // Young modulus
    double nu = 0; // Poisson ratio
    double rho = 0; // volumic mass
    double bx = 0; // volumic force along x axis
    double by = 0; // volumic force along y axis
    getPhysicalProperties(E, nu, rho, bx, by);

    /*-----------------H matrix involved in stress computation------------------*/
    Eigen::Matrix<double, 3, 3> H_matrix = computeHmatrix(E, nu);
    double start2 = omp_get_wtime();
    std::string step_name = "initialization of global variables";
    displayTime(start, start2, step_name, display_time);

    /*------------computing full f vector------------------------*/
    // we first consider volumic part of the f vector, computed simultaneously with local K matrices
    std::vector<double> full_f_vector(2*FEM_nodeTags.size());
    for (size_t i = 0; i < full_f_vector.size(); i++)
    {
        full_f_vector[i] = 0;
    }

    // Looping over entities of FEM_dimain and computing local K matrices and volumic f vector
    fillElementalKandF(H_matrix, rho, bx, by, dim, tags, nodeIndexMap, full_f_vector, FEM_kinematicsMap);

    double start3 = omp_get_wtime();
    step_name = "volumic f and elemental Kl";
    displayTime(start2, start3, step_name, display_time);


    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    std::vector<std::string> BCkeys;
    gmsh::onelab::getNames(BCkeys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    includeNeumannBC(BCkeys, integration_rule, nodeIndexMap, groups, full_f_vector);

    /* taking ELECTROSTATIC PRESSURE into account */
    if(electrostaticPressure.size() != 0) // size == 0 if FEM solver used alone
        applyElecPressure(integration_rule, FEM_nodeCoordsMap, nodeIndexMap, groups, electrostaticPressure, full_f_vector);
    double start4 = omp_get_wtime();
    step_name = "Neumann BC and elec pressure";
    displayTime(start3, start4, step_name, display_time);

    /*------------FOCUS ON DIRICHLET BC'S---------------*/
    // new_DOFindices vector will contain the new indices of the K matrix after the rows/columns corresponding
    // to dirichlet boundary conditions have been removed
    // IN THIS IMPLEMENTATION WE DECIDE TO REMOVE EVERY ROW/COL ASSOCIATED TO DIR. BC's
    // non-homo. dir BC will induce a correction in the RHS term f of the other rows, while homogeneous not,
    // if the given index has been removed, it will be set to -1 in new_DOFindices
    std::vector<int> new_DOFindices(2*FEM_nodeTags.size());

    // dir_BC[index] will contain the dirichlet boundary condition prescribed to this node
    std::vector<double> dir_BC(2*FEM_nodeTags.size());

    // Fills the two vectors defined just above
    int number_removed_lines = 0;
    fillDOFindicesDirBC(BCkeys, nodeIndexMap, groups, new_DOFindices, dir_BC, number_removed_lines);

    // initializing the full_d_vector with the dirichlet boundary conditions
    std::vector<double> full_d_vector(2*FEM_nodeTags.size());
    for (std::size_t j = 0; j < full_d_vector.size(); ++j){
        if(new_DOFindices[j] < 0){ // line was previously removed due to dirichlet boundary conditions
            full_d_vector[j] = dir_BC[j];
        }
        else{
            full_d_vector[j] = 0;
        }
    }
    // initializing the applied forces vector (will remain constant) and the global p vector (p_g)
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

    // useful during the iterative process
    std::vector<double> global_f_vector(nb_free_DOFs); // f_g
    Eigen::Matrix<double, Eigen::Dynamic, 1> residual(nb_free_DOFs, 1);

    // useful for the stopping criterion
    std::vector<double> previousTotalNodalForces(2*FEM_nodeTags.size());
    std::vector<double> currentTotalNodalForces(2*FEM_nodeTags.size());
    std::vector<double> relativeDifference(2*FEM_nodeTags.size());
    double relativeForces = 1e17;
    double previousRelativeForces = 1e17;
    double currentMaxNodalForce;
    for(size_t i = 0; i < 2*FEM_nodeTags.size(); i++)
    {
        previousTotalNodalForces[i] = 0;
        currentTotalNodalForces[i] = 0;
    }
    double residualTolerance = 1e-10; // empirical

    double start5 = omp_get_wtime();
    step_name = "further initialization of global variables for the iterative algorithm";
    displayTime(start4, start5, step_name, display_time);


    /*---------ITERATIVE PART---------------*/
    int max_number_of_steps = 200; // can be tuned
    int step = 0;

    double num_diverged_steps = 0; // stores the number of times the residual has increased between two steps, after 5 times,
    // we take the assumption it has converged towards an other value
    for(step = 0; step < max_number_of_steps; step++)
    {
        /*------------STEP 2 of the algorithm----------------*/
        /* 
        - update the kinematics from last iteration (or from dirichlet boundary conditions)
        - compute the global nodal forces f_g based on previous results
        - compute the residual = norm(f_g - f_app)*/

        /*--------- UPDATING KINEMATICS---------------*/
        double start_update = omp_get_wtime();
        updateKinematics(FEM_kinematicsMap, full_d_vector, nodeIndexMap, FEM_elementTags);
        double end_update = omp_get_wtime();
        step_name = "update kinematics";
        displayTime(start_update, end_update, step_name, display_time);


        // stopping criterion: focussing on the iterative difference of the full nodal forces vector
        for(size_t i = 0; i < 2*FEM_nodeTags.size(); i++)
        {
            previousTotalNodalForces[i] = currentTotalNodalForces[i];
        }
        retrieveTotalNodalForces(currentTotalNodalForces, FEM_kinematicsMap, nodeIndexMap, FEM_elementTags);
        double end_retrieveTotalForces = omp_get_wtime();
        step_name = "retrieve total nodal forces";
        displayTime(end_update, end_retrieveTotalForces, step_name, display_time);

        if(step > 1)
        {
            previousRelativeForces = relativeForces;
            relativeForces = 0;
            currentMaxNodalForce = 0;
            for(size_t i = 0; i < 2*FEM_nodeTags.size(); i++)
            {
                if(previousTotalNodalForces[i] != 0.0)
                    relativeDifference[i] = abs( (currentTotalNodalForces[i] - previousTotalNodalForces[i]));
                else
                    relativeDifference[i] = 0;
                
                if(abs(currentTotalNodalForces[i]) > currentMaxNodalForce)
                    currentMaxNodalForce = abs(currentTotalNodalForces[i]);

                relativeForces += relativeDifference[i]*relativeDifference[i];
            }
            if(currentMaxNodalForce > 0)
                relativeForces = sqrt(relativeForces/(2*FEM_nodeTags.size())/currentMaxNodalForce); // norm made relative
            else //to avoid division by zero
                relativeForces = 1e17; 
            //std::cout << "step: " << step << ", relativeForces: " << relativeForces << "\n";

            if(previousRelativeForces < relativeForces) //the current step induces an increase of nodal difference
            {
                num_diverged_steps++;
                if(num_diverged_steps == 5)
                {
                    std::cout << "-----------------------------------------\n";
                    std::cout << "Pay attention, the differential residual has not reached the tolerance of " << residualTolerance << ", it has converged towards: " << relativeForces << ", after " << step <<" steps\n";
                }
            }
        }

        /*--------COMPUTING AND ASSEMBLING THE GLOBAL F VECTOR BASED ON PREVIOUS ITERATION-----------*/
        for(int i = 0; i < nb_free_DOFs; i++)
        {
            global_f_vector[i] = 0;
        }
        assembleGlobalF(global_f_vector, FEM_kinematicsMap, nodeIndexMap, new_DOFindices, FEM_elementTags);

        // computing the residual
        for(size_t i = 0; i < global_f_vector.size(); i++)
        {
            residual(i,0) = global_f_vector[i] - f_app_vector[i];
        }

        double end_assembleglobalf = omp_get_wtime();
        step_name = "assemble global f and residual norm";
        displayTime(end_retrieveTotalForces, end_assembleglobalf, step_name, display_time);

        /*------------STEP 3: EVALUATING THE RESIDUAL, have we converged ??----------*/
        if ((relativeForces < residualTolerance && step > 0) || num_diverged_steps == 5)
           break;

        /*------------STEP 4: ASSEMBLING THE GLOBAL TANGENT STIFFNESS MATRIX----------*/
        // it is assembled as a sparse matrix
        Eigen::SparseMatrix<double> global_tangent_matrix(nb_free_DOFs, nb_free_DOFs);
        assembleGlobalKg(global_tangent_matrix, FEM_kinematicsMap, nodeIndexMap, new_DOFindices, FEM_elementTags);

        double end_assembleglobalKg = omp_get_wtime();
        step_name = "assemble global stiffness matrix";
        displayTime(end_assembleglobalf, end_assembleglobalKg, step_name, display_time);

        //std::cout << "global tangent matrix: " << global_tangent_matrix << "\n";

        /*-----------STEP 5: SOLVING THE SYSTEM----------*/
        Eigen::Matrix<double, Eigen::Dynamic, 1> global_delta_p(nb_free_DOFs, 1);
        
        // solve Kg*Deltapg=-res, using a direct sparse Cholesky factorization,
        // efficient for sparse positive semi-definite matrix
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(global_tangent_matrix);

        global_delta_p = solver.solve(-residual);

        double end_solve = omp_get_wtime();
        step_name = "solving the system";
        displayTime(end_assembleglobalKg, end_solve, step_name, display_time);

        /*----------STEP 6: updating the global p vector and the full d vector of nodal values-----------*/
        for(size_t i = 0; i < global_p_vector.size(); i++)
        {
            global_p_vector[i] += global_delta_p(i,0);
        }
        updateNodalDispVector(full_d_vector, global_p_vector, new_DOFindices);

        double end_updateNodalDisp = omp_get_wtime();
        step_name = "updating the nodal disp vector";
        displayTime(end_solve, end_updateNodalDisp, step_name, display_time);
    }

    if(step == max_number_of_steps)
    {
        std::cout << "-----------------------------------------\n";
        std::cout << "Pay attention, the differential residual has not reached the tolerance of " << residualTolerance << ", final value: " << relativeForces << "\n";
    }

    double start_postPro = omp_get_wtime();
    
    //updating the kinematics according to the results of last iteration
    updateKinematics(FEM_kinematicsMap, full_d_vector, nodeIndexMap, FEM_elementTags);

    /*--------------COMPUTING NODAL FORCES-------------*/
    std::vector<double> final_nodal_forces(2*FEM_nodeTags.size());
    retrieveTotalNodalForces(final_nodal_forces, FEM_kinematicsMap, nodeIndexMap, FEM_elementTags);

    std::vector<std::string> names;
    gmsh::model::list(names);

    /*-------------POST PROCESSING---------------------*/
    if(postProcessing)
    {
        // retrieving and displaying maximal nodal displacement
        double max_disp = 0;
        double max_u = 0;
        double max_v = 0;
        int max_disp_nodeTag = 0;

        for(size_t i = 0; i < FEM_nodeTags.size(); i++)
        {
            double tmp_index = nodeIndexMap[FEM_nodeTags[i]];
            double tmp_disp = full_d_vector[tmp_index]*full_d_vector[tmp_index] + full_d_vector[tmp_index + 1]*full_d_vector[tmp_index + 1];
            if (tmp_disp > max_disp)
            {
                max_disp = tmp_disp;
                max_disp_nodeTag = FEM_nodeTags[i];
                max_u = full_d_vector[tmp_index];
                max_v = full_d_vector[tmp_index + 1];
            }
        }
        std::cout << "-----------------------------------------\n";
        std::cout << "maximal nodal displacement (u,v) = (" << max_u << "," << max_v << ") [m], at node " << max_disp_nodeTag << "\n";

        /*--------------REACTION FORCES COMPUTATION------------------*/
        computeReactionForces(groups, nodeIndexMap, final_nodal_forces, BCkeys);

        /*--------------PLOTTING NODAL VALUES------------------*/
        plotUxUy(FEM_nodeTags, full_d_vector, names);

        /*-------------POST PROCESSING---------------*/
        /*-------------Computing and saving elemental-nodal strains/stresses-----------------*/
        double max_VM_stress = 0;
        int max_VM_nodeTag;
        computeStressStrainElNodal(dim, tags, H_matrix, nu, E, FEM_kinematicsMap, nbelems, max_VM_stress, max_VM_nodeTag);
        std::cout << "-----------------------------------------\n";
        std::cout << "maximal nodal von Mises stress: " << max_VM_stress << " [Pa], at node " << max_VM_nodeTag << "\n";

        // retrieving work done by external forces
        double externalWork = 0;
        computeWorkDoneByExternalForces(externalWork, full_d_vector, final_nodal_forces);
        std::cout << "-----------------------------------------\n";
        std::cout << "work done by external forces: " << externalWork << " [J/m].\n";

        // strain energy
        double strainEnergy = 0;
        computeStrainEnergy(strainEnergy, dim, tags, FEM_kinematicsMap, H_matrix, integration_rule);
        // peut etre pas OK en non-linear -> repartir des local displacements pour calculer tout ça

        std::cout << "-----------------------------------------\n";
        std::cout << "strain energy: " << strainEnergy << " [J/m].\n";

        std::cout << "-----------------------------------------\n";
        std::cout << "Total potential energy: " << strainEnergy - externalWork << " [J/m].\n";
    }
    // representing vector displacement field and nodal forces as vectorial quantities
    std::vector<std::vector<double>> dataU(FEM_nodeTags.size());
    std::vector<std::vector<double>> dataF(FEM_nodeTags.size());
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        dataU[i].resize(3);
        dataF[i].resize(3);
        dataU[i][0] = full_d_vector[2*i]; //ux, here use NodeIndexMap
        dataU[i][1] = full_d_vector[2*i+1]; //uy
        dataU[i][2] = 0.; //uz
        dataF[i][0] = final_nodal_forces[2*i]; //fx
        dataF[i][1] = final_nodal_forces[2*i+1]; //fy
        dataF[i][2] = 0.; //fz
    }

    gmsh::view::addModelData(viewTagU, iteration, names[0], "NodeData",
                             FEM_nodeTags, dataU, iteration); 
    gmsh::view::addModelData(viewTagF, iteration, names[0], "NodeData",
                             FEM_nodeTags, dataF, iteration); 


    if(postProcessing)
    {
        nbViews = 12 + nbViews; // FEM code adds 12 views
        // avoiding that all the views are superposed
        for( int i = 0; i < nbViews; ++i){
            std::string my_string = "View[" + std::to_string(i) + "].Visible";
            gmsh::option::setNumber(my_string, 0);
        }
        gmsh::option::setNumber("View[0].VectorType", 5); // plot the displacement in deformed configuration
        gmsh::option::setNumber("View[0].RangeType", 3);
    }

    double end_postPro = omp_get_wtime();
    step_name = "post processing";
    displayTime(start_postPro, end_postPro, step_name, display_time);

    /*--------------RETRIEVING NODAL DISPLACEMENTS ON BEM-FEM BOUNDARY------------*/
    // those displacements will be sent to the BEM code in the non-linear iterative solver
    // mapping nodeTag (of node on the boundary) into (u,v) displacements of the node
    if(electrostaticPressure.size())
    {
        retrieveBoundaryDisplacements(boundaryDisplacementMap, full_d_vector, nodeIndexMap, groups);
    }

    // if mesh untangler is activated, all nodal coordinates are updated according to their final
    // nodal displacements. It is useful when coupled to the BEM solver in order to untangle
    // the BEM domain and compute the potential field in the whole BEM domain adapted to the
    // displacement of the FEM domain.
    if(untangleMesh)
    {
        for(size_t i = 0; i < FEM_nodeTags.size(); i++)
        {
            std::vector<double> coord, parametricCoord;
            int dim, tag;
            gmsh::model::mesh::getNode(FEM_nodeTags[i], coord, parametricCoord, dim, tag);
            std::vector<double> total_disp = {coord[0] + full_d_vector[2*i], coord[1] + full_d_vector[2*i + 1], coord[2]};
            gmsh::model::mesh::setNode(FEM_nodeTags[i], total_disp, {});
        }
    }

    // it is only useful when using "large_rotation_validation.geo" file.
    // printing the displacement of "point A" (see report for graph).
    if(validation)
    {
        std::string groupname = "point_A";
        int dim = groups[groupname].first;
        int tag = groups[groupname].second;
        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
        std::vector<std::size_t> A_nodeTag;
        std::vector<double> A_nodecoord;
        std::vector<double> A_nodeparametricCoord;
        gmsh::model::mesh::getNodes(A_nodeTag, A_nodecoord, A_nodeparametricCoord, dim, tags[0], true); 
        std::cout << "-----------------------------------------\n";
        std::cout << "point A displacement: u = " << full_d_vector[nodeIndexMap[A_nodeTag[0]]] << ", v = " << full_d_vector[nodeIndexMap[A_nodeTag[0]] + 1] << "\n";
    }
}

// solves the linear elastic problem in the FEM domain.
// "electrostaticPressure" contains the pressure resulting from the electric field computation on the 
// BEM_FEM_boundary. "nbViews" is the current number of views already allocated in the gmsh window.
void solverFEM(std::map<int, double> &electrostaticPressure, int &nbViews)
{
    /*--------Integration rule used in the code for surface integration (neumann BC)----------*/
    std::string integration_rule = "CompositeGauss4";

    bool display_time = false; //displays the time performance of the different parts of the code

    double start_init = omp_get_wtime();

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

    // Getting entities of the domain to fill FEM_nodeTags
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

    std::vector<std::size_t> FEM_nodeTags;
    std::map<int,std::pair<double,double>> FEM_nodeCoordsMap; // nodeTag -> (x,y) (useless in this linear case)
    RetrieveFEMNodeTagsAndCoords(FEM_nodeTags, FEM_nodeCoordsMap, dim, tags);

    /*-------map between nodeTag and index of first nodal displacement associated to the node (second one is the same +1)------*/
    std::map<int, int> nodeIndexMap;
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        int nodeTag = FEM_nodeTags[i];
        int node_first_index = 2*i; 
        nodeIndexMap[nodeTag] = node_first_index;
    }

    /*-----------get physical properties of domain----------------------*/
    double E = 0; // Young modulus
    double nu = 0; // Poisson ratio
    double rho = 0; // volumic mass
    double bx = 0; // volumic force along x axis
    double by = 0; // volumic force along y axis
    getPhysicalProperties(E, nu, rho, bx, by);

    /*-----------------H matrix involved in stress computation------------------*/
    Eigen::Matrix<double, 3, 3> H_matrix = computeHmatrix(E, nu);

    double end_init = omp_get_wtime();
    std::string step_name = "initialization";
    displayTime(start_init, end_init, step_name, display_time);

    /*------------computing K matrix and f vector------------------------*/
    // we first consider the K matrix and the volumic part of the f vector
    // K matrix declared as sparse matrix for efficiency.
    Eigen::SparseMatrix<double> full_K_matrix(2*FEM_nodeTags.size(),2*FEM_nodeTags.size());
    std::vector<double> full_f_vector(2*FEM_nodeTags.size());
    for (size_t i = 0; i < full_f_vector.size(); i++)
    {
        full_f_vector[i] = 0;
    }

    // Looping over entities and compute full K and volumic F
    assembleKVolumicF(H_matrix, rho, bx, by, dim, tags, nodeIndexMap, full_K_matrix, full_f_vector);

    double end_assembleK = omp_get_wtime();
    step_name = "K assembly";
    displayTime(end_init, end_assembleK, step_name, display_time);
    

    /*-----------------BOUNDARY CONDITIONS------------*/
    /*-----------get group_names related to BC's----------------------*/
    std::vector<std::string> BCkeys;
    gmsh::onelab::getNames(BCkeys, "(Boundary Conditions).+");
    /*-----------first focus on neumann BC's (surface traction) into full_f_vector------------*/
    includeNeumannBC(BCkeys, integration_rule, nodeIndexMap, groups, full_f_vector);

    /* ELECTROSTATIC PRESSURE */
    if(electrostaticPressure.size() != 0) // size is only greater than zero if there was already a BEM solver
    {
        applyElecPressure(integration_rule, FEM_nodeCoordsMap, nodeIndexMap, groups, electrostaticPressure, full_f_vector);
    }

    double end_neumannBC = omp_get_wtime();
    step_name = "neumann BC and elec pressure";
    displayTime(end_assembleK, end_neumannBC, step_name, display_time);

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
    int number_removed_lines = 0;
    fillDOFindicesDirBC(BCkeys, nodeIndexMap, groups, new_DOFindices, dir_BC, number_removed_lines);

    double end_dirichletBC = omp_get_wtime();
    step_name = "dirichlet BC";
    displayTime(end_neumannBC, end_dirichletBC, step_name, display_time);

    /*------------FILLING THE REDUCED K MATRIX AND F VECTOR--------------*/
    //the correction term on the f vector comes from the "Finite Element Method" theoretical course (JP Ponthot, ULiege), Chap2 Slide 38/50
    Eigen::SparseMatrix<double> reduced_K_matrix(2*FEM_nodeTags.size() - number_removed_lines,2*FEM_nodeTags.size() - number_removed_lines);
    std::vector<double> reduced_f_vector(2*FEM_nodeTags.size() - number_removed_lines);

    fillReducedKf(full_f_vector, full_K_matrix, new_DOFindices, dir_BC, reduced_K_matrix, reduced_f_vector);

    double end_reducedFilling = omp_get_wtime();
    step_name = "filling reduced K and f";
    displayTime(end_dirichletBC, end_reducedFilling, step_name, display_time);

    /*--------------SOLVING THE REDUCED SYSTEM--------------------*/
    Eigen::VectorXd final_reduced_d(2*FEM_nodeTags.size() - number_removed_lines);

    Eigen::VectorXd final_reduced_f(2*FEM_nodeTags.size() - number_removed_lines);
    for (std::size_t j = 0; j < reduced_f_vector.size(); ++j){
        final_reduced_f(j) = reduced_f_vector[j];
    }

    // solve Kd=f
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(reduced_K_matrix);

    final_reduced_d = solver.solve(final_reduced_f);

    double end_solveSystem = omp_get_wtime();
    step_name = "solving the linear system";
    displayTime(end_reducedFilling, end_solveSystem, step_name, display_time);

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

    double end_reconstructing = omp_get_wtime();
    step_name = "reconstructing full nodal disp and forces";
    displayTime(end_solveSystem, end_reconstructing, step_name, display_time);

    /*--------------PLOTTING AND SAVING NODAL VALUES------------------*/
    std::vector<std::string> names;
    gmsh::model::list(names);
    plotUxUy(FEM_nodeTags, full_d_vector, names);

    /*-------------POST PROCESSING---------------*/
    /*-------------Computing and saving strains/stresses-----------------*/

    groupname = "FEM_domain"; // we consider the FEM domain as containing the FEM 2d model
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

    std::cout << "The FEM domain contains " << nbelems << " element(s) and " << FEM_nodeTags.size() << " nodes\n";

    // retrieving and displaying maximal nodal displacement
        double max_disp = 0;
        double max_u = 0;
        double max_v = 0;
        int max_disp_nodeTag = 0;

        for(size_t i = 0; i < FEM_nodeTags.size(); i++)
        {
            double tmp_index = nodeIndexMap[FEM_nodeTags[i]];
            double tmp_disp = full_d_vector[tmp_index]*full_d_vector[tmp_index] + full_d_vector[tmp_index + 1]*full_d_vector[tmp_index + 1];
            if (tmp_disp > max_disp)
            {
                max_disp = tmp_disp;
                max_disp_nodeTag = FEM_nodeTags[i];
                max_u = full_d_vector[tmp_index];
                max_v = full_d_vector[tmp_index + 1];
            }
        }
        std::cout << "-----------------------------------------\n";
        std::cout << "maximal nodal displacement (u,v) = (" << max_u << "," << max_v << ") [m], at node " << max_disp_nodeTag << "\n";

    /*--------------REACTION FORCES COMPUTATION------------------*/
    computeReactionForces(groups, nodeIndexMap, nodal_forces, BCkeys);

    /*-------------------ELEMENTAL-NODAL VALUES-----------------------*/
    double max_VM_stress = 0;
    int max_VM_nodeTag;
    computeStressStrainLinearWay(dim, tags, full_d_vector, H_matrix, nu, E, nodeIndexMap, nbelems, max_VM_stress, max_VM_nodeTag);
    std::cout << "-----------------------------------------\n";
    std::cout << "maximal nodal von Mises stress: " << max_VM_stress << " [Pa], at node " << max_VM_nodeTag << "\n";

    // representing vector displacement field and nodal forces
    int viewTagU = gmsh::view::add("u");
    int viewTagF = gmsh::view::add("f");
    int nodeIndex;
    std::vector<std::vector<double>> dataU(FEM_nodeTags.size());
    std::vector<std::vector<double>> dataF(FEM_nodeTags.size());
    for (std::size_t i = 0; i < FEM_nodeTags.size(); ++i)
    {
        nodeIndex = nodeIndexMap[FEM_nodeTags[i]];
        dataU[i].resize(3);
        dataF[i].resize(3);
        dataU[i][0] = full_d_vector[nodeIndex]; //ux
        dataU[i][1] = full_d_vector[nodeIndex+1]; //uy
        dataU[i][2] = 0.; //uz
        dataF[i][0] = nodal_forces[nodeIndex]; //fx
        dataF[i][1] = nodal_forces[nodeIndex+1]; //fy
        dataF[i][2] = 0.; //fz
    }

    gmsh::view::addModelData(viewTagU, 0, names[0], "NodeData",
                             FEM_nodeTags, dataU); //independent of time here
    gmsh::view::addModelData(viewTagF, 0, names[0], "NodeData",
                             FEM_nodeTags, dataF); //independent of time here

    int nb_views = 12 + nbViews; //FEM code adds 12 views
    for( int i = 0; i < nb_views; ++i){
        std::string my_string = "View[" + std::to_string(i) + "].Visible"; // avoids that all views are represented on top of each other
        gmsh::option::setNumber(my_string, 0);
    }

    // retrieving work done by external forces
    double externalWork = 0;
    computeWorkDoneByExternalForces(externalWork, full_d_vector, nodal_forces);
    std::cout << "-----------------------------------------\n";
    std::cout << "work done by external forces: " << externalWork << " [J/m].\n";

    // strain energy
    double strainEnergy = 0;
    computeStrainEnergyLinearWay(strainEnergy, dim, tags, full_d_vector, H_matrix, nodeIndexMap, integration_rule);

    std::cout << "-----------------------------------------\n";
    std::cout << "strain energy: " << strainEnergy << " [J/m].\n";

    std::cout << "-----------------------------------------\n";
    std::cout << "Total potential energy: " << strainEnergy - externalWork << " [J/m].\n";

    double end_postPro = omp_get_wtime();
    step_name = "post processing";
    displayTime(end_reconstructing, end_postPro, step_name, display_time);
}