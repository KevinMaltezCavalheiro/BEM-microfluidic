#include "functionsBEM.hpp"
#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include <omp.h>
#include <list>

int solverBEM(int argc, char **argv, std::map<int, double> &electrostaticPressure, std::map<int,std::pair<double,double>> & boundaryDisplacementMap, int & nbViews, bool postProcessing, const int iteration)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // Parameters for displaying variables (for debugging).
    bool dispNodeCoord = 0;  // Coordinates of the nodes.
    bool dispHierarchy = 0;  // Hierarchy (physical groups -> entities -> elements -> nodes).
    bool dispMatrices = 0;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
    bool dispValues = 0;     // Elements and their values of u and u_n.
    bool show_element_node_data = 0;
    bool show_time = 0; //display the elapsed time for the internal potential calculus in the terminal
    int visualisation = 1;  //0: no visualisation, 1: element node data, 2: element data
    bool convergence = false;

    //to display the time of the different functions
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::time_point<std::chrono::system_clock> start_total, end_total;
    std::string function_name;
    start_total = std::chrono::system_clock::now();

    // Initialises a vector of elements (that will contain all the elements of interrest)
    std::vector<elementStruct> elementVector;
    std::string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.

    // Gets the tags, the coordinates and the parametric coordinates of the nodes.
    std::map<size_t, int> allNodeIndex; // nodeTag -> position in allNodeTags and allCoord
    std::vector<size_t> allNodeTags;
    std::vector<double> allCoord;
    std::vector<double> allParametricCoord;
    gmsh::model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);

    for(size_t i = 0; i < allNodeTags.size(); i++)
    {
        allNodeIndex[allNodeTags[i]] = i;
    }

    // Print the coordinates of all nodes if dispNodeCoord != 0
    if(dispNodeCoord)
    {
        dispNodeCoordFun(allNodeTags, allCoord);
    }

    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    int dirichletCount = 0;

    //pour compute_internal_phi
    std::vector<size_t> NodeTags2d_domain;
    std::vector<double> NodeCoord2d_domain;
    //    vector<double> NodeParametricCoord2d_domain;
    //    model::mesh::getNodes(NodeTags2d_domain, NodeCoord2d_domain, NodeParametricCoord2d_domain,2);

    //pour elementnodedata
    std::vector<int> elementTypes2d_domain_BEM;
    std::vector<std::vector<size_t> > elementTags2d_domain_BEM;
    std::vector<std::vector<size_t> > nodeTags2d_domain_BEM;
    int entities_tags_for_physical_group_domain_BEM;
    std::vector<int> entities_tags_2d_domain_BEM;

    //  model::mesh::getElements(elementTypes2d,elementTags2d,nodeTags2d,2);

    start = std::chrono::system_clock::now();
    //calculates the potential at the barycentre of each element of the domain
    getBC(allNodeIndex, allCoord, elementVector, dirichletCount, dispHierarchy, NodeTags2d_domain, NodeCoord2d_domain,elementTypes2d_domain_BEM,elementTags2d_domain_BEM,nodeTags2d_domain_BEM,entities_tags_for_physical_group_domain_BEM, entities_tags_2d_domain_BEM, boundaryDisplacementMap);
    end = std::chrono::system_clock::now();
    function_name = "getBC";
    display_time(start, end, function_name, show_time);

    std::vector<std::string> names;
    gmsh::model::list(names);
    //std::cout << "the file contains " << names.size() << " model names\n";

    // Declaration of matrices A, B and vector c (the matrix type is perhaps not the most efficient).
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(elementVector.size(), 1);

    //  Double loop over each element in order to compute the element h and g that will be given to the matrices A and B.
    start = std::chrono::system_clock::now();
    double h;
    double g;
#pragma omp parallel for schedule(guided) shared(elementVector,A,B,c,integrationType,dirichletCount) private(h,g)
    for(int i = 0; i < (int) elementVector.size(); i++)
    {
        for(int j = 0; j < (int) elementVector.size(); j++)
        {
            // Will be computed and stored in A and B.
            h = 0;
            g = 0;

            if(i == j)  // Formula for diagonal elements.
            {
                h = -0.5;
                g = (elementVector[i].length/(2*M_PI)) * (log(elementVector[i].length/2)-1);
            }
            else    // General case
            {
                h = computeH(elementVector[i].midX, elementVector[i].midY, elementVector[j]);
                g = computeG(elementVector[i].midX, elementVector[i].midY, elementVector[j], integrationType);
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
    end = std::chrono::system_clock::now();
    function_name = "compute H and G";
    display_time(start, end, function_name, show_time);


    //If we want to display the time it takes for eigen to solve the linear system
    start = std::chrono::system_clock::now();
    Eigen::MatrixXd b = B*c;
    Eigen::MatrixXd x = A.fullPivLu().solve(b);
    end = std::chrono::system_clock::now();
    function_name = "solve linear system";
    display_time(start, end, function_name, show_time);

    // Matrices verification
    if(dispMatrices)
    {
        dispMatricesFun(A, B, c, b, x);
    }

    // Assign each value of x to the right element.
    for(int i = 0; i < (int) elementVector.size(); i++)
    {
        if( i < dirichletCount)
            elementVector[i].bcValue[1] = x(i);
        else
            elementVector[i].bcValue[0] = x(i);
    }

    if(dispValues)
    {
        dispValuesFun(elementVector);
    }

    /*
     *  Visualisation
     */

    //visualisation by elementnodedata
    if(visualisation == 1 && postProcessing)
    {
        //The first index indicates the elementType, the second the elementTag and the third the nodesTags associated with this element. For example, it says whether it is a square or a triangle, the element tag and finally the tag of the 3 or 4 nodes associated with it.
        std::vector<std::vector<std::vector<double>>> data;
        //This map associates the value of the potential to all nodes strictly included in the domain and not on the boundary.
        std::map<int, double> nodeMapDomain;

        //To displays the time of the this code section
        //the purpose of this function is to fill the nodemapdomain map and thus to calculate the value of the potential at each node in the domain
        start = std::chrono::system_clock::now();
        //calculates the potential at the barycentre of each element of the domain
        compute_internal_phi(allCoord,NodeTags2d_domain,allNodeIndex,elementVector,nodeMapDomain,integrationType);
        end = std::chrono::system_clock::now();
        function_name = "compute_internal_phi";
        display_time(start, end, function_name, show_time);

        //this function allows to correctly associate the value of the potential to each node of an element in order to display the ElementNodeData view
        start = std::chrono::system_clock::now();
        element_node_data(data,elementTags2d_domain_BEM,elementTypes2d_domain_BEM,nodeTags2d_domain_BEM,elementVector, nodeMapDomain,show_element_node_data,entities_tags_for_physical_group_domain_BEM);
        end = std::chrono::system_clock::now();
        std::string function_name = "ElementNodeData";
        display_time(start, end, function_name, show_time);

        int viewNodeTag = gmsh::view::add("phi");

        for(size_t i = 0; i < data.size(); i++)
        {
            gmsh::view::addModelData(viewNodeTag,0,names[0],"ElementNodeData",elementTags2d_domain_BEM[i],data[i],0.0,1);
        }
        gmsh::view::write(viewNodeTag, "results.msh");

        gmsh::option::setNumber("View[0].IntervalsType", 3);
        //option::setNumber("View[0].AdaptVisualizationGrid", 1);
        gmsh::option::setNumber("View[0].MaxRecursionLevel", 3);
        gmsh::option::setNumber("View[0].TargetError", -0.0001);
    }
    //visualisation by element barycenter
    if(visualisation == 2 && postProcessing)
    {
        std::vector<std::string> names;
        gmsh::model::list(names);

        //To get the total number of 2d elements that make up the domain, we simply add the number of squares to the number of triangles,...
        int nb2dElements = 0;
        for(size_t i = 0; i < elementTypes2d_domain_BEM.size(); i++)
            nb2dElements += elementTags2d_domain_BEM[i].size();

        //Tag and associated potential of all 2d elements. In this view the potential is calculated at the barycentre and is constant per element.
        std::vector<std::vector<double>> phi(nb2dElements,std::vector<double>(1,0));
        std::vector<size_t> tags(nb2dElements);

        //To displays the time of the this code section
        start = std::chrono::system_clock::now();
        //calculates the potential at the barycentre of each element of the domain
        computeData(elementVector, tags, phi, integrationType, entities_tags_2d_domain_BEM);
        end = std::chrono::system_clock::now();
        function_name = "computeData";
        display_time(start, end, function_name, show_time);

        /*
         *  Convergence
         */
        
        if(convergence)
        {
            double l = 1.0;
            double phiLeft = 100.0;
            double phiRight = 200.0;
            
            double error = 0.0;

            std::vector<double> barycenters;

            for(std::size_t i = 0; i < entities_tags_2d_domain_BEM.size(); i++)
            {
                std::vector<int> elementTypes_for_entity;
                std::vector<std::vector<size_t> > elementTags_for_entity;
                std::vector<std::vector<size_t> > nodeTags_for_entity;
                gmsh::model::mesh::getElements(elementTypes_for_entity, elementTags_for_entity, nodeTags_for_entity, 2, entities_tags_2d_domain_BEM[i]);
                
                for(std::size_t j = 0; j < elementTypes_for_entity.size(); j++)
                {
                    gmsh::model::mesh::getBarycenters(elementTypes_for_entity[j], entities_tags_2d_domain_BEM[i], false, false, barycenters);
                    for(std::size_t k = 0; k < elementTags_for_entity[i].size(); k++)
                    {
                        std::size_t tag = elementTags_for_entity[i][k];
                        int index = 0;
                        for(std::size_t l = 0; l < tags.size(); l++)
                        {
                            if(tag == tags[l])
                            {
                                index = l;
                                break;
                            }
                        }
                        double numericalPhi = phi[index][0];

                        double x = barycenters[3*k];
                        double exactPhi = phiLeft + (phiRight - phiLeft) * x / l;

                        double relError = abs( (exactPhi - numericalPhi)/exactPhi );
                        error += relError*relError;
                    }
                }
            }
            error = sqrt(error/tags.size());
            std::cout << "\nNumber of elements = " << tags.size() << "\t-> Error = " << error << "\n\n";
        }

        //For ease of use, a function handles the code for the display.
        displayData("phi", names[0], "results.msh", tags, phi);
    }

    //Calculates the electrostatic pressure per element in order to provide these results to the FEM code that will handle the rest of the algorithm
    double epsilon = 8.8541878128e-12;
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        electrostaticPressure[elementVector[i].tag] = epsilon * pow(elementVector[i].bcValue[1], 2.0) / 2; // not forget factor 2
    }

    end_total = std::chrono::system_clock::now();
    function_name = "total time of the code BEM";
    display_time(start_total, end_total, function_name, show_time);

    nbViews = 1; // number of views generated (useful for FEM display of results)

    return 0;
}
