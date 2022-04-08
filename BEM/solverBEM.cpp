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

using namespace gmsh;
using namespace std;
using namespace Eigen;

int solverBEM(int argc, char **argv, map<int, double> &electrostaticPressure, int & nbViews)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // Parameters for displaying variables (for debugging).
    bool dispNodeCoord = 0;  // Coordinates of the nodes.
    bool dispHierarchy = 0;  // Hierarchy (physical groups -> entities -> elements -> nodes).
    bool dispMatrices = 0;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
    bool dispValues = 0;     // Elements and their values of u and u_n.
    bool openGUI = 0;        // Graphical interface of GMSH.
    bool show_element_node_data = 0;
    bool show_time = 1; //display the elapsed time for the internal potential calculus in the terminal
    int visualisation = 1;  //0: no visualisation, 1: element node data, 2: element data

    //to display the time of the different functions
    chrono::time_point<chrono::system_clock> start, end;
    chrono::time_point<chrono::system_clock> start_total, end_total;
    string function_name;
    start_total = chrono::system_clock::now();

    // Initialises a vector of elements (that will contain all the elements of interrest)
    vector<elementStruct> elementVector;
    string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.

    // Gets the tags, the coordinates and the parametric coordinates of the nodes.
    map<size_t, int> allNodeIndex; // nodeTag -> position in allNodeTags and allCoord
    vector<size_t> allNodeTags;
    vector<double> allCoord;
    vector<double> allParametricCoord;
    model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);

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
    vector<size_t> NodeTags2d_domain;
    vector<double> NodeCoord2d_domain;
    //    vector<double> NodeParametricCoord2d_domain;
    //    model::mesh::getNodes(NodeTags2d_domain, NodeCoord2d_domain, NodeParametricCoord2d_domain,2);

    //pour elementnodedata
    vector<int> elementTypes2d_domain_BEM;
    vector<vector<size_t> > elementTags2d_domain_BEM;
    vector<vector<size_t> > nodeTags2d_domain_BEM;
    int entities_tags_for_physical_group_domain_BEM;
    //  model::mesh::getElements(elementTypes2d,elementTags2d,nodeTags2d,2);

    start = chrono::system_clock::now();
    //calculates the potential at the barycentre of each element of the domain
    getBC(allNodeIndex,allCoord, elementVector, dirichletCount, dispHierarchy, NodeTags2d_domain, NodeCoord2d_domain,elementTypes2d_domain_BEM,elementTags2d_domain_BEM,nodeTags2d_domain_BEM,entities_tags_for_physical_group_domain_BEM);
    end = chrono::system_clock::now();
    function_name = "getBC";
    display_time(start, end, function_name, show_time);

    vector<string> names;
    model::list(names);
    cout << "the file contains " << names.size() << " model names\n";

    // Declaration of matrices A, B and vector c (the matrix type is perhaps not the most efficient).
    MatrixXd A = MatrixXd::Zero(elementVector.size(), elementVector.size());
    MatrixXd B = MatrixXd::Zero(elementVector.size(), elementVector.size());
    MatrixXd c = MatrixXd::Zero(elementVector.size(), 1);

    //  Double loop over each element in order to compute the element h and g that will be given to the matrices A and B.
    start = chrono::system_clock::now();
    double h;
    double g;
#pragma omp parallel for schedule(guided) shared(elementVector,A,B,c,integrationType,dirichletCount) private(h,g)
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        for(size_t j = 0; j < elementVector.size(); j++)
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
    end = chrono::system_clock::now();
    function_name = "compute H and G";
    display_time(start, end, function_name, show_time);


    //If we want to display the time it takes for eigen to solve the linear system
    start = chrono::system_clock::now();
    MatrixXd b = B*c;
    MatrixXd x = A.fullPivLu().solve(b);
    end = chrono::system_clock::now();
    function_name = "solve linear system";
    display_time(start, end, function_name, show_time);

    // Matrices verification
    if(dispMatrices)
    {
        dispMatricesFun(A, B, c, b, x);
    }

    // Assign each value of x to the right element.
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        if(i < dirichletCount)
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
    if(visualisation == 1)
    {
        //The first index indicates the elementType, the second the elementTag and the third the nodesTags associated with this element. For example, it says whether it is a square or a triangle, the element tag and finally the tag of the 3 or 4 nodes associated with it.
        vector<vector<vector<double>>> data;
        //This map associates the value of the potential to all nodes strictly included in the domain and not on the boundary.
        map<int, double> nodeMapDomain;

        //To displays the time of the this code section
        //the purpose of this function is to fill the nodemapdomain map and thus to calculate the value of the potential at each node in the domain
        start = chrono::system_clock::now();
        //calculates the potential at the barycentre of each element of the domain
        compute_internal_phi(allCoord,NodeTags2d_domain,elementVector,nodeMapDomain,integrationType);
        end = chrono::system_clock::now();
        function_name = "compute_internal_phi";
        display_time(start, end, function_name, show_time);

        //this function allows to correctly associate the value of the potential to each node of an element in order to display the ElementNodeData view
        start = chrono::system_clock::now();
        element_node_data(data,elementTags2d_domain_BEM,elementTypes2d_domain_BEM,nodeTags2d_domain_BEM,elementVector, nodeMapDomain,show_element_node_data,entities_tags_for_physical_group_domain_BEM);
        end = chrono::system_clock::now();
        string function_name = "ElementNodeData";
        display_time(start, end, function_name, show_time);

        int viewNodeTag = view::add("ElementNodeData");

        for(size_t i = 0; i < data.size(); i++)
        {
            view::addModelData(viewNodeTag,0,names[0],"ElementNodeData",elementTags2d_domain_BEM[i],data[i],0.0,1);
        }
        view::write(viewNodeTag, "results.msh");

        option::setNumber("View[0].IntervalsType", 3);
        option::setNumber("View[0].AdaptVisualizationGrid", 1);
        option::setNumber("View[0].MaxRecursionLevel", 3);
        option::setNumber("View[0].TargetError", -0.0001);
    }
    //visualisation by element barycenter
    if(visualisation == 2)
    {
        vector<string> names;
        model::list(names);

        //To get the total number of 2d elements that make up the domain, we simply add the number of squares to the number of triangles,...
        int nb2dElements = 0;
        for(size_t i = 0; i < elementTypes2d_domain_BEM.size(); i++)
            nb2dElements += elementTags2d_domain_BEM[i].size();

        //Tag and associated potential of all 2d elements. In this view the potential is calculated at the barycentre and is constant per element.
        vector<vector<double>> phi(nb2dElements);
        vector<size_t> tags(nb2dElements);

        //To displays the time of the this code section
        start = chrono::system_clock::now();
        //calculates the potential at the barycentre of each element of the domain
        computeData(elementVector, elementTypes2d_domain_BEM, elementTags2d_domain_BEM, tags, phi, integrationType);
        end = chrono::system_clock::now();
        function_name = "computeData";
        display_time(start, end, function_name, show_time);

        //For ease of use, a function handles the code for the display.
        displayData("phi", names[0], "results.msh", tags, phi);


        /*---------ajout Louis--------*/
        //        vector<string> names;
        //        model::list(names);
        //        //vector<vector<double>> phi(nb2dElements);
        //        //vector<size_t> tags(nb2dElements);
        //        vector<vector<double>> phi;
        //        vector<size_t> tags;
        //
        //        // get physical groups
        //        map<string, pair<int, int>> groups;
        //        vectorpair dimTags;
        //        model::getPhysicalGroups(dimTags);
        //        for (size_t i = 0; i < dimTags.size(); ++i)
        //        {
        //            int dim = dimTags[i].first;
        //            int tag = dimTags[i].second;
        //            string name;
        //            model::getPhysicalName(dim, tag, name);
        //            groups[name] = {dim, tag};
        //        }
        //
        //        // get all the 2D elements in the BEM domain
        //        string groupname = "BEM_domain";
        //        int dim = groups[groupname].first;
        //        int tag = groups[groupname].second;
        //        vector<int> entities_tags;
        //        model::getEntitiesForPhysicalGroup(dim, tag, entities_tags);
        //        vector<int> elementTypes;
        //        vector<vector<size_t>> elementTags;
        //        vector<vector<size_t>> elnodeTags;
        //
        //
        //        //To displays the time of the this code section
        //        start = chrono::system_clock::now();
        //        //calculates the potential at the barycentre of each element of the domain
        //        int nb2dElements = 0;
        //        for(size_t i=0; i< entities_tags.size(); ++i)
        //        {
        //            model::mesh::getElements(elementTypes, elementTags, elnodeTags, dim, entities_tags[i]);
        //            for(size_t j=0; j< elementTags.size(); ++j)
        //            {
        //                nb2dElements += elementTags[j].size();
        //            }
        //            computeData(elementVector, elementTypes, elementTags, tags, phi, integrationType);
        //        }
        //        end = chrono::system_clock::now();
        //        function_name = "computeData";
        //        display_time(start, end, function_name, show_time);
        //
        //        displayData("phi", names[0], "results.msh", tags, phi);
    }

    //Calculates the electrostatic pressure per element in order to provide these results to the FEM code that will handle the rest of the algorithm
    double epsilon = 8.8541878128e-12;
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        electrostaticPressure[elementVector[i].tag] = epsilon * pow(elementVector[i].bcValue[1], 2.0);
    }

    end_total = chrono::system_clock::now();
    function_name = "total time of the code BEM";
    display_time(start_total, end_total, function_name, show_time);
    if(openGUI)
    {
        fltk::run();
    }

    nbViews = 1; // number of views generated (useful for FEM display of results)

    return 0;
}
