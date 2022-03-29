#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include "code_kevin.hpp"

//todo: séparer en différentes fonctions
//todo: comparer au livre BEM
//todo: voir pour gérer des trous dans les surfaces
//todo: champs electrostatique
//todo: combiner avec la FEM
//todo: faire plein de tests différent (dans les .geo)
//todo: afficher et calculer la force electrostatique
//todo: test surface avec des trous
//todo: afficher champ en fonction de la distance
//todo: convergence
//todo: on a besoin de la pression electrostatique à la frontière mais on peut essayer de l'afficher partout dans le domaine
//todo: hybrid mesh

using namespace gmsh;
using namespace std;
using namespace Eigen;

int code_py(int argc, char **argv);

void show_mat(MatrixXd &A);

void show_vec(VectorXd &A);

void get_bc_from_onelab();

double compute_angle(double &deltaX1, double &deltaX2, double &deltaY1, double &deltaY2);

void get_bc_from_onelab(vector<elementStruct> &elementVector, int &dirichletCount, bool &dispHierarchy, vector<double> &allCoord);

void compute_g(double &jMidX, double &jDeltaX, double &jMidY, double &jDeltaY, double &iMidX, double &iMidY, double &jLength, double &g);

void compute_internal_phi(vector<double> &allCoord,vector<size_t> &NodeTags2d, vector<elementStruct> &elementVector, vector<vector<double>> &data_boundary_nodes, vector<vector<double>> &data2d, map<int,double> &nodeMapDomain);

void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d, vector<int> &elementTypes2d, vector<vector<size_t>> &nodeTags2d, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data);

int code_py(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // Parameters for displaying variables (for debugging).
    bool dispNodeCoord = 0;  // Coordinates of the nodes.
    bool dispHierarchy = 0;  // Hierarchy (physical groups -> entities -> elements -> nodes).
    bool dispAngles = 0;     // Angles used for the computation of h. /!\ Displays ~n^2 angles.
    bool dispMatrices = 0;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
    bool dispValues = 0;     // Elements and their values of u and u_n.
    bool openGUI = 1;        // Graphical interface of GMSH.
    bool show_potential = 0;
    bool show_element_node_data = 0;

    initialize();

    open(argv[1]);
    model::mesh::generate(2);

    // Initialises a vector of elements (that will contain all the elements of interrest)
    vector<elementStruct> elementVector;

    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    int dirichletCount = 0;

    // Gets the tags, the coordinates and the parametric coordinates of the nodes.
    vector<size_t> allNodeTags;
    vector<double> allCoord;
    vector<double> allParametricCoord;
    model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);

    // Print the coordinates of all nodes if dispNodeCoord != 0
    if(dispNodeCoord)
    {
        for(int i = 0; i < allNodeTags.size(); i++)
        {
            cout << "Node '" << allNodeTags[i] << "': (" << allCoord[3*i] << "; " << allCoord[3*i + 1] << "; " << allCoord[3*i + 2] << ")\n";
        }
        cout << endl;
    }

    vector<size_t> NodeTags2d_domain;
    vector<double> NodeCoord2d_domain;
    vector<double> NodeParametricCoord2d_domain;
    model::mesh::getNodes(NodeTags2d_domain, NodeCoord2d_domain, NodeParametricCoord2d_domain,2);

    vector<int> elementTypes_boundary;
    vector<vector<size_t> > elementTags_boundary;
    vector<vector<size_t> > nodeTags_boundary;
    model::mesh::getElements(elementTypes_boundary,elementTags_boundary,nodeTags_boundary,1);

    vector<string> names;
    model::list(names);
    cout << "the file contains " << names.size() << " model names\n";

    vector<int> elementTypes2d;
    vector<vector<size_t> > elementTags2d;
    vector<vector<size_t> > nodeTags2d;
    model::mesh::getElements(elementTypes2d,elementTags2d,nodeTags2d,2);

    get_bc_from_onelab(elementVector, dirichletCount, dispHierarchy, allCoord);

    int nb_element_boundary = elementVector.size();

    // Declaration of matrices A, B and vector c (the matrix type is perhaps not the most efficient).
    MatrixXd A = MatrixXd::Zero(nb_element_boundary, nb_element_boundary);
    MatrixXd B = MatrixXd::Zero(nb_element_boundary, nb_element_boundary);
    MatrixXd c = MatrixXd::Zero(nb_element_boundary, 1);

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
                double jMidX = (j1X + j2X) / 2.0;
                double jMidY = (j1Y + j2Y) / 2.0;
                double jDeltaX = j2X - j1X;
                double jDeltaY = j2Y - j1Y;
                double jLength = sqrt(pow(jDeltaX, 2.0) + pow(jDeltaY, 2.0));

                //compute the angle.
                double deltaX1 = j1X - iMidX;
                double deltaX2 = j2X - iMidX;
                double deltaY1 = j1Y - iMidY;
                double deltaY2 = j2Y - iMidY;

                double angle = compute_angle(deltaX1, deltaX2, deltaY1, deltaY2);

                // Angle verification
                if(dispAngles)
                {
                    double angleConst = 180.0 / M_PI;
                    cout << angle*angleConst << endl;
                }


                h = angle / (2*M_PI);

                compute_g(jMidX, jDeltaX, jMidY, jDeltaY, iMidX, iMidY, jLength, g);

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
    MatrixXd b = B*c;
    //MatrixXd x = A.householderQr().solve(b);
    MatrixXd x = A.fullPivLu().solve(b);

    // Matrices verification
    if(dispMatrices)
    {
        cout << "A matrix:\n";
        cout << A << endl;
        cout << "\nB matrix:\n";
        cout << B << endl;
        cout << "\nc vector:\n";
        cout << c << endl;
        cout << "\nb vector:\n";
        cout << b << endl;
        cout << "\nx vector:\n";
        cout << x << endl;
    }


    // Assign each value of x to the right element.
    for(int i = 0; i < elementVector.size(); i++)
    {
        if(i < dirichletCount)
            elementVector[i].bcValue[1] = x(i);
        else
            elementVector[i].bcValue[0] = x(i);
        if(dispValues)
        {
            cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
        }
    }

    /*
     *  Visualisation
     */

    vector<vector<vector<double>>> data;
    vector<vector<double>> data_boundary_nodes(allNodeTags.size(),vector<double> (1,0));
    vector<vector<double>> data2d(NodeTags2d_domain.size(),vector<double> (1,0));
    map<int, double> nodeMapDomain;
    double phi;
    vector<double> phi_value;
    for(int i = 0; i < elementVector.size(); i++)
    {
        data_boundary_nodes[elementVector[i].nodeTags[0]-1][0] = elementVector[i].bcValue[0];
    }

    compute_internal_phi(allCoord,NodeTags2d_domain,elementVector,data_boundary_nodes,data2d,nodeMapDomain);

    if(show_potential)
    {
        for(int i = 0; i < data_boundary_nodes.size(); i++)
        {
            cout << data_boundary_nodes[i][0] << endl;
        }
    }

    element_node_data(data,elementTags2d,elementTypes2d,nodeTags2d,elementVector, nodeMapDomain,show_element_node_data);

    int viewNodeTag = view::add("ElementNodeData");

    //view::addModelData(viewNodeTag, 0, names[0], "NodeData", allNodeTags, data_boundary_nodes, 0.0, 1);
    //view::addModelData(viewNodeTag, 0, names[0], "NodeData", NodeTags2d, data2d, 0.0, 1);
    for(int i = 0; i < data.size(); i++)
    {
        view::addModelData(viewNodeTag, 0, names[0], "ElementNodeData", elementTags2d[i], data[i], 0.0, 1);
    }
    view::write(viewNodeTag, "results.msh");

    option::setNumber("View[0].IntervalsType", 3);
    option::setNumber("View[0].AdaptVisualizationGrid", 1);
    option::setNumber("View[0].MaxRecursionLevel", 3);
    option::setNumber("View[0].TargetError", -0.0001);

    if(openGUI)
    {
        fltk::run();
    }

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

double compute_angle(double &deltaX1, double &deltaX2, double &deltaY1, double &deltaY2)
{
    double cosAlpha = (deltaX1*deltaX2 + deltaY1*deltaY2) / (sqrt((pow(deltaX1,2.0) + pow(deltaY1,2.0)) * (pow(deltaX2,2.0) + pow(deltaY2,2.0))));
    if(cosAlpha < -1.0)
    {
        cosAlpha = -1.0;
    }
    else if(cosAlpha > 1.0)
    {
        cosAlpha = 1.0;
    }

    double alpha = acos(cosAlpha);

    double crossProduct = deltaX1*deltaY2 - deltaY1*deltaX2;
    if(crossProduct < 0)
    {
        alpha = -alpha;
    }
    return alpha;
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

void get_bc_from_onelab(vector<elementStruct> &elementVector, int &dirichletCount, bool &dispHierarchy, vector<double> &allCoord)
{
    // Gets the dim and the tags of the physical groups
    vectorpair physicalGroupsDimTags;
    model::getPhysicalGroups(physicalGroupsDimTags, 1);

    // Gets the boundary conditions from onelab and create a map: "name of the physical group" -> ("dirichlet"/"neumann", value)
    map<string, pair<string, double>> bcMap;
    vector<string> bcKeys;
    onelab::getNames(bcKeys, "Boundary Conditions.+");
    for(auto &key : bcKeys)
    {
        vector<double> values;
        onelab::getNumber(key, values);

        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/'))
            words.push_back(word);

        if(words.size() == 3)
            bcMap[words[1]] = {words[2], values[0]};
    }

    for (int i = 0; i < physicalGroupsDimTags.size(); ++i)
    {
        int physicalGroupsDim = physicalGroupsDimTags[i].first;
        int physicalGroupsTag = physicalGroupsDimTags[i].second;
        string name;
        model::getPhysicalName(physicalGroupsDim, physicalGroupsTag, name);

        string bcType = bcMap[name].first;
        double bcValue = bcMap[name].second;

        // Prints the physical groups name, dim, tag, bcType and bcValue (if dispHierarchy != 0).
        if(dispHierarchy)
        {
            cout << "Physical group '" << name << "' (" << physicalGroupsDim << "D) with tag '" << physicalGroupsTag << "' has a '" << bcType << "' type and the value " << bcValue << "\n";
        }

        vector<int> entitiesTags;
        model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entitiesTags);

        for (int k = 0; k < entitiesTags.size(); k++)
        {
            // Print the tag of the entities of each physical group (if dispHierarchy != 0).
            if(dispHierarchy)
            {
                cout << "\t contains entity '" << entitiesTags[k] << "' ";
            }

            vector<size_t> elementTags, nodeTags;
            model::mesh::getElementsByType(1,elementTags,nodeTags,entitiesTags[k]);

            // Print the number of element of the given entity (if dispHierarchy != 0).
            if(dispHierarchy)
            {
                cout << "which has " << elementTags.size() << " elements:\n";
            }


            for(int l = 0; l < elementTags.size(); l++)
            {
                elementStruct element;
                element.tag = elementTags[l];
                element.bcType = bcType;
                element.nodeTags[0] = nodeTags[2*l];
                element.nodeTags[1] = nodeTags[2*l+1];

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
                    cout << "\t\t-Element '" << elementTags[l] << "' (" << element.bcType << "):\n";
                    cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*(element.nodeTags[0]-1)] << "; " << allCoord[3*(element.nodeTags[0]-1) + 1] << ")\n";
                    cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*(element.nodeTags[1]-1)] << "; " << allCoord[3*(element.nodeTags[1]-1) + 1] << ")\n";
                }
            }
            if(dispHierarchy)
            {
                cout << endl;
            }
        }
    }
}

void compute_g(double &jMidX, double &jDeltaX, double &jMidY, double &jDeltaY, double &iMidX, double &iMidY, double &jLength, double &g)
{
    // Computation of the integration points (parametric) and weights.
    string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.
    vector<double> localCoord;
    vector<double> weights;
    model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

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

void compute_internal_phi(vector<double> &allCoord,vector<size_t> &NodeTags2d_domain, vector<elementStruct> &elementVector, vector<vector<double>> &data_boundary_nodes, vector<vector<double>> &data2d, map<int,double> &nodeMapDomain)
{
    for(int i = 0; i < NodeTags2d_domain.size(); i++)
    {
        // Coordinates of the node
        double x = allCoord[3*(NodeTags2d_domain[i]-1)];
        double y = allCoord[3*(NodeTags2d_domain[i]-1) + 1];

        double phi = 0;
        for(int j = 0; j < elementVector.size(); j++)
        {
            double h;
            double g;

            // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element i here).
            int j1Tag = elementVector[j].nodeTags[0];
            int j2Tag = elementVector[j].nodeTags[1];

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

            // Compute the angle
            double deltaX1 = j1X - x;
            double deltaX2 = j2X - x;
            double deltaY1 = j1Y - y;
            double deltaY2 = j2Y - y;

            double angle = compute_angle(deltaX1, deltaX2, deltaY1, deltaY2);

            h = angle / (2*M_PI);

            // Computation of the integration points (parametric) and weights.
            string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.
            vector<double> localCoord;
            vector<double> weights;
            model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

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
        data_boundary_nodes[NodeTags2d_domain[i]-1][0] = phi;
        data2d[i][0] = phi;
        nodeMapDomain[NodeTags2d_domain[i]] = phi;
    }
}

void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d, vector<int> &elementTypes2d, vector<vector<size_t>> &nodeTags2d, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data)
{
    data.resize(elementTags2d.size());
    for(int i = 0; i < elementTags2d.size(); i++)
    {
        data[i].resize(elementTags2d[i].size());
        for(int j = 0; j < elementTags2d[i].size(); j++)
        {
            int nb_nodes_element2d = elementTypes2d[i]+1;
            vector<int> elements_commun_tag(0,0);
            vector<int> elements_commun_indice(0,0);
            vector<int> nodes_commun_tag(0,0);
            data[i][j].resize(nb_nodes_element2d,HUGE_VAL);
            for(int k = 0; k < nb_nodes_element2d; k++)
            {
                int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                for(int l = 0; l < elementVector.size(); l++)
                {
                    if(elementVector[l].nodeTags[0] == loop_element_2d_nodes_tag||elementVector[l].nodeTags[1] == loop_element_2d_nodes_tag)
                    {
                        int compteur_node = 0;
                        for(int m = 0; m < nodes_commun_tag.size(); m++)
                        {
                            if(nodes_commun_tag[m] == loop_element_2d_nodes_tag)
                            {
                                compteur_node++;
                            }
                        }
                        if(compteur_node == 0)
                        {
                            nodes_commun_tag.push_back(loop_element_2d_nodes_tag);
                        }
                        int compteur_element = 0;
                        for(int m = 0; m < elements_commun_tag.size(); m++)
                        {
                            if(elements_commun_tag[m] == elementVector[l].tag)
                            {
                                compteur_element++;
                            }
                        }
                        if(compteur_element == 0)
                        {
                            elements_commun_tag.push_back(elementVector[l].tag);
                            elements_commun_indice.push_back(l);
                        }
                    }
                }
            }
            if(show_element_node_data)
            {
                cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
                for(int m = 0; m < nodes_commun_tag.size(); m++)
                {
                    cout << nodes_commun_tag[m] << " ";
                }
                cout << endl;
                cout << "elements commun: " << elements_commun_tag.size() << ": ";
                for(int m = 0; m < elements_commun_tag.size(); m++)
                {
                    cout << elements_commun_tag[m] << " ";
                }
                cout << endl;
                cout << endl;
            }
            if(nodes_commun_tag.size() == 0)
            {
                for(int k = 0; k < nb_nodes_element2d; k++)
                {
                    int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                    data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                }
            }
            if(nodes_commun_tag.size() == 1)
            {
                if(elements_commun_tag.size() == 2)
                {
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                        if(loop_element_2d_nodes_tag == nodes_commun_tag[0])
                        {
                            double phi = 0;
                            phi = elementVector[elements_commun_indice[0]].bcValue[0];
                            phi = phi + elementVector[elements_commun_indice[1]].bcValue[0];
                            phi = phi/2.0;
                            data[i][j][k] = phi;
                        }
                        else if(data[i][j][k] == HUGE_VAL)
                        {
                            data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                        }
                    }
                }
                else
                {
                    cout << "problème dans la géométrie node = 1" << endl;
                    cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
                    for(int m = 0; m < nodes_commun_tag.size(); m++)
                    {
                        cout << nodes_commun_tag[m] << " ";
                    }
                    cout << endl;
                    cout << "elements commun: " << elements_commun_tag.size() << ": ";
                    for(int m = 0; m < elements_commun_tag.size(); m++)
                    {
                        cout << elements_commun_tag[m] << " ";
                    }
                    cout << endl;
                    cout << endl;
                }
            }
            if(nodes_commun_tag.size() == 2)
            {
                if(elements_commun_tag.size() == 3)
                {
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                        for(int l = 0; l < elements_commun_tag.size(); l++)
                        {
                            int compteur = 0;
                            for(int m = 0; m < nodes_commun_tag.size(); m++)
                            {
                                if(elementVector[elements_commun_indice[l]].nodeTags[0] == nodes_commun_tag[m]||elementVector[elements_commun_indice[l]].nodeTags[1] == nodes_commun_tag[m])
                                {
                                    compteur++;
                                }
                            }
                            if(compteur == 2)
                            {
                                if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]||loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                                {

                                    data[i][j][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                }
                            }
                        }
                        if(data[i][j][k] == HUGE_VAL)
                        {
                            data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                        }
                    }
                }
                else if(elements_commun_tag.size() == 4)
                {
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                        for(int l = 0; l < elements_commun_tag.size(); l++)
                        {
                            if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]||loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                            {
                                if(data[i][j][k] == HUGE_VAL)
                                {
                                    data[i][j][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                }
                                else
                                {
                                    data[i][j][k] = (data[i][j][k] + elementVector[elements_commun_indice[l]].bcValue[0])/2.0;
                                }
                            }
                        }
                        if(data[i][j][k] == HUGE_VAL)
                        {
                            data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                        }
                    }
                }
                else
                {
                    cout << "problème dans la géométrie node = 2" << endl;
                    cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
                    for(int m = 0; m < nodes_commun_tag.size(); m++)
                    {
                        cout << nodes_commun_tag[m] << " ";
                    }
                    cout << endl;
                    cout << "elements commun: " << elements_commun_tag.size() << ": ";
                    for(int m = 0; m < elements_commun_tag.size(); m++)
                    {
                        cout << elements_commun_tag[m] << " ";
                    }
                    cout << endl;
                    cout << endl;
                }
            }
            if(nodes_commun_tag.size() == 3)
            {
                if(elements_commun_tag.size() == 4)
                {
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                        for(int l = 0; l < elements_commun_tag.size(); l++)
                        {
                            int compteur = 0;
                            for(int m = 0; m < nodes_commun_tag.size(); m++)
                            {
                                if(elementVector[elements_commun_indice[l]].nodeTags[0] == nodes_commun_tag[m]||elementVector[elements_commun_indice[l]].nodeTags[1] == nodes_commun_tag[m])
                                {
                                    compteur++;
                                }
                            }
                            if(compteur == 2)
                            {
                                if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]||loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                                {
                                    if(data[i][j][k] == HUGE_VAL)
                                    {
                                        data[i][j][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                    }
                                    else
                                    {
                                        data[i][j][k] = (data[i][j][k] + elementVector[elements_commun_indice[l]].bcValue[0])/2.0;
                                    }
                                }
                            }
                        }
                        if(data[i][j][k] == HUGE_VAL)
                        {
                            data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                        }
                    }
                }
                else
                {
                    cout << "problème dans la géométrie node = 3" << endl;
                    cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
                    for(int m = 0; m < nodes_commun_tag.size(); m++)
                    {
                        cout << nodes_commun_tag[m] << " ";
                    }
                    cout << endl;
                    cout << "elements commun: " << elements_commun_tag.size() << ": ";
                    for(int m = 0; m < elements_commun_tag.size(); m++)
                    {
                        cout << elements_commun_tag[m] << " ";
                    }
                    cout << endl;
                    cout << endl;
                }
                for(int k = 0; k < nb_nodes_element2d; k++)
                {
                    int loop_element_2d_nodes_tag = nodeTags2d[i][nb_nodes_element2d*j+k];
                    if(data[i][j][k] == HUGE_VAL)
                    {
                        data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                    }
                }
            }
            elements_commun_tag.clear();
            elements_commun_indice.clear();
            nodes_commun_tag.clear();
        }
    }
}
