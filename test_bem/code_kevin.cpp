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
//todo: combiner avec la FEM
//todo: faire plein de tests différent (dans les .geo)
//todo: afficher et calculer la force electrostatique et champs electriques
//todo: on a besoin de la pression electrostatique à la frontière mais on peut essayer de l'afficher partout dans le domaine => representation formula pour le domain ou par elementnodedata
//todo: fusionner avec PY
//todo: potentiel2d fonction
//todo: gérer permuttation matrice plus efficacement (en vrai mnt, osef)
//todo: hybrid mesh + éléments 1d 3 noeuds

//todo: convergence
//todo: mettre en parallèle
//todo: afficher les déplacements, le potentiel en fonction de la distance dans l'élément (une coupe 1D dans le domaine)
//todo: elementtype pour optimiser elementnodedata
//todo: commenter le code

using namespace gmsh;
using namespace std;
using namespace Eigen;

int solverBEMFun(int argc, char **argv, vector<size_t> &finalElementTags, vector<double> &electrostaticPressure);
void dispNodeCoordFun(const vector<size_t> &allNodeTags, const vector<double> &allCoord);
void dispMatricesFun(const MatrixXd &A, const MatrixXd &B, const MatrixXd &c, const MatrixXd &b, const MatrixXd &x);
void dispValuesFun(const vector<elementStruct> &elementVector);
void getBC(map<std::size_t, int> allNodeIndex, const vector<double> &allCoord, vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy, vector<size_t> &NodeTags2d_domain, vector<double> &NodeCoord2d_domain);
void computeGeometricParam(elementStruct &element, map<std::size_t, int> allNodeIndex, const vector<double> &allCoord);
double computeH(const double x, const double y, const elementStruct &element);
double computeG(const double x, const double y, const elementStruct &element, const string integrationType);
void computeData(const vector<elementStruct> &elementVector, const vector<int> &elementTypes, const vector<vector<size_t>> &elementTags, vector<size_t> &tags, vector<vector<double>> &data, const string integrationType);
void displayData(const string viewName, const string modelName, const string fileName, const vector<size_t> &tags, const vector<vector<double>> &data);
void compute_internal_phi(vector<double> &allCoord,vector<size_t> &NodeTags2d_domain, vector<elementStruct> &elementVector, map<int,double> &nodeMapDomain, const string integrationType);
void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d, vector<int> &elementTypes2d, vector<vector<size_t>> &nodeTags2d, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data);
void disp_element_node_data(vector<int> &nodes_commun_tag, vector<int> &elements_commun_tag);

int solverBEMFun(int argc, char **argv, vector<size_t> &finalElementTags, vector<double> &electrostaticPressure)
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
    bool openGUI = 1;        // Graphical interface of GMSH.
    bool show_element_node_data = 0;
    int visualisation = 1;  //0: no visualisation, 1: element node data, 2: element data

    // Initialises a vector of elements (that will contain all the elements of interrest)
    vector<elementStruct> elementVector;
    string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.

    // Gets the tags, the coordinates and the parametric coordinates of the nodes.
    map<size_t, int> allNodeIndex; // nodeTag -> position in allNodeTags and allCoord
    vector<size_t> allNodeTags;
    vector<double> allCoord;
    vector<double> allParametricCoord;
    model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);

    for(int i = 0; i < allNodeTags.size(); i++)
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

    vector<size_t> NodeTags2d_domain;
    vector<double> NodeCoord2d_domain;
    vector<double> NodeParametricCoord2d_domain;
    model::mesh::getNodes(NodeTags2d_domain, NodeCoord2d_domain, NodeParametricCoord2d_domain,2);

    getBC(allNodeIndex,allCoord, elementVector, dirichletCount, dispHierarchy, NodeTags2d_domain, NodeCoord2d_domain);

    //pour optimiser element node data
//    double posx = elementVector[0].x1;
//    double posy = elementVector[0].y1;
//    double posz = 0;
//    vector<size_t> elementTags;
//    model::mesh::getElementsByCoordinates(posx,posy,posz,elementTags);

    vector<string> names;
    model::list(names);
    cout << "the file contains " << names.size() << " model names\n";

    vector<int> elementTypes2d;
    vector<vector<size_t> > elementTags2d;
    vector<vector<size_t> > nodeTags2d;
    model::mesh::getElements(elementTypes2d,elementTags2d,nodeTags2d,2);

    // Declaration of matrices A, B and vector c (the matrix type is perhaps not the most efficient).
    MatrixXd A = MatrixXd::Zero(elementVector.size(), elementVector.size());
    MatrixXd B = MatrixXd::Zero(elementVector.size(), elementVector.size());
    MatrixXd c = MatrixXd::Zero(elementVector.size(), 1);

    //  Double loop over each element in order to compute the element h and g that will be given to the matrices A and B.
    for(int i = 0; i < elementVector.size(); i++)
    {
        for(int j = 0; j < elementVector.size(); j++)
        {
            // Will be computed and stored in A and B.
            double h;
            double g;

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

    MatrixXd b = B*c;
    MatrixXd x = A.fullPivLu().solve(b);

    // Matrices verification
    if(dispMatrices)
    {
        dispMatricesFun(A, B, c, b, x);
    }

    // Assign each value of x to the right element.
    for(int i = 0; i < elementVector.size(); i++)
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

    if(visualisation == 1)
    {
        vector<vector<vector<double>>> data;
        map<int, double> nodeMapDomain;

        compute_internal_phi(allCoord,NodeTags2d_domain,elementVector,nodeMapDomain,integrationType);

        element_node_data(data,elementTags2d,elementTypes2d,nodeTags2d,elementVector, nodeMapDomain,show_element_node_data);

        int viewNodeTag = view::add("ElementNodeData");

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
    }
    if(visualisation == 2)
    {
        vector<string> names;
        model::list(names);

        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);

        int nb2dElements = 0;
        for(int i = 0; i < elementTypes.size(); i++)
            nb2dElements += elementTags[i].size();

        vector<vector<double>> phi(nb2dElements);
        vector<size_t> tags(nb2dElements);

        computeData(elementVector, elementTypes, elementTags, tags, phi, integrationType);
        displayData("phi", names[0], "results.msh", tags, phi);

        if(openGUI)
        {
            fltk::run();
        }
    }
    double epsilon = 8.8541878128e-12;
    for(int i = 0; i < elementVector.size(); i++)
    {
        finalElementTags[i] = elementVector[i].tag;
        electrostaticPressure[i] = epsilon * pow(elementVector[i].bcValue[1], 2.0);
    }

    return 0;
}

void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d, vector<int> &elementTypes2d, vector<vector<size_t>> &nodeTags2d, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data)
{
    data.resize(elementTags2d.size());
    for(int i = 0; i < elementTags2d.size(); i++)
    {
        data[i].resize(elementTags2d[i].size());
        string elementName;
        int dim;
        int order;
        int numNodes;
        vector<double> localNodeCoord;
        int numPrimaryNodes;
        model::mesh::getElementProperties(elementTypes2d[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
        int nb_nodes_element2d = numNodes;
        for(int j = 0; j < elementTags2d[i].size(); j++)
        {
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
                    disp_element_node_data(nodes_commun_tag, elements_commun_tag);
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
                                if(elementTypes2d[i] == 9 || elementTypes2d[i] == 10)
                                {
                                    if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].midNodeTags)
                                    {
                                        data[i][j][k] = elementVector[elements_commun_indice[l]].bcValue[0];
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
                    disp_element_node_data(nodes_commun_tag, elements_commun_tag);
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
                                if(elementTypes2d[i] == 9 || elementTypes2d[i] == 10)
                                {
                                    if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].midNodeTags)
                                    {
                                        data[i][j][k] = elementVector[elements_commun_indice[l]].bcValue[0];
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
                    disp_element_node_data(nodes_commun_tag, elements_commun_tag);
                }
            }
            elements_commun_tag.clear();
            elements_commun_indice.clear();
            nodes_commun_tag.clear();
        }
    }
}

void dispNodeCoordFun(const vector<size_t> &allNodeTags, const vector<double> &allCoord)
{
    for(int i = 0; i < allNodeTags.size(); i++)
    {
        cout << "Node '" << allNodeTags[i] << "': (" << allCoord[3*i] << "; " << allCoord[3*i + 1] << "; " << allCoord[3*i + 2] << ")\n";
    }
    cout << endl;
}
void dispMatricesFun(const MatrixXd &A, const MatrixXd &B, const MatrixXd &c, const MatrixXd &b, const MatrixXd &x)
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
void dispValuesFun(const vector<elementStruct> &elementVector)
{
    for(int i = 0; i < elementVector.size(); i++)
    {
        cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
    }
}

void getBC(map<std::size_t, int> allNodeIndex, const vector<double> &allCoord, vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy, vector<size_t> &NodeTags2d_domain, vector<double> &NodeCoord2d_domain)
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

    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    for (int i = 0; i < physicalGroupsDimTags.size(); i++)
    {
        int physicalGroupsDim = physicalGroupsDimTags[i].first;
        int physicalGroupsTag = physicalGroupsDimTags[i].second;
        string name;
        model::getPhysicalName(physicalGroupsDim, physicalGroupsTag, name);

        if(name == "between")
        {
            vector<size_t> nodeTags;
            vector<double> coord;
            model::mesh::getNodesForPhysicalGroup(physicalGroupsDim,physicalGroupsTag,nodeTags,coord);

            NodeTags2d_domain.insert(NodeTags2d_domain.end(), nodeTags.begin(), nodeTags.end());
            NodeCoord2d_domain.insert(NodeCoord2d_domain.end(), coord.begin(), coord.end());
        }

        if(bcMap.find(name) == bcMap.end())
            continue;

        // Retrieves the type and the value of the boundary condition of this physical group.
        // /!\ Should check that this physical groups was well defined in the list of "Boundary Condition" !
        string bcType = bcMap[name].first;
        double bcValue = bcMap[name].second;

        // Prints the physical groups name, dim, tag, bcType and bcValue (if dispHierarchy != 0).
        if(dispHierarchy)
            cout << "Physical group '" << name << "' (" << physicalGroupsDim << "D) with tag '" << physicalGroupsTag << "' has a '" << bcType << "' type and the value " << bcValue << "\n";

        vector<int> entitiesTags;
        model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entitiesTags);

        for (int k = 0; k < entitiesTags.size(); k++)
        {
            // Print the tag of the entities of each physical group (if dispHierarchy != 0).
            if(dispHierarchy)
                cout << "\t contains entity '" << entitiesTags[k] << "' ";


            vector<size_t> elementTags1;
            vector<size_t> nodeTags1;
            vector<size_t> elementTags8;
            vector<size_t> nodeTags8;

            model::mesh::getElementsByType(1, elementTags1, nodeTags1, entitiesTags[k]);
            model::mesh::getElementsByType(8, elementTags8, nodeTags8, entitiesTags[k]);

            // Print the number of element of the given entity (if dispHierarchy != 0).
            if(dispHierarchy)
                cout << "which has " << elementTags1.size() + elementTags8.size() << " elements:\n";

            for(int el = 0; el < elementTags1.size(); el++)
            {
                int elTag = elementTags1[el];

                elementStruct element;
                element.tag = elTag;
                element.nodeTags[0] = nodeTags1[2*el];
                element.nodeTags[1] = nodeTags1[2*el + 1];

                computeGeometricParam(element, allNodeIndex, allCoord);
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
                    cout << "\t\t-Element '" << elTag << "' (" << element.bcType << "):\n";
                    cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[0]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[0]] + 1] << ")\n";
                    cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[1]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[1]] + 1] << ")\n";
                }
            }
            for(int el = 0; el < elementTags8.size(); el++)
            {
                int elTag = elementTags8[el];

                elementStruct element;
                element.tag = elTag;
                element.nodeTags[0] = nodeTags8[3*el];
                element.nodeTags[1] = nodeTags8[3*el + 1];
                element.midNodeTags = nodeTags8[3*el + 2];
                computeGeometricParam(element, allNodeIndex, allCoord);
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
                    cout << "\t\t-Element '" << elTag << "' (" << element.bcType << "):\n";
                    cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[0]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[0]] + 1] << ")\n";
                    cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[1]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[1]] + 1] << ")\n";
                }
            }
            if(dispHierarchy)
                cout << "\n";
        }
    }

}


void computeGeometricParam(elementStruct &element, map<std::size_t, int> allNodeIndex, const vector<double> &allCoord)
{
    int i1Tag = element.nodeTags[0];
    int i2Tag = element.nodeTags[1];

    int i1Index = allNodeIndex[i1Tag];
    int i2Index = allNodeIndex[i2Tag];

    // Get the coordinates of the nodes from their tag.
    element.x1 = allCoord[3*i1Index];
    element.y1 = allCoord[3*i1Index + 1];
    element.x2 = allCoord[3*i2Index];
    element.y2 = allCoord[3*i2Index + 1];

    element.midX = (element.x1 + element.x2) / 2;
    element.midY = (element.y1 + element.y2) / 2;
    element.deltaX = element.x2 - element.x1;
    element.deltaY = element.y2 - element.y1;
    element.length = sqrt(pow(element.deltaX, 2.0) + pow(element.deltaY, 2.0));
}


double computeH(const double x, const double y, const elementStruct &element)
{
    double deltaX1 = element.x1 - x;
    double deltaX2 = element.x2 - x;
    double deltaY1 = element.y1 - y;
    double deltaY2 = element.y2 - y;
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

    return alpha / (2*M_PI);
}


double computeG(const double x, const double y, const elementStruct &element, const string integrationType)
{
    // Computation of the integration points (parametric) and weights.
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
        double gX = element.midX + element.deltaX*xi/2;
        double gY = element.midY + element.deltaY*xi/2;

        // Distance between mid point of element 'i' and integration point.
        double r = sqrt(pow(x - gX, 2.0) + pow(y - gY, 2.0));
        gSum += log(r) * weights[k];
    }

    return gSum * element.length / (4*M_PI);
}

void displayData(const string viewName, const string modelName, const string fileName, const vector<size_t> &tags, const vector<vector<double>> &data)
{
    int viewTag = view::add(viewName);
    view::addModelData(viewTag, 0, modelName, "ElementData", tags, data);
    view::write(viewTag, fileName);

    option::setNumber("View[0].IntervalsType", 3);
    // option::setNumber("View[0].AdaptVisualizationGrid", 1); // Problem with this command !
    option::setNumber("View[0].MaxRecursionLevel", 3);
    option::setNumber("View[0].TargetError", -0.0001);
}

void compute_internal_phi(vector<double> &allCoord,vector<size_t> &NodeTags2d_domain, vector<elementStruct> &elementVector, map<int,double> &nodeMapDomain, const string integrationType)
{
    for(int i = 0; i < NodeTags2d_domain.size(); i++)
    {
        // Coordinates of the node
        double x = allCoord[3*(NodeTags2d_domain[i]-1)];
        double y = allCoord[3*(NodeTags2d_domain[i]-1) + 1];

        double phi = 0;
        for(int j = 0; j < elementVector.size(); j++)
        {
            double h = computeH(x, y, elementVector[j]);
            double g = computeG(x, y, elementVector[j], integrationType);

            phi = phi + h*elementVector[j].bcValue[0] - g*elementVector[j].bcValue[1];
        }
        // Potential at the node is phi.
        nodeMapDomain[NodeTags2d_domain[i]] = phi;
    }
}

void computeData(const vector<elementStruct> &elementVector, const vector<int> &elementTypes, const vector<vector<size_t>> &elementTags, vector<size_t> &tags, vector<vector<double>> &data, const string integrationType)
{
    int elementCount = 0;
    for(int i = 0; i < elementTypes.size(); i++)
    {
        int elType = elementTypes[i];
        vector<double> barycenters;
        model::mesh::getBarycenters(elType, -1, false, false, barycenters);

        for(int j = 0; j < elementTags[i].size(); j++)
        {
            double x = barycenters[3*j];
            double y = barycenters[3*j + 1];

            double phi = 0;
            for(int k = 0; k < elementVector.size(); k++)
            {
                double h = computeH(x, y, elementVector[k]);
                double g = computeG(x, y, elementVector[k], integrationType);
                phi += h*elementVector[k].bcValue[0] - g*elementVector[k].bcValue[1];
            }
            // Potential at the node is phi.
            tags[elementCount] = elementTags[i][j];
            data[elementCount] = vector<double> (1);
            data[elementCount][0] = phi;
            elementCount++;
        }
    }
}

void disp_element_node_data(vector<int> &nodes_commun_tag, vector<int> &elements_commun_tag)
{
    cout << "problème dans la géométrie, noeuds communs = " << nodes_commun_tag.size() << endl;
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

