#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include "functionsBEM.hpp"
#include <omp.h>
#include <algorithm>
#include <list>

/*-----------------todo déjà fait mais il faudra peut être en refaire certains par après donc je les garde------------------*/
//todo: séparer en différentes fonctions
//todo: comparer au livre BEM
//todo: voir pour gérer des trous dans les surfaces
//todo: combiner avec la FEM
//todo: faire plein de tests différent (dans les .geo)
//todo: afficher et calculer la force electrostatique et champs electriques
//todo: on a besoin de la pression electrostatique à la frontière mais on peut essayer de l'afficher partout dans le domaine => representation formula pour le domain ou par elementnodedata et gradient
//todo: fusionner avec PY
//todo: potentiel2d fonction
//todo: gérer permuttation matrice plus efficacement (en vrai mnt, osef)
//todo: hybrid mesh + éléments 1d 3 noeuds
//todo: mettre en parallèle
//todo: elementtype pour optimiser elementnodedata
//todo: commenter le code

/*-----------------todo qu'il reste vraiment à faire------------------*/
//todo: convergence
//todo: afficher les déplacements, le potentiel en fonction de la distance dans l'élément (une coupe 1D dans le domaine)
//todo: representation formula pour le champ electrique
//todo: faire une fonction pour chrono afin de mesurer le temps de toutes les différentes parties du code
//todo: créer d'autres .geo avec des géométries plus complexes
//gérer des matériaux composites

using namespace gmsh;
using namespace std;
using namespace Eigen;

void element_node_data(vector<vector<vector<double>>> &data, vector<vector<size_t>> &elementTags2d_domain_BEM, vector<int> &elementTypes2d_domain_BEM, vector<vector<size_t>> &nodeTags2d_domain_BEM, vector<elementStruct> &elementVector, map<int, double> &nodeMapDomain, bool &show_element_node_data, int &entities_tags_for_physical_group_domain_BEM)
{
    //The first part of the algorithm consists in recovering the tag of the 2d elements of the domain which have at least one node in common with the boundary
    //We initialise the first dimension of the vector that will contain all 2d elements in common with the boundary with the number of element types present in the domain.
    vector<vector<size_t>> elementTags2d_boundary(elementTypes2d_domain_BEM.size());

    //this loop goes through all the nodes of the boundary (node 1 of each element). We then use getElementsByCoordinates to obtain all the elements in common with the node considered. We only take the 2d elements.
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        double posx = elementVector[i].x1;
        double posy = elementVector[i].y1;
        double posz = 0;
        vector<size_t> elementTags;
        model::mesh::getElementsByCoordinates(posx,posy,posz,elementTags,2);

        //we loop on the 2d elements which are in common with the considered node of the boundary. Depending on their type, they are added to the vector comprising all the 2d elements common with the boundary.
        for(size_t j = 0; j < elementTags.size(); j++)
        {
            int elementType;
            vector<size_t> nodeTags;
            int dim;
            int tag;
            vector<int> physicalTags;
            model::mesh::getElement(elementTags[j],elementType,nodeTags,dim,tag);
            model::getPhysicalGroupsForEntity(dim,tag,physicalTags);

            for(size_t k = 0; k < elementTypes2d_domain_BEM.size(); k++)
            {
                for(size_t l = 0; l < physicalTags.size(); l++)
                {
                    if(elementType == elementTypes2d_domain_BEM[k] && physicalTags[l] == entities_tags_for_physical_group_domain_BEM)
                    {
                        elementTags2d_boundary[k].push_back(elementTags[j]);
                    }
                }
            }
        }
    }

    //This loop allows duplicate tags to be removed. If an element has several nodes with the border, it has been counted several times before.
    for(size_t i = 0; i < elementTypes2d_domain_BEM.size(); i++)
    {
        sort( elementTags2d_boundary[i].begin(), elementTags2d_boundary[i].end() );
        elementTags2d_boundary[i].erase( unique( elementTags2d_boundary[i].begin(), elementTags2d_boundary[i].end() ), elementTags2d_boundary[i].end() );
    }

    //This map may be replaced in the future. It allows to retain the position of the elements of the elementTag2d loop and linked their position to the loop of the elementTag2d_boundary (the tags of the elements in common with the boundary). This map is the link between the array of tags of the 2d elements of the boundary and that of the tags of all the 2d elements.
    map<int, int> tagposition;

    //The first dimension of data corresponds to the number of element types.
    data.resize(elementTags2d_domain_BEM.size());
    for(size_t i = 0; i < elementTags2d_domain_BEM.size(); i++)
    {
        //The second dimension of data corresponds to the number of element tags for on element type.
        data[i].resize(elementTags2d_domain_BEM[i].size());
        string elementName;
        int dim;
        int order;
        int numNodes;
        vector<double> localNodeCoord;
        int numPrimaryNodes;
        model::mesh::getElementProperties(elementTypes2d_domain_BEM[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);

        //getElementProperties was used to get the number of nodes for a specific element
        int nb_nodes_element2d = numNodes;
        for(size_t j = 0; j < elementTags2d_domain_BEM[i].size(); j++)
        {
            //The third dimension of data corresponds to the number of nodes(third: k loop) for an element tags(second: j loop) for an element type(first: i loop).
            data[i][j].resize(nb_nodes_element2d,HUGE_VAL);
            tagposition[elementTags2d_domain_BEM[i][j]] = j;
        }
    }

    //The second part of the algorithm will consist of assigning the correct potential value to the nodes of the elements in common with the boundary.
    //The first loop is on the elementType in 2d in common with the boundary.
    for(size_t i = 0; i < elementTags2d_boundary.size(); i++)
    {
        string elementName;
        int dim;
        int order;
        int numNodes;
        vector<double> localNodeCoord;
        int numPrimaryNodes;
        model::mesh::getElementProperties(elementTypes2d_domain_BEM[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
        //getElementProperties was used to get the number of nodes for a specific element
        int nb_nodes_element2d = numNodes;

        //The second loop is on the elementTag in 2d in common with the boundary.
        for(size_t j = 0; j < elementTags2d_boundary[i].size(); j++)
        {
            vector<int> elements_commun_tag(0,0);
            vector<int> elements_commun_indice(0,0);
            vector<int> nodes_commun_tag(0,0);
            int elementType;
            vector<size_t> nodeTags_element2d;
            int dim;
            int tag;
            model::mesh::getElement(elementTags2d_boundary[i][j],elementType,nodeTags_element2d,dim,tag);
            //The third loop is on the nodes of a specific elementTag in 2d in common with the boundary.
            for(size_t k = 0; k < nb_nodes_element2d; k++)
            {
                int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                //For each node of the 2d elements, we loop over all the 1d elements of the boundary and look for both the node in common and the elements in common. Indeed, on the boundary, a specific node belongs to two 1d elements.
                for(size_t l = 0; l < elementVector.size(); l++)
                {
                    //condition vérifié si le noeud de l’élément 2d est égal au noeud d’un élément 1d de la frontière.
                    if(elementVector[l].nodeTags[0] == loop_element_2d_nodes_tag||elementVector[l].nodeTags[1] == loop_element_2d_nodes_tag)
                    {
                        //this loop allows the node to be added only once if it is present several times. As the 1d elements of the boundary are traversed, 1 node encounters 2 elements and can therefore be counted multiple times.
                        int compteur_node = 0;
                        for(size_t m = 0; m < nodes_commun_tag.size(); m++)
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
                        //The same procedure is applied for the elements. You only want to count them once if they are in common
                        int compteur_element = 0;
                        for(size_t m = 0; m < elements_commun_tag.size(); m++)
                        {
                            if(elements_commun_tag[m] == elementVector[l].tag)
                            {
                                compteur_element++;
                            }
                        }
                        if(compteur_element == 0)
                        {
                            //Both the tag of the element and its position in the loop on the elements of the boundary are retained.
                            elements_commun_tag.push_back(elementVector[l].tag);
                            elements_commun_indice.push_back(l);
                        }
                    }
                }
            }
            //For each 2d element in common with the boundary, displays the tag of the nodes and the tag of the common 1d elements of the boundary.
            if(show_element_node_data)
            {
                cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
                for(size_t m = 0; m < nodes_commun_tag.size(); m++)
                {
                    cout << nodes_commun_tag[m] << " ";
                }
                cout << endl;
                cout << "elements commun: " << elements_commun_tag.size() << ": ";
                for(size_t m = 0; m < elements_commun_tag.size(); m++)
                {
                    cout << elements_commun_tag[m] << " ";
                }
                cout << endl;
                cout << endl;
            }
            //Several geometric configurations are to be considered. Here, the 2d element has 1 node in common with the boundary and must therefore touch 2 of its elements.
            if(nodes_commun_tag.size() == 1)
            {
                if(elements_commun_tag.size() == 2)
                {
                    //We loop over the nodes of the 2d element, find the one that is common to the boundary and assign it the average value of the potential of the 2 elements it touches. Indeed the potential is discontinuous from one element of the boundary to the other.
                    for(size_t k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        if(loop_element_2d_nodes_tag == nodes_commun_tag[0])
                        {
                            double phi = 0;
                            phi = elementVector[elements_commun_indice[0]].bcValue[0];
                            phi = phi + elementVector[elements_commun_indice[1]].bcValue[0];
                            phi = phi/2.0;
                            data[i][tagposition[elementTags2d_boundary[i][j]]][k] = phi;
                        }
                    }
                }
                //displays an error in the case of geometry where 1 node is in contact with more than 2 elements.
                else
                {
                    disp_element_node_data(nodes_commun_tag, elements_commun_tag);
                }
            }
            //The principle is the same as before. For the case where there are 2 nodes in common, the 2d element can touch either 3 or 4 1d elements of the boundary. In the first case, it has one of its sides in common with the boundary, the other case can be found in the corners of the domain (two of its corners touch the boundary in 1 node and meet 2 elements like case 1).
            if(nodes_commun_tag.size() == 2)
            {
                if(elements_commun_tag.size() == 3)
                {
                    for(size_t k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        for(size_t l = 0; l < elements_commun_tag.size(); l++)
                        {
                            //Among the 3 elements of the boundary in common, we look for the one which has 2 nodes in common with the 2d element. We look for the side and therefore the 1d element that they have in common.
                            int compteur_nodes = 0;
                            for(size_t m = 0; m < nodes_commun_tag.size(); m++)
                            {
                                if(elementVector[elements_commun_indice[l]].nodeTags[0] == nodes_commun_tag[m]||elementVector[elements_commun_indice[l]].nodeTags[1] == nodes_commun_tag[m])
                                {
                                    compteur_nodes++;
                                }
                            }
                            if(compteur_nodes == 2)
                            {
                                //Once we have found the element in common, the node of the 2d element on which we are looping must be the right one
                                if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]||loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                                {
                                    //The constant value of the boundary element in common is assigned to the 2 nodes of the 2d element that touches the boundary.
                                    data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                }
                                //In the case of second order elements, the value of the common 1d boundary element is assigned to the middle node.
                                if(order == 2)
                                {
                                    if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].midNodeTags)
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                    }
                                }
                            }
                        }
                    }
                }
                else if(elements_commun_tag.size() == 4)
                {
                    for(size_t k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        for(size_t l = 0; l < elements_commun_tag.size(); l++)
                        {
                            //This geometry behaves almost like the first case (1 node, 2 elements) except that it is true for 2 different nodes of the 2d element. We loop over the elements of the boundary. If one of its nodes is common with the boundary, it is assigned the value of the potential of this first element. When the second one comes, it is averaged with the potential already present on the node.
                            if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]||loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                            {
                                if(data[i][tagposition[elementTags2d_boundary[i][j]]][k] == HUGE_VAL)
                                {
                                    data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                }
                                else
                                {
                                    data[i][tagposition[elementTags2d_boundary[i][j]]][k] = (data[i][tagposition[elementTags2d_boundary[i][j]]][k] + elementVector[elements_commun_indice[l]].bcValue[0])/2.0;
                                }
                            }
                        }
                    }
                }
                else
                {
                    disp_element_node_data(nodes_commun_tag, elements_commun_tag);
                }
            }
            //The last case is that of 3 common nodes which therefore affect 4 elements. This case is mainly found in the corners. The average potential is attributed to the node that touches 2 elements (in the corner of a square for example). The other 2 nodes must simply match both the side in common and the node in common. The principle of assigning the potential value of the node is the superposition of the previous methods.
            if(nodes_commun_tag.size() == 3)
            {
                if(elements_commun_tag.size() == 4)
                {
                    for(size_t k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        for(size_t l = 0; l < elements_commun_tag.size(); l++)
                        {
                            int compteur_nodes = 0;
                            for(size_t m = 0; m < nodes_commun_tag.size(); m++)
                            {
                                if(elementVector[elements_commun_indice[l]].nodeTags[0] == nodes_commun_tag[m]||elementVector[elements_commun_indice[l]].nodeTags[1] == nodes_commun_tag[m])
                                {
                                    compteur_nodes++;
                                }
                            }
                            if(compteur_nodes == 2)
                            {
                                if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]||loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                                {
                                    if(data[i][tagposition[elementTags2d_boundary[i][j]]][k] == HUGE_VAL)
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                    }
                                    else
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = (data[i][tagposition[elementTags2d_boundary[i][j]]][k] + elementVector[elements_commun_indice[l]].bcValue[0])/2.0;
                                    }
                                }
                                //In the case of second order elements, the value of the common 1d boundary element is assigned to the middle node.
                                if(order == 2)
                                {
                                    if(loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].midNodeTags)
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    disp_element_node_data(nodes_commun_tag, elements_commun_tag);
                }
            }
            //We clear the vectors that we push at each loop
            elements_commun_tag.clear();
            elements_commun_indice.clear();
            nodes_commun_tag.clear();
        }
    }


    //The last step is to loop over all the 2d elements in the domain. If the HUGE_VAL value is still present, this means that the element node does not belong to the boundary (since the value of the potential has not been modified) and is therefore part of the domain. We then use the map which associates the value of the potential with the tag of the node inside the domain.
    for(size_t i = 0; i < elementTags2d_domain_BEM.size(); i++)
    {
        string elementName;
        int dim;
        int order;
        int numNodes;
        vector<double> localNodeCoord;
        int numPrimaryNodes;
        model::mesh::getElementProperties(elementTypes2d_domain_BEM[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
        int nb_nodes_element2d = numNodes;
        for(size_t j = 0; j < elementTags2d_domain_BEM[i].size(); j++)
        {
            int elementType;
            vector<size_t> nodeTags_element2d;
            int dim;
            int tag;
            model::mesh::getElement(elementTags2d_domain_BEM[i][j],elementType,nodeTags_element2d,dim,tag);
            for(size_t k = 0; k < nb_nodes_element2d; k++)
            {
                int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                if(data[i][j][k] == HUGE_VAL)
                {
                    data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                }
            }
        }
    }
}

void dispNodeCoordFun(const vector<size_t> &allNodeTags, const vector<double> &allCoord)
{
    for(size_t i = 0; i < allNodeTags.size(); i++)
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
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
    }
}

void getBC(map<size_t, int> allNodeIndex, const vector<double> &allCoord, vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy, vector<size_t> &NodeTags2d_domain, vector<double> &NodeCoord2d_domain,            vector<int> &elementTypes2d_domain_BEM, vector<vector<size_t>> &elementTags2d_domain_BEM, vector<vector<size_t>> &elnodeTags2d_domain_BEM, int &entities_tags_for_physical_group_domain_BEM)
{
    // Gets the dim and the tags of the physical groups
    vectorpair physicalGroupsDimTags;
    model::getPhysicalGroups(physicalGroupsDimTags);

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

        // focus only on dirichlet or neumann boundary conditions -> associated to BEM
        if(words.size() == 3 && (words[2].compare("neumann") == 0 || words[2].compare("dirichlet") == 0))
        {
            bcMap[words[1]] = {words[2], values[0]};
        }
    }

    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    for (int i = 0; i < physicalGroupsDimTags.size(); i++)
    {
        int physicalGroupsDim = physicalGroupsDimTags[i].first;
        int physicalGroupsTag = physicalGroupsDimTags[i].second;
        string name;
        model::getPhysicalName(physicalGroupsDim, physicalGroupsTag, name);

        if(name == "BEM_domain")
        {
            vector<size_t> nodeTags;
            vector<double> coord;
            model::mesh::getNodesForPhysicalGroup(physicalGroupsDim,physicalGroupsTag,nodeTags,coord);

            NodeTags2d_domain.insert(NodeTags2d_domain.end(), nodeTags.begin(), nodeTags.end());
            NodeCoord2d_domain.insert(NodeCoord2d_domain.end(), coord.begin(), coord.end());

            // c'est un peu difficile d'expliquer pourquoi je souvegarde le tag du domaine BEM. C'est pour un algo qui doit récupérer tous les éléments 2d de la frontière mais doit faire la distinction entre les éléments du domaine FEM et BEM
            entities_tags_for_physical_group_domain_BEM = physicalGroupsTag;
            vector<int> entities_tags_2d_domain_BEM;
            model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entities_tags_2d_domain_BEM);

            for(size_t j=0; j< entities_tags_2d_domain_BEM.size(); j++)
            {
                vector<int> elementTypes;
                vector<vector<size_t> > elementTags;
                vector<vector<size_t> > nodeTags;
                model::mesh::getElements(elementTypes2d_domain_BEM, elementTags2d_domain_BEM, elnodeTags2d_domain_BEM, physicalGroupsDim, entities_tags_2d_domain_BEM[j]);
                elementTypes2d_domain_BEM.insert(elementTypes2d_domain_BEM.end(), elementTypes.begin(), elementTypes.end());
                elementTags2d_domain_BEM.insert(elementTags2d_domain_BEM.end(), elementTags.begin(), elementTags.end());
            }
        }
        //Boundaries that separate 2 types of elements within the domain are assumed to have the name = "between". The nodes and their coordinates are added to the vectors containing all the tags and coordinates of the nodes in the domain.
        else if(name == "between")
        {
            vector<size_t> nodeTags;
            vector<double> coord;
            model::mesh::getNodesForPhysicalGroup(physicalGroupsDim,physicalGroupsTag,nodeTags,coord);

            NodeTags2d_domain.insert(NodeTags2d_domain.end(), nodeTags.begin(), nodeTags.end());
            NodeCoord2d_domain.insert(NodeCoord2d_domain.end(), coord.begin(), coord.end());
        }

        if(bcMap.find(name) == bcMap.end())
        {
            continue;
        }

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

            for(size_t el = 0; el < elementTags1.size(); el++)
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
            for(size_t el = 0; el < elementTags8.size(); el++)
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


void computeGeometricParam(elementStruct &element, map<size_t, int> allNodeIndex, const vector<double> &allCoord)
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


double computeH(const double &x, const double &y, const elementStruct &element)
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


double computeG(const double &x, const double &y, const elementStruct &element, const string &integrationType)
{
    // Computation of the integration points (parametric) and weights.
    vector<double> localCoord;
    vector<double> weights;
    model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Computation of the sum (of integrand evaluated at integration point * weight)
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
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
    double x = 0;
    double y = 0;
    double phi = 0;
    double h = 0;
    double g = 0;

#pragma omp parallel for schedule(static) private(x,y,phi,h,g) shared(NodeTags2d_domain,allCoord,integrationType,elementVector,nodeMapDomain)
    for(size_t i = 0; i < NodeTags2d_domain.size(); i++)
    {
        // Coordinates of the node
        x = allCoord[3*(NodeTags2d_domain[i]-1)];
        y = allCoord[3*(NodeTags2d_domain[i]-1) + 1];

        phi = 0;
        for(size_t j = 0; j < elementVector.size(); j++)
        {
            h = computeH(x, y, elementVector[j]);
            g = computeG(x, y, elementVector[j], integrationType);

            phi = phi + h*elementVector[j].bcValue[0] - g*elementVector[j].bcValue[1];
        }
        // Potential at the node is phi.
#pragma omp critical
        nodeMapDomain[NodeTags2d_domain[i]] = phi;
    }
}

void computeData(const vector<elementStruct> &elementVector, const vector<int> &elementTypes, const vector<vector<size_t>> &elementTags, vector<size_t> &tags, vector<vector<double>> &data, const string integrationType)
{
        double x = 0;
        double y = 0;
        double h = 0;
        double g = 0;
        double phi = 0;
        int elType = 0;
        vector<double> barycenters;
        size_t elementCount = 0;

        for(size_t i = 0; i < elementTypes.size(); i++)
        {
            elType = elementTypes[i];
            model::mesh::getBarycenters(elType, -1, false, false, barycenters);
#pragma omp parallel for schedule(guided) private(x,y,phi,h,g) shared(elementTags,elementVector,tags,data,elementCount,barycenters)
            for(size_t j = 0; j < elementTags[i].size(); j++)
            {
                x = barycenters[3*j];
                y = barycenters[3*j+1];

                phi = 0;
                for(size_t k = 0; k < elementVector.size(); k++)
                {
                    h = computeH(x, y, elementVector[k]);
                    g = computeG(x, y, elementVector[k], integrationType);
                    phi += h*elementVector[k].bcValue[0] - g*elementVector[k].bcValue[1];
                }
                // Potential at the node is phi.
    #pragma omp critical
                {
                    tags[elementCount + j] = elementTags[i][j];
                    data[elementCount + j].resize(1,phi);
                }
            }
            elementCount = elementCount + elementTags[i].size();
        }

    /*int elementCount = 0;
     for(size_t i = 0; i < elementTypes.size(); i++)
     {
     int elType = elementTypes[i];
     vector<double> barycenters;
     model::mesh::getBarycenters(elType, -1, false, false, barycenters);

     //#pragma omp parallel for
     for(size_t j = 0; j < elementTags[i].size(); j++)
     {
     double x = barycenters[3*j];
     double y = barycenters[3*j + 1];

     double phi = 0;
     for(size_t k = 0; k < elementVector.size(); k++)
     {
     double h = computeH(x, y, elementVector[k]);
     double g = computeG(x, y, elementVector[k], "Gauss4");
     phi += h*elementVector[k].bcValue[0] - g*elementVector[k].bcValue[1];
     }
     // Potential at the node is phi.
     tags[elementCount] = elementTags[i][j];
     data[elementCount] = vector<double> (1);
     data[elementCount][0] = phi;
     elementCount++; // cette ligne empêche de paralléliser
     }
     }
     cout << elementCount;*/

    /*--------ajout Louis-------*/
//    int max_threads = omp_get_max_threads();
//    omp_set_num_threads(max_threads);
//    //each thread will have its own data vectors (lists to simplify the push_back operation), then the vectors will be merged
//    vector<int> elementCountThreads(max_threads); // elementCountPerThread
//    vector<list<size_t>> tagsThreads(max_threads);
//    vector<list<vector<double>>> dataThreads(max_threads);
//    for(int i = 0; i < max_threads; i++)
//    {
//        elementCountThreads[i] = 0;
//    }
//    for(size_t i = 0; i < elementTypes.size(); i++)
//    {
//        int elType = elementTypes[i];
//        vector<double> barycenters;
//        model::mesh::getBarycenters(elType, -1, false, false, barycenters);
//
//#pragma omp parallel for
//        for(size_t j = 0; j < elementTags[i].size(); j++)
//        {
//            double x = barycenters[3*j];
//            double y = barycenters[3*j + 1];
//
//            double phi = 0;
//            for(size_t k = 0; k < elementVector.size(); k++)
//            {
//                double h = computeH(x, y, elementVector[k]);
//                double g = computeG(x, y, elementVector[k], "Gauss4");
//                phi += h*elementVector[k].bcValue[0] - g*elementVector[k].bcValue[1];
//            }
//            // Potential at the node is phi.
//            int my_thread = omp_get_thread_num();
//            tagsThreads[my_thread].push_back(elementTags[i][j]);
//            vector<double> v(1);
//            v[0] = phi;
//            dataThreads[my_thread].push_back(v);
//            elementCountThreads[my_thread] = elementCountThreads[my_thread] + 1;
//        }
//    }
//    // merging the results
//    for(int i = 0; i < max_threads; i++)
//    {
//        vector<size_t> myTags(tagsThreads[i].begin(), tagsThreads[i].end());
//        vector<vector<double>> myData(dataThreads[i].begin(), dataThreads[i].end());
//
//        tags.insert(tags.end(), myTags.begin(), myTags.end());
//        data.insert(data.end(), myData.begin(), myData.end());
//    }
}

void disp_element_node_data(vector<int> &nodes_commun_tag, vector<int> &elements_commun_tag)
{
    cout << "problème dans la géométrie, noeuds communs = " << nodes_commun_tag.size() << endl;
    cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
    for(size_t m = 0; m < nodes_commun_tag.size(); m++)
    {
        cout << nodes_commun_tag[m] << " ";
    }
    cout << endl;
    cout << "elements commun: " << elements_commun_tag.size() << ": ";
    for(size_t m = 0; m < elements_commun_tag.size(); m++)
    {
        cout << elements_commun_tag[m] << " ";
    }
    cout << endl;
    cout << endl;
}

void display_time(chrono::time_point<chrono::system_clock> &start, chrono::time_point<chrono::system_clock> &end, string &function_name, bool &show_time)
{
    //displays the time of the desired code section
    if(show_time)
    {
        cout << endl;
        cout << function_name << endl;
        cout << "Elapsed time in nanoseconds: "
        << chrono::duration_cast<chrono::nanoseconds>(end - start).count()
        << " ns" << endl;

        cout << "Elapsed time in microseconds: "
        << chrono::duration_cast<chrono::microseconds>(end - start).count()
        << " µs" << endl;

        cout << "Elapsed time in milliseconds: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << endl;

        cout << "Elapsed time in seconds: "
        << chrono::duration_cast<chrono::seconds>(end - start).count()
        << " sec";
        cout << endl;
        cout << endl;
    }
}
