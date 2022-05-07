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
//todo: faire une fonction pour chrono afin de mesurer le temps de toutes les différentes parties du code
//todo: convergence => utiliser une erreur L2
//todo: gérer des matériaux composites
//todo: implémenter solver dynamique
//todo: demander pour le critical d'openmp

/*-----------------todo qu'il reste vraiment à faire------------------*/
//todo: afficher les déplacements, le potentiel en fonction de la distance dans l'élément (une coupe 1D dans le domaine)
//todo: representation formula pour le champ electrique
//todo: créer d'autres .geo avec des géométries plus complexes

//astuce: pour reverse le sens de parcours d'une courbe et inverser les tangentes dans le .geo, utiliser reversemesh

//implémenter openmp pour getBC + résoudre problème temps de calcul (section critical à enlever)
//calcul analytique de G et des gradients
//coupe 2d
//pourquoi bem_domain_1 fonctionne
//comment utiliser les plugins dans le code (continuous color map)
//convergence par type d'élément
//convergence de la solution numérique par rapport à l'annalytique pour le champs éléectrique en fonction du nombre de points de gauss
//temps de calculs
//refaire graphique vectoriel


void element_node_data(std::vector<std::vector<std::vector<double>>> &data, std::vector<std::vector<size_t>> &elementTags2d_domain_BEM, std::vector<int> &elementTypes2d_domain_BEM, std::vector<std::vector<size_t>> &nodeTags2d_domain_BEM, std::vector<elementStruct> &elementVector, std::map<int, double> &nodeMapDomain, bool &show_element_node_data, int &entities_tags_for_physical_group_domain_BEM, std::map<std::size_t, std::pair<bool,bool>> &nodeOnBoundary)
{
    //The first part of the algorithm consists in recovering the tag of the 2d elements of the domain which have at least one node in common with the boundary
    //We initialise the first dimension of the vector that will contain all 2d elements in common with the boundary with the number of element types present in the domain.
    std::vector<std::vector<size_t>> elementTags2d_boundary(elementTypes2d_domain_BEM.size());

    //this loop goes through all the nodes of the boundary (node 1 of each element). We then use getElementsByCoordinates to obtain all the elements in common with the node considered. We only take the 2d elements.
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        double posx = elementVector[i].x1;
        double posy = elementVector[i].y1;
        double posz = 0;
        std::vector<size_t> elementTags;
        gmsh::model::mesh::getElementsByCoordinates(posx,posy,posz,elementTags,2);

        //we loop on the 2d elements which are in common with the considered node of the boundary. Depending on their type, they are added to the vector comprising all the 2d elements common with the boundary.
        for(size_t j = 0; j < elementTags.size(); j++)
        {
            int elementType;
            std::vector<size_t> nodeTags;
            int dim;
            int tag;
            std::vector<int> physicalTags;
            gmsh::model::mesh::getElement(elementTags[j],elementType,nodeTags,dim,tag);
            gmsh::model::getPhysicalGroupsForEntity(dim,tag,physicalTags);

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
    std::map<int, int> tagposition;

    //The first dimension of data corresponds to the number of element types.
    data.resize(elementTags2d_domain_BEM.size());
    for(size_t i = 0; i < elementTags2d_domain_BEM.size(); i++)
    {
        //The second dimension of data corresponds to the number of element tags for on element type.
        data[i].resize(elementTags2d_domain_BEM[i].size());
        std::string elementName;
        int dim;
        int order;
        int numNodes;
        std::vector<double> localNodeCoord;
        int numPrimaryNodes;
        gmsh::model::mesh::getElementProperties(elementTypes2d_domain_BEM[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);

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
        std::string elementName;
        int dim;
        int order;
        int numNodes;
        std::vector<double> localNodeCoord;
        int numPrimaryNodes;
        gmsh::model::mesh::getElementProperties(elementTypes2d_domain_BEM[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
        //getElementProperties was used to get the number of nodes for a specific element
        int nb_nodes_element2d = numNodes;

        //The second loop is on the elementTag in 2d in common with the boundary.
        for(size_t j = 0; j < elementTags2d_boundary[i].size(); j++)
        {
            std::vector<int> elements_commun_tag(0,0);
            std::vector<int> elements_commun_indice(0,0);
            std::vector<int> nodes_commun_tag(0,0);
            int elementType;
            std::vector<size_t> nodeTags_element2d;
            int dim;
            int tag;
            gmsh::model::mesh::getElement(elementTags2d_boundary[i][j],elementType,nodeTags_element2d,dim,tag);
            //The third loop is on the nodes of a specific elementTag in 2d in common with the boundary.
            for(int k = 0; k < nb_nodes_element2d; k++)
            {
                int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                //For each node of the 2d elements, we loop over all the 1d elements of the boundary and look for both the node in common and the elements in common. Indeed, on the boundary, a specific node belongs to two 1d elements.
                for(size_t l = 0; l < elementVector.size(); l++)
                {
                    //condition vérifié si le noeud de l’élément 2d est égal au noeud d’un élément 1d de la frontière.
                    if(elementVector[l].nodeTags[0] == (size_t) loop_element_2d_nodes_tag || elementVector[l].nodeTags[1] == (size_t) loop_element_2d_nodes_tag)
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
                std::cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
                for(size_t m = 0; m < nodes_commun_tag.size(); m++)
                {
                    std::cout << nodes_commun_tag[m] << " ";
                }
                std::cout << std::endl;
                std::cout << "elements commun: " << elements_commun_tag.size() << ": ";
                for(size_t m = 0; m < elements_commun_tag.size(); m++)
                {
                    std::cout << elements_commun_tag[m] << " ";
                }
                std::cout << std::endl;
                std::cout << std::endl;
            }
            //Several geometric configurations are to be considered. Here, the 2d element has 1 node in common with the boundary and must therefore touch 2 of its elements.
            if(nodes_commun_tag.size() == 1)
            {
                if(elements_commun_tag.size() == 2)
                {
                    //We loop over the nodes of the 2d element, find the one that is common to the boundary and assign it the average value of the potential of the 2 elements it touches. Indeed the potential is discontinuous from one element of the boundary to the other.
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        if(loop_element_2d_nodes_tag == nodes_commun_tag[0])
                        {
                            double phi = 0;
                            phi = elementVector[elements_commun_indice[0]].bcValue[0];
                            phi = phi + elementVector[elements_commun_indice[1]].bcValue[0];
                            phi = phi/2.0;
                            data[i][tagposition[elementTags2d_boundary[i][j]]][k] = phi;
                            if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 0)
                            {
                                std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from common nodes: " << nodes_commun_tag.size() << " and common elements: " << elements_commun_tag.size() << std::endl;
                            }
                            else
                            {
                                nodeOnBoundary[loop_element_2d_nodes_tag].second = 1;
                            }
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
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        for(size_t l = 0; l < elements_commun_tag.size(); l++)
                        {
                            //Among the 3 elements of the boundary in common, we look for the one which has 2 nodes in common with the 2d element. We look for the side and therefore the 1d element that they have in common.
                            int compteur_nodes = 0;
                            for(size_t m = 0; m < nodes_commun_tag.size(); m++)
                            {
                                if(elementVector[elements_commun_indice[l]].nodeTags[0] == (size_t) nodes_commun_tag[m]||elementVector[elements_commun_indice[l]].nodeTags[1] == (size_t) nodes_commun_tag[m])
                                {
                                    compteur_nodes++;
                                }
                            }
                            if(compteur_nodes == 2)
                            {
                                //Once we have found the element in common, the node of the 2d element on which we are looping must be the right one
                                if((size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]|| (size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                                {
                                    //The constant value of the boundary element in common is assigned to the 2 nodes of the 2d element that touches the boundary.
                                    data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                    if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 0)
                                    {
                                        std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from common nodes: " << nodes_commun_tag.size() << " and common elements: " << elements_commun_tag.size() << std::endl;
                                    }
                                    else
                                    {
                                        nodeOnBoundary[loop_element_2d_nodes_tag].second = 1;
                                    }
                                }
                                //In the case of second order elements, the value of the common 1d boundary element is assigned to the middle node.
                                if(order == 2)
                                {
                                    if((size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].midNodeTags)
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                        if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 0)
                                        {
                                            std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from common nodes: " << nodes_commun_tag.size() << " and common elements: " << elements_commun_tag.size() << std::endl;
                                        }
                                        else
                                        {
                                            nodeOnBoundary[loop_element_2d_nodes_tag].second = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else if(elements_commun_tag.size() == 4)
                {
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        for(size_t l = 0; l < elements_commun_tag.size(); l++)
                        {
                            //This geometry behaves almost like the first case (1 node, 2 elements) except that it is true for 2 different nodes of the 2d element. We loop over the elements of the boundary. If one of its nodes is common with the boundary, it is assigned the value of the potential of this first element. When the second one comes, it is averaged with the potential already present on the node.
                            if((size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]|| (size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                            {
                                if(data[i][tagposition[elementTags2d_boundary[i][j]]][k] == HUGE_VAL)
                                {
                                    data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                }
                                else
                                {
                                    data[i][tagposition[elementTags2d_boundary[i][j]]][k] = (data[i][tagposition[elementTags2d_boundary[i][j]]][k] + elementVector[elements_commun_indice[l]].bcValue[0])/2.0;
                                }
                                if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 0)
                                {
                                    std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from common nodes: " << nodes_commun_tag.size() << " and common elements: " << elements_commun_tag.size() << std::endl;
                                }
                                else
                                {
                                    nodeOnBoundary[loop_element_2d_nodes_tag].second = 1;
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
                    for(int k = 0; k < nb_nodes_element2d; k++)
                    {
                        int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                        for(size_t l = 0; l < elements_commun_tag.size(); l++)
                        {
                            int compteur_nodes = 0;
                            for(size_t m = 0; m < nodes_commun_tag.size(); m++)
                            {
                                if(elementVector[elements_commun_indice[l]].nodeTags[0] == (size_t) nodes_commun_tag[m]||elementVector[elements_commun_indice[l]].nodeTags[1] == (size_t) nodes_commun_tag[m])
                                {
                                    compteur_nodes++;
                                }
                            }
                            if(compteur_nodes == 2)
                            {
                                if((size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[0]|| (size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].nodeTags[1])
                                {
                                    if(data[i][tagposition[elementTags2d_boundary[i][j]]][k] == HUGE_VAL)
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                    }
                                    else
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = (data[i][tagposition[elementTags2d_boundary[i][j]]][k] + elementVector[elements_commun_indice[l]].bcValue[0])/2.0;
                                    }
                                    if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 0)
                                    {
                                        std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from common nodes: " << nodes_commun_tag.size() << " and common elements: " << elements_commun_tag.size() << std::endl;
                                    }
                                    else
                                    {
                                        nodeOnBoundary[loop_element_2d_nodes_tag].second = 1;
                                    }
                                }
                                //In the case of second order elements, the value of the common 1d boundary element is assigned to the middle node.
                                if(order == 2)
                                {
                                    if((size_t) loop_element_2d_nodes_tag == elementVector[elements_commun_indice[l]].midNodeTags)
                                    {
                                        data[i][tagposition[elementTags2d_boundary[i][j]]][k] = elementVector[elements_commun_indice[l]].bcValue[0];
                                        if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 0)
                                        {
                                            std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from common nodes: " << nodes_commun_tag.size() << " and common elements: " << elements_commun_tag.size() << std::endl;
                                        }
                                        else
                                        {
                                            nodeOnBoundary[loop_element_2d_nodes_tag].second = 1;
                                        }
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
        std::string elementName;
        int dim;
        int order;
        int numNodes;
        std::vector<double> localNodeCoord;
        int numPrimaryNodes;
        gmsh::model::mesh::getElementProperties(elementTypes2d_domain_BEM[i], elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
        int nb_nodes_element2d = numNodes;
        for(size_t j = 0; j < elementTags2d_domain_BEM[i].size(); j++)
        {
            int elementType;
            std::vector<size_t> nodeTags_element2d;
            int dim;
            int tag;
            gmsh::model::mesh::getElement(elementTags2d_domain_BEM[i][j],elementType,nodeTags_element2d,dim,tag);
            for(int k = 0; k < nb_nodes_element2d; k++)
            {
                int loop_element_2d_nodes_tag = nodeTags_element2d[k];
                //je laisse la condition hugeval pour l'instant car elle est plus "gentille". Le code se plantera moins facilement mais il y aura quand même un message dans la console si on a des noeuds qui ont un potentiel anormal
                if(data[i][j][k] == HUGE_VAL)
                {
                    if(nodeMapDomain.count(loop_element_2d_nodes_tag))
                    {
                        data[i][j][k] = nodeMapDomain[loop_element_2d_nodes_tag];
                    }
                    else
                    {
                        std::cout << "error: node not found in nodeMapDomain code BEM" << std::endl;
                    }
                    if(nodeOnBoundary[loop_element_2d_nodes_tag].first == 1 || nodeOnBoundary[loop_element_2d_nodes_tag].second == 1)
                    {
                        std::cerr << "error node: " << loop_element_2d_nodes_tag << " is considered on the boundary and in the BEM domain" << " output from the use of nodemapdomain" << std::endl;
                    }
                }
            }
        }
    }
}

void dispNodeCoordFun(const std::vector<size_t> &allNodeTags, const std::vector<double> &allCoord)
{
    for(size_t i = 0; i < allNodeTags.size(); i++)
    {
        std::cout << "Node '" << allNodeTags[i] << "': (" << allCoord[3*i] << "; " << allCoord[3*i + 1] << "; " << allCoord[3*i + 2] << ")\n";
    }
    std::cout << std::endl;
}

void dispMatricesFun(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const Eigen::MatrixXd &c, const Eigen::MatrixXd &b, const Eigen::MatrixXd &x)
{
    std::cout << "A matrix:\n";
    std::cout << A << std::endl;
    std::cout << "\nB matrix:\n";
    std::cout << B << std::endl;
    std::cout << "\nc vector:\n";
    std::cout << c << std::endl;
    std::cout << "\nb vector:\n";
    std::cout << b << std::endl;
    std::cout << "\nx vector:\n";
    std::cout << x << std::endl;
}

void dispValuesFun(const std::vector<elementStruct> &elementVector)
{
    for(size_t i = 0; i < elementVector.size(); i++)
    {
        std::cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
    }
}

void getBC(std::map<size_t, int> &allNodeIndex, const std::vector<double> &allCoord, std::vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy, std::vector<size_t> &NodeTags2d_domain, std::vector<double> &NodeCoord2d_domain, std::vector<int> &elementTypes2d_domain_BEM, std::vector<std::vector<size_t>> &elementTags2d_domain_BEM, std::vector<std::vector<size_t>> &elnodeTags2d_domain_BEM, int &entities_tags_for_physical_group_domain_BEM, std::vector<int> &entities_tags_2d_domain_BEM, std::map<int, std::pair<double, double>> &boundaryDisplacementMap, std::map<std::size_t, std::pair<bool,bool>> &nodeOnBoundary, const int domainTag)
{
//    start = std::chrono::system_clock::now();

//    end = std::chrono::system_clock::now();
//    function_name = "getBC";
//    display_time(start, end, function_name, show_time);

    // pour moi nodeTags2d_domain inutile, ainsi que nodeCoords2d_domain
    // Gets the dim and the tags of the physical groups
    //std::chrono::time_point<std::chrono::system_clock> start, end;
    //std::string function_name;
    //bool show_time = 1;

    //start = std::chrono::system_clock::now();
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    // Gets the boundary conditions from onelab and create a map: "name of the physical group" -> ("dirichlet"/"neumann", value)
    std::map<std::string, std::pair<std::string, double>> bcMap;

    std::vector<std::string> bcKeys;
    gmsh::onelab::getNames(bcKeys, "Boundary Conditions.+");

    std::string domainName = "BEM_domain_" + std::to_string(domainTag);
    for(auto &key : bcKeys)
    {
        std::vector<double> values;
        gmsh::onelab::getNumber(key, values);

        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(getline(ss, word, '/'))
            words.push_back(word);

        // focus only on dirichlet or neumann boundary conditions -> associated to BEM // CHECK BON DOMAIN
        if(words.size() == 4 && (words[3].compare("neumann") == 0 || words[3].compare("dirichlet") == 0) && words[2].compare(domainName) == 0)
        {
            bcMap[words[1]] = {words[3], values[0]}; // doit etre particularisée
        }
    }

//    end = std::chrono::system_clock::now();
//    function_name = "get boundary conditions";
//    display_time(start, end, function_name, show_time);

    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    for (size_t i = 0; i < physicalGroupsDimTags.size(); i++)
    {
        int physicalGroupsDim = physicalGroupsDimTags[i].first;
        int physicalGroupsTag = physicalGroupsDimTags[i].second;
        std::string name;
        gmsh::model::getPhysicalName(physicalGroupsDim, physicalGroupsTag, name);

        if(name == domainName)
        {
            //start = std::chrono::system_clock::now();
            gmsh::model::mesh::getNodesForPhysicalGroup(physicalGroupsDim,physicalGroupsTag,NodeTags2d_domain,NodeCoord2d_domain);

            entities_tags_for_physical_group_domain_BEM = physicalGroupsTag;
            gmsh::model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entities_tags_2d_domain_BEM);

            for(size_t j=0; j< entities_tags_2d_domain_BEM.size(); j++)
            {
                std::vector<int> elementTypes;
                std::vector<std::vector<size_t> > elementTags;
                std::vector<std::vector<size_t> > nodeTags;
                gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, physicalGroupsDim, entities_tags_2d_domain_BEM[j]);
                elementTypes2d_domain_BEM.insert(elementTypes2d_domain_BEM.end(), elementTypes.begin(), elementTypes.end());
                elementTags2d_domain_BEM.insert(elementTags2d_domain_BEM.end(), elementTags.begin(), elementTags.end());
            }
//            end = std::chrono::system_clock::now();
//            function_name = "get domain mesh";
//            display_time(start, end, function_name, show_time);
        }

        //s'il ne trouve pas name dans la map, find retourne map.end() => condition pour continue
        if(bcMap.find(name) == bcMap.end())
        {
            continue;
        }

        // Retrieves the type and the value of the boundary condition of this physical group.
        // /!\ Should check that this physical groups was well defined in the list of "Boundary Condition" !
        //start = std::chrono::system_clock::now();
        std::string bcType = bcMap[name].first;
        double bcValue = bcMap[name].second;

        // Prints the physical groups name, dim, tag, bcType and bcValue (if dispHierarchy != 0).
        if(dispHierarchy)
            std::cout << "Physical group '" << name << "' (" << physicalGroupsDim << "D) with tag '" << physicalGroupsTag << "' has a '" << bcType << "' type and the value " << bcValue << "\n";

        std::vector<int> entitiesTags;
        gmsh::model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entitiesTags);

        // devrait passer seulement sur les physical groups qui sont à la frontière de BEM_domain_i + Neumann ou Dirichlet
        for (size_t k = 0; k < entitiesTags.size(); k++)
        {
            // Print the tag of the entities of each physical group (if dispHierarchy != 0).
            if(dispHierarchy)
                std::cout << "\t contains entity '" << entitiesTags[k] << "' ";

            std::vector<size_t> elementTags1;
            std::vector<size_t> nodeTags1;
            std::vector<size_t> elementTags8;
            std::vector<size_t> nodeTags8;

            gmsh::model::mesh::getElementsByType(1, elementTags1, nodeTags1, entitiesTags[k]);
            gmsh::model::mesh::getElementsByType(8, elementTags8, nodeTags8, entitiesTags[k]);

            // Print the number of element of the given entity (if dispHierarchy != 0).
            if(dispHierarchy)
                std::cout << "which has " << elementTags1.size() + elementTags8.size() << " elements:\n";

            //start = std::chrono::system_clock::now();
            for(size_t el = 0; el < elementTags1.size(); el++)
            {

                int elTag = elementTags1[el];

                elementStruct element;
                element.tag = elTag;
                element.nodeTags[0] = nodeTags1[2*el];
                element.nodeTags[1] = nodeTags1[2*el + 1];
                element.midNodeTags = 0;

                std::pair<double, double> nodalDisp1 = {0.0, 0.0}, nodalDisp2 = {0.0, 0.0};
                if(boundaryDisplacementMap.count(element.nodeTags[0]))
                    nodalDisp1 = boundaryDisplacementMap[element.nodeTags[0]];
                if(boundaryDisplacementMap.count(element.nodeTags[1]))
                    nodalDisp2 = boundaryDisplacementMap[element.nodeTags[1]];
                computeGeometricParam(element, allNodeIndex, allCoord, nodalDisp1, nodalDisp2);

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
                    std::cout << "\t\t-Element '" << elTag << "' (" << element.bcType << "):\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[0]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[0]] + 1] << ")\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[1]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[1]] + 1] << ")\n";
                }
            }
//            end = std::chrono::system_clock::now();
//            function_name = "get elements";
//            display_time(start, end, function_name, show_time);
            for(size_t el = 0; el < elementTags8.size(); el++)
            {
                int elTag = elementTags8[el];

                elementStruct element;
                element.tag = elTag;
                element.nodeTags[0] = nodeTags8[3*el];
                element.nodeTags[1] = nodeTags8[3*el + 1];
                element.midNodeTags = nodeTags8[3*el + 2];

                std::pair<double, double> nodalDisp1 = {0.0, 0.0}, nodalDisp2 = {0.0, 0.0};
                if(boundaryDisplacementMap.count(element.nodeTags[0]))
                    nodalDisp1 = boundaryDisplacementMap[element.nodeTags[0]];
                if(boundaryDisplacementMap.count(element.nodeTags[1]))
                    nodalDisp2 = boundaryDisplacementMap[element.nodeTags[1]];
                computeGeometricParam(element, allNodeIndex, allCoord, nodalDisp1, nodalDisp2);

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
                    std::cout << "\t\t-Element '" << elTag << "' (" << element.bcType << "):\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[0]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[0]] + 1] << ")\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*allNodeIndex[element.nodeTags[1]]] << "; " << allCoord[3*allNodeIndex[element.nodeTags[1]] + 1] << ")\n";
                }
            }
            if(dispHierarchy)
                std::cout << "\n";
        }
//        end = std::chrono::system_clock::now();
//        function_name = "get boundary elements";
//        display_time(start, end, function_name, show_time);
    }
    for(std::size_t i = 0; i < NodeTags2d_domain.size(); i++)
    {
        nodeOnBoundary[NodeTags2d_domain[i]].first = 0;
        nodeOnBoundary[NodeTags2d_domain[i]].second = 0;
    }
    for(std::size_t i = 0; i < elementVector.size(); i++)
    {
        nodeOnBoundary[elementVector[i].nodeTags[0]].first = 1;
        nodeOnBoundary[elementVector[i].nodeTags[1]].first = 1;
        if(elementVector[i].midNodeTags != 0)
        {
            nodeOnBoundary[elementVector[i].midNodeTags].first = 1;
        }
    }
}

void computeGeometricParam(elementStruct &element, std::map<std::size_t, int> &allNodeIndex, const std::vector<double> &allCoord, std::pair<double, double> &nodalDisp1, std::pair<double, double> &nodalDisp2)
{
    int i1Tag = element.nodeTags[0];
    int i2Tag = element.nodeTags[1];

    int i1Index = allNodeIndex[i1Tag];
    int i2Index = allNodeIndex[i2Tag];

    // Get the coordinates of the nodes from their tag.
    element.x1 = allCoord[3*i1Index] + nodalDisp1.first;
    element.y1 = allCoord[3*i1Index + 1] + nodalDisp1.second;
    element.x2 = allCoord[3*i2Index] + nodalDisp2.first;
    element.y2 = allCoord[3*i2Index + 1] + nodalDisp2.second;

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


double computeG(const double &x, const double &y, const elementStruct &element, const std::string &integrationType)
{
    // Computation of the integration points (parametric) and weights.
    std::vector<double> localCoord;
    std::vector<double> weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

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

    double g = gSum * element.length / (4*M_PI);

    return g;
}

void displayData(const std::string viewName, const std::string modelName, const std::string fileName, const std::vector<size_t> &tags, const std::vector<std::vector<double>> &data)
{
    int viewTag = gmsh::view::add(viewName);
    gmsh::view::addModelData(viewTag, 0, modelName, "ElementData", tags, data);
    gmsh::view::write(viewTag, fileName);

    gmsh::option::setNumber("View[0].IntervalsType", 3);
    // option::setNumber("View[0].AdaptVisualizationGrid", 1); // Problem with this command !
    gmsh::option::setNumber("View[0].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[0].TargetError", -0.0001);
}

void compute_internal_phi(std::vector<double> &allCoord, std::vector<size_t> &NodeTags2d_domain, std::map<size_t, int> &allNodeIndex, std::vector<elementStruct> &elementVector, std::map<int,double> &nodeMapDomain, const std::string integrationType, std::map<std::size_t, std::pair<bool,bool>> &nodeOnBoundary, bool &solution_type)
{
    double x = 0;
    double y = 0;
    double phi = 0;
    double h = 0;
    double g = 0;
    std::vector<double> vectorForNodeMap(NodeTags2d_domain.size(),0);

#pragma omp parallel for schedule(static) private(x,y,phi,h,g) shared(NodeTags2d_domain,allCoord,elementVector,nodeMapDomain,nodeOnBoundary,solution_type)
    for(std::size_t i = 0; i < NodeTags2d_domain.size(); i++)
    {
        if(nodeOnBoundary[NodeTags2d_domain[i]].first == 0)
        {
            // Coordinates of the node
            x = allCoord[3*allNodeIndex[NodeTags2d_domain[i]]];
            y = allCoord[3*allNodeIndex[NodeTags2d_domain[i]] + 1];

            phi = 0;
            for(size_t j = 0; j < elementVector.size(); j++)
            {
                h = computeH(x, y, elementVector[j]);
                if(solution_type == 0)
                {
                    g = computeG(x, y, elementVector[j], integrationType);
                }
                else
                {
                    g = computeG_analytique(x, y, elementVector[j]);
                }

                phi = phi + h*elementVector[j].bcValue[0] - g*elementVector[j].bcValue[1];
            }
            vectorForNodeMap[i] = phi;
        }
    }
    for(std::size_t i = 0; i < NodeTags2d_domain.size(); i++)
    {
        nodeMapDomain[NodeTags2d_domain[i]] = vectorForNodeMap[i];
    }
}

void computeData(const std::vector<elementStruct> &elementVector, std::vector<size_t> &tags, std::vector<std::vector<double>> &data, const std::string integrationType, std::vector<int> &entities_tags_2d_domain_BEM, bool &solution_type)
{
    double x = 0;
    double y = 0;
    double h = 0;
    double g = 0;
    double phi = 0;
    std::vector<double> barycenters;
    size_t elementCount = 0;

    for(size_t i = 0; i < entities_tags_2d_domain_BEM.size(); i++)
    {
        std::vector<int> elementTypes_for_entity;
        std::vector<std::vector<size_t> > elementTags_for_entity;
        std::vector<std::vector<size_t> > nodeTags_for_entity;
        gmsh::model::mesh::getElements(elementTypes_for_entity, elementTags_for_entity, nodeTags_for_entity, 2, entities_tags_2d_domain_BEM[i]);

        for(size_t l = 0; l < elementTypes_for_entity.size(); l++)
        {
            gmsh::model::mesh::getBarycenters(elementTypes_for_entity[l], entities_tags_2d_domain_BEM[i], false, false, barycenters);

#pragma omp parallel for schedule(static) private(x,y,phi,h,g) shared(elementTags_for_entity,elementVector,tags,data,elementCount,barycenters,solution_type)
            for(size_t j = 0; j < elementTags_for_entity[l].size(); j++)
            {
                x = barycenters[3*j];
                y = barycenters[3*j+1];

                phi = 0;
                for(size_t k = 0; k < elementVector.size(); k++)
                {
                    h = computeH(x, y, elementVector[k]);
                    if(solution_type == 0)
                    {
                        g = computeG(x, y, elementVector[k], integrationType);
                    }
                    else
                    {
                        g = computeG_analytique(x, y, elementVector[k]);
                    }
                    phi += h*elementVector[k].bcValue[0] - g*elementVector[k].bcValue[1];
                }
                // Potential at the node is phi.
                tags[elementCount + j] = elementTags_for_entity[l][j];
                data[elementCount + j][0] = phi;
            }
            elementCount = elementCount + elementTags_for_entity[l].size();
        }
    }
}

void disp_element_node_data(std::vector<int> &nodes_commun_tag, std::vector<int> &elements_commun_tag)
{
    std::cout << "problème dans la géométrie, noeuds communs = " << nodes_commun_tag.size() << std::endl;
    std::cout << "noeuds commun: "<< nodes_commun_tag.size() << ": ";
    for(size_t m = 0; m < nodes_commun_tag.size(); m++)
    {
        std::cout << nodes_commun_tag[m] << " ";
    }
    std::cout << std::endl;
    std::cout << "elements commun: " << elements_commun_tag.size() << ": ";
    for(size_t m = 0; m < elements_commun_tag.size(); m++)
    {
        std::cout << elements_commun_tag[m] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void display_time(std::chrono::time_point<std::chrono::system_clock> &start, std::chrono::time_point<std::chrono::system_clock> &end, std::string &function_name, bool &show_time)
{
    //displays the time of the desired code section
    if(show_time)
    {
        std::cout << std::endl;
        std::cout << function_name << std::endl;
        std::cout << "Elapsed time in nanoseconds: "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << " ns" << std::endl;

        std::cout << "Elapsed time in microseconds: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " µs" << std::endl;

        std::cout << "Elapsed time in milliseconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;

        std::cout << "Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " sec";
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

void convergenceFun(std::vector<int> &entities_tags_2d_domain_BEM, std::vector<size_t> &tags, std::vector<std::vector<double>> &phi)
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

double computeGradHx(const double &x, const double &y, const elementStruct &element, const std::string &integrationType)
{
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;

    double dx = xjp1-xj;
    double dy = yjp1-yj;
    double nx = dy;
    double ny = -dx;
    double norm = sqrt(nx*nx+ny*ny);
    nx = nx / norm;
    ny = ny / norm;

    // Computation of the integration points (parametric) and weights.
    std::vector<double> localCoord;
    std::vector<double> weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Computation of the sum (of integrand evaluated at integration point * weight)
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        // Parametric coordinate
        double integration_variable = localCoord[3*k];
        // Global coordinates
        double hX = element.midX + element.deltaX*integration_variable/2.0;
        double hY = element.midY + element.deltaY*integration_variable/2.0;

        // Distance between mid point of element 'i' and integration point.
        double r = sqrt(pow(x - hX, 2.0) + pow(y - hY, 2.0));
        gSum += (-((1.0/2.0*pow(r, 2) - pow((xi - hX),2))/(pow(r, 4)))*nx + ((xi - hX)*(yi - hY)*ny)/(pow(r, 4))) * weights[k];
    }

    double gradHx = gSum * element.length / (2*M_PI);

    return gradHx;
}
double computeGradHy(const double &x, const double &y, const elementStruct &element, const std::string &integrationType)
{
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;

    double dx = xjp1-xj;
    double dy = yjp1-yj;
    double nx = dy;
    double ny = -dx;
    double norm = sqrt(nx*nx+ny*ny);
    nx = nx / norm;
    ny = ny / norm;

    // Computation of the integration points (parametric) and weights.
    std::vector<double> localCoord;
    std::vector<double> weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Computation of the sum (of integrand evaluated at integration point * weight)
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        // Parametric coordinate
        double integration_variable = localCoord[3*k];
        // Global coordinates
        double hX = element.midX + element.deltaX*integration_variable/2.0;
        double hY = element.midY + element.deltaY*integration_variable/2.0;

        // Distance between mid point of element 'i' and integration point.
        double r = sqrt(pow(xi - hX, 2.0) + pow(yi - hY, 2.0));
        gSum += (-((1.0/2.0*pow(r, 2) - pow((yi - hY),2))/(pow(r, 4)))*ny + ((xi - hX)*(yi - hY)*nx)/(pow(r, 4))) * weights[k];
    }

    double gradHy = gSum * element.length / (2*M_PI);

    return gradHy;
}
double computeGradGx(const double &x, const double &y, const elementStruct &element, const std::string &integrationType)
{
    // Computation of the integration points (parametric) and weights.
    std::vector<double> localCoord;
    std::vector<double> weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Computation of the sum (of integrand evaluated at integration point * weight)
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        // Parametric coordinate
        double integration_variable = localCoord[3*k];
        // Global coordinates
        double gX = element.midX + element.deltaX*integration_variable/2.0;
        double gY = element.midY + element.deltaY*integration_variable/2.0;

        // Distance between mid point of element 'i' and integration point.
        double r = sqrt(pow(x - gX, 2.0) + pow(y - gY, 2.0));
        gSum += ((x-gX)/(pow(r, 2))) * weights[k];
    }

    double gradGx = gSum * element.length / (4*M_PI);

    return gradGx;
}
double computeGradGy(const double &x, const double &y, const elementStruct &element, const std::string &integrationType)
{
    // Computation of the integration points (parametric) and weights.
    std::vector<double> localCoord;
    std::vector<double> weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Computation of the sum (of integrand evaluated at integration point * weight)
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        // Parametric coordinate
        double integration_variable = localCoord[3*k];
        // Global coordinates
        double gX = element.midX + element.deltaX*integration_variable/2.0;
        double gY = element.midY + element.deltaY*integration_variable/2.0;

        // Distance between mid point of element 'i' and integration point.
        double r = sqrt(pow(x - gX, 2.0) + pow(y - gY, 2.0));
        gSum += ((y-gY)/(pow(r, 2))) * weights[k];
    }

    double gradGy = gSum * element.length / (4*M_PI);

    return gradGy;
}

double computeGradHx_analytique(const double &x, const double &y, const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double dx = xjp1-xj;
    double dy = yjp1-yj;
    double nx = dy;
    double ny = -dx;
    double norm = sqrt(nx*nx+ny*ny);
    nx = nx / norm;
    ny = ny / norm;

    double gradHx_analytique = (lj*(nx*(pow(xi,2) + xj*xjp1 - xi*(xj + xjp1) - (yi - yj)*(yi - yjp1)) - ny*(xj*yi + xjp1*yi - xjp1*yj - xj*yjp1 + xi*(-2*yi + yj + yjp1))))/(2.*Pi*(pow(xi - xj,2) + pow(yi - yj,2))*(pow(xi - xjp1,2) + pow(yi - yjp1,2)));
    return gradHx_analytique;
}
double computeGradHy_analytique(const double &x, const double &y, const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double dx = xjp1-xj;
    double dy = yjp1-yj;
    double nx = dy;
    double ny = -dx;
    double norm = sqrt(nx*nx+ny*ny);
    nx = nx / norm;
    ny = ny / norm;

    double gradHy_analytique = -(lj*(ny*(pow(xi,2) + xj*xjp1 - xi*(xj + xjp1) - (yi - yj)*(yi - yjp1)) + nx*(xj*yi + xjp1*yi - xjp1*yj - xj*yjp1 + xi*(-2*yi + yj + yjp1))))/(2.*Pi*(pow(xi - xj,2) + pow(yi - yj,2))*(pow(xi - xjp1,2) + pow(yi - yjp1,2)));
    return gradHy_analytique;
}
double computeGradGx_analytique(const double &x, const double &y, const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double gradGx_analytique = (lj*(-2*(yj - yjp1)*(atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1))) + atan(((xi - xjp1)*(-xj + xjp1) + (yi - yjp1)*(-yj + yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1)))) - (xj - xjp1)*(log(pow(xi - xj,2) + pow(yi - yj,2)) - log(pow(xi - xjp1,2) + pow(yi - yjp1,2)))))/(4.*Pi*(pow(xj - xjp1,2) + pow(yj - yjp1,2)));
    return gradGx_analytique;
}
double computeGradGy_analytique(const double &x, const double &y, const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double gradGy_analytique = (lj*(2*(xj - xjp1)*(atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1))) + atan(((xi - xjp1)*(-xj + xjp1) + (yi - yjp1)*(-yj + yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1)))) - (yj - yjp1)*(log(pow(xi - xj,2) + pow(yi - yj,2)) - log(pow(xi - xjp1,2) + pow(yi - yjp1,2)))))/(4.*Pi*(pow(xj - xjp1,2) + pow(yj - yjp1,2)));
    return gradGy_analytique;
}

double computeG_analytique(const double &x, const double &y, const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double g_analytique = (lj*(-2*(pow(xj - xjp1,2) + pow(yj - yjp1,2)) + 2*(-(xj*yi) + xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*yjp1)* atan(((xi - xjp1)*(xj - xjp1) + (yi - yjp1)*(yj - yjp1))/(-(xj*yi) + xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*yjp1)) + 2*(-(xjp1*yi) - xi*yj + xjp1*yj + xj*(yi - yjp1) + xi*yjp1)* atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1))/ (xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1))) + (-((xi - xj)*(xj - xjp1)) - (yi - yj)*(yj - yjp1))* log(pow(xi - xj,2) + pow(yi - yj,2)) + ((xi - xjp1)*(xj - xjp1) + (yi - yjp1)*(yj - yjp1))* log(pow(xi - xjp1,2) + pow(yi - yjp1,2))))/(4.*Pi*(pow(xj - xjp1,2) + pow(yj - yjp1,2)));
    return g_analytique;
}

void electricField(int &nb2dElements, std::vector<int> &entities_tags_2d_domain_BEM, const std::vector<elementStruct> &elementVector, std::vector<size_t> &tags, const std::string &integrationType, std::vector<std::vector<double>> &data_Electric_field, bool &solution_type)
{
    double x = 0;
    double y = 0;
    double grad_phi_x = 0;
    double grad_phi_y = 0;
    double h_x = 0;
    double g_x = 0;
    double h_y = 0;
    double g_y = 0;
    std::vector<double> barycenters;
    size_t elementCount = 0;

    for(size_t i = 0; i < entities_tags_2d_domain_BEM.size(); i++)
    {
        std::vector<int> elementTypes_for_entity;
        std::vector<std::vector<size_t> > elementTags_for_entity;
        std::vector<std::vector<size_t> > nodeTags_for_entity;
        gmsh::model::mesh::getElements(elementTypes_for_entity, elementTags_for_entity, nodeTags_for_entity, 2, entities_tags_2d_domain_BEM[i]);

        for(size_t l = 0; l < elementTypes_for_entity.size(); l++)
        {
            gmsh::model::mesh::getBarycenters(elementTypes_for_entity[l], entities_tags_2d_domain_BEM[i], false, false, barycenters);
#pragma omp parallel for schedule(static) private(x,y,grad_phi_x,grad_phi_y,h_x,g_x,h_y,g_y) shared(elementTags_for_entity,elementVector,tags,data_Electric_field,elementCount,barycenters,solution_type,integrationType)
            for(size_t j = 0; j < elementTags_for_entity[l].size(); j++)
            {
                x = barycenters[3*j];
                y = barycenters[3*j+1];

                grad_phi_x = 0;
                grad_phi_y = 0;
                for(size_t k = 0; k < elementVector.size(); k++)
                {
                    if(solution_type == 0)
                    {
                        h_x = computeGradHx(x, y, elementVector[k], integrationType);
                        g_x = computeGradGx(x, y, elementVector[k], integrationType);
                        h_y = computeGradHy(x, y, elementVector[k], integrationType);
                        g_y = computeGradGy(x, y, elementVector[k], integrationType);
                    }
                    else
                    {
                        h_x = computeGradHx_analytique(x, y, elementVector[k]);
                        g_x = computeGradGx_analytique(x, y, elementVector[k]);
                        h_y = computeGradHy_analytique(x, y, elementVector[k]);
                        g_y = computeGradGy_analytique(x, y, elementVector[k]);
                    }
                    grad_phi_x = grad_phi_x + h_x*elementVector[k].bcValue[0] - g_x*elementVector[k].bcValue[1];
                    grad_phi_y = grad_phi_y + h_y*elementVector[k].bcValue[0] - g_y*elementVector[k].bcValue[1];
                }
                tags[elementCount + j] = elementTags_for_entity[l][j];
                data_Electric_field[elementCount + j][0] = -grad_phi_x;
                data_Electric_field[elementCount + j][1] = -grad_phi_y;
                data_Electric_field[elementCount + j][2] = 0;
            }
            elementCount = elementCount + elementTags_for_entity[l].size();
        }
    }
}
