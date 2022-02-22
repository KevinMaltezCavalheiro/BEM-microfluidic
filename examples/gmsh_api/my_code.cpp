// This program takes a .geo file as argument, generates the mesh of the 
// surfaces and displays the result in the Gmsh window.
// 
// This example introduces:
// - the gmsh.h header
// - initialize()/finalize() functions
// - C++ namespaces / standard library (std::cout)
//
// How to build/run? (windows) 
//   [in examples/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code1_loadgeo.exe ..\rectangle.geo

#include <gmsh.h>       // gmsh main header (READ THIS FILE)
#include <iostream>     // for std::cout
#include <map>
#include <sstream>  // for std::stringstream
#include <Eigen/Dense> // Eigen library

int main(int argc, char **argv)
{
    /*--------------INIT----------------*/
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    gmsh::initialize(); // first thing to do before using any gmsh SDK function

    gmsh::open(argv[1]); // similar to "File / Open" in the GUI
    gmsh::model::mesh::generate(2);  // similar to "Mesh / 2D" in the GUI

    /*--------get physical groups (more particularly: the domain)--------------*/
    std::map<std::string, std::pair<int, int>> groups;
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    /*-----------get physical properties of domain----------------------*/
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Materials).+");
    double E;
    double nu;
    for (auto &key : keys)
    {
        // get corresponding value
        std::vector<double> value;
        gmsh::onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("Young") == 0){
                E = value[0];
            }
            if(words[2].compare("Poisson") == 0){
                nu = value[0];
            }
        }
    }
    // std::cout << "E = " << E << "\n";
    // std::cout << "nu = " << nu << "\n";

    /*-----------------H matrix------------------*/
    Eigen::Matrix<double, 3, 3> H_matrix;
    H_matrix(0,0) = 1.;
    H_matrix(1,0) = nu;
    H_matrix(2,0) = 0.;
    H_matrix(0,1) = nu;
    H_matrix(1,1) = 1.;
    H_matrix(2,1) = 0.;
    H_matrix(0,2) = 0.;
    H_matrix(1,2) = 0.;
    H_matrix(2,2) = (1-nu)/2;
    H_matrix = E/(1-nu*nu)*H_matrix;
    // std::cout << H_matrix << '\n'; 

    /*------------loop over elements of domain------------------------*/
    std::string groupname = "domain";

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        std::cerr << "Group '" << groupname << "' does not exist!\n";
        return 1;
    }

    // get entities of the domain
    std::vector<int> tags;
    gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);

    std::cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        std::cout << tags[i] << " ";
    std::cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        std::cout << "Entity (" << dim << "," << tags[k] << "):\n";
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            std::cout << "\telement type " << elementTypes[i] << ":\n";
            std::cout << "\t\telements: ";
            for (int j = 0; j < elementTags[i].size(); ++j)
                std::cout << elementTags[i][j] << " ";
            std::cout << '\n';

            // getting element properties
            std::string elementName;
            int dim, order, numNodes, numPrimaryNodes;
            std::vector<double> localNodeCoord;
            gmsh::model::mesh::getElementProperties(elementTypes[i],
                                            elementName, dim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);
            
            // get Gauss points coordinates and weights
            std::vector<double> localCoords, weights;
            gmsh::model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    "CompositeGauss2", // "Gauss2"
                                                    localCoords, weights);
            std::cout << "\t" << weights.size() << " integration points\n";

            std::cout << "\tweights: ";
            for (int j = 0; j < weights.size(); ++j)
                std::cout << weights[j] << " ";
            std::cout << '\n';

            std::cout << "\tlocalCoords (gmsh): ";
            for (int j = 0; j < localCoords.size(); ++j)
                std::cout << localCoords[j] << " ";
            std::cout << '\n';

            // convert to Eigen format
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
                localCoords_E(&localCoords[0], localCoords.size()/3, 3);
            std::cout << "\tlocalCoords (Eigen):\n";
            Eigen::IOFormat fmt(4, 0, ", ", "\n", "\t\t[", "]");
            std::cout << localCoords_E.format(fmt) << '\n';

            // get determinants and Jacobians
            std::vector<double> jacobians, determinants, coords;
            gmsh::model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords);
            std::cout << "\tGot " << determinants.size() 
                    << " Jacobians for type " << elementTypes[i] << "\n";

            std::cout << "\tdeterminants (gmsh): ";
            for (int j = 0; j < determinants.size(); ++j)
                std::cout << determinants[j] << " ";
            std::cout << '\n';

            std::cout << "\tjacobians (gmsh): ";
            for (int j = 0; j < jacobians.size(); ++j)
                std::cout << jacobians[j] << " ";
            std::cout << '\n';

            // derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            int numComponents, numOrientations;
            std::vector<double> basisFunctions;
            gmsh::model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);
            std::cout << "\tbasis functions (gmsh): ";
            for (int j = 0; j < basisFunctions.size(); ++j)
                std::cout << basisFunctions[j] << " ";
            std::cout << '\n';

            /*---------computing the B matrix of the first element's (i.e elementTags[i][0]) first GP:------*/
            // ne serait-ça pas plutôt utile de directement convertir les Jacobiens en 2D ? travailler en full 2D
            // et faire pareil pour l'évaluation des gradients des basis functions ci-dessus ?
            std::vector<double> b_matrix;

            // convert first jacobian to Eigen format
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
                jacob(&jacobians[0], 3, 3);
            // passage en 2D xD
            Eigen::Matrix<double, 2, 2> jacob2D;
            jacob2D(0,0) = jacob(0,0);
            jacob2D(1,0) = jacob(1,0);
            jacob2D(0,1) = jacob(0,1);
            jacob2D(1,1) = jacob(1,1);
            std::cout << "\tfirst jacobian (Eigen):\n" 
                      << jacob2D.format(fmt) << '\n';
                    
            // compute the inverse with Eigen
            Eigen::Matrix<double, 2, 2> jacobinv = jacob2D.inverse();
            std::cout << "\tinverse of the first jacobian (Eigen):\n" 
                      << jacobinv.format(fmt) << "\n";

            // compute the transpose of the inverse with Eigen
            Eigen::Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();
            std::cout << "\ttranspose of the inverse of the first jacobian (Eigen):\n" 
                      << jacobinvtrans.format(fmt) << "\n";
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B_matrix(3,2*numNodes);
            B_matrix(0,0) = 0.;
            std::cout << "\ttranspose of the inverse of the first jacobian (Eigen):\n" 
                      << B_matrix.size() << "\n";
            for (int j = 0; j < weights.size(); ++j){ // looping over the number of shape functions
                
            }
        }
    }
    
    gmsh::finalize(); // clean-up things before ending the program
    return 0;
}
