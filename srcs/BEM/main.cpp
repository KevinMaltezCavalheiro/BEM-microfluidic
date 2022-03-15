#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <math.h>
#include <Eigen/Dense>

// Each element has all its information gathered as in this structure.
struct elementStruct{
    int tag;            // tag of the element
    std::string bcType; // "dirichlet" or "neumann"
    double bcValue [2]; // values of [u, u_n]
    int nodeTags [2];   // tags of the [first, second] nodes of the element
};

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // Parameters for displaying variables (for debugging).
    int dispNodeCoord = 0;  // Coordinates of the nodes.
    int dispHierarchy = 0;  // Hierarchy (physical groups -> entities -> elements -> nodes).
    int dispAngles = 0;     // Angles used for the computation of h. /!\ Displays ~n^2 angles.
    int dispMatrices = 0;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
    int dispValues = 0;     // Elements and their values of u and u_n.
    int openGUI = 1;        // Graphical interface of GMSH.

    gmsh::initialize();

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);
    
    // Gets the dim and the tags of the physical groups
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags, 1);

    // Gets the boundary conditions from onelab and create a map: "name of the physical group" -> ("dirichlet"/"neumann", value)
    std::map<std::string, std::pair<std::string, double>> bcMap;

    std::vector<std::string> bcKeys;
    gmsh::onelab::getNames(bcKeys, "Boundary Conditions.+");
    for(auto &key : bcKeys)
    {
        std::vector<double> values;
        gmsh::onelab::getNumber(key, values);

        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/'))
            words.push_back(word);

        if(words.size() == 3)
            bcMap[words[1]] = {words[2], values[0]};
    }

    // Initialises a vector of elements (that will contain all the elements of interrest)
    std::vector<elementStruct> elementVector;
    
    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    int dirichletCount = 0;

    // Gets the tags, the coordinates and the parametric coordinates of the nodes.
    std::vector<std::size_t> allNodeTags;
    std::vector<double> allCoord;
    std::vector<double> allParametricCoord;
    gmsh::model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);
    
    // Print the coordinates of all nodes if dispNodeCoord != 0
    if(dispNodeCoord)
    {
        for(int i = 0; i < allNodeTags.size(); i++)
            std::cout << "Node '" << allNodeTags[i] << "': (" << allCoord[3*i] << "; " << allCoord[3*i + 1] << "; " << allCoord[3*i + 2] << ")\n";
        std::cout << "\n";
    }

    for (int i = 0; i < physicalGroupsDimTags.size(); ++i)
    {
        int physicalGroupsDim = physicalGroupsDimTags[i].first;
        int physicalGroupsTag = physicalGroupsDimTags[i].second;
        std::string name;
        gmsh::model::getPhysicalName(physicalGroupsDim, physicalGroupsTag, name);

        // Retrieves the type and the value of the boundary condition of this physical group.
        // /!\ Should check that this physical groups was well defined in the list of "Boundary Condition" !
        std::string bcType = bcMap[name].first;
        double bcValue = bcMap[name].second;
        
        // Prints the physical groups name, dim, tag, bcType and bcValue (if dispHierarchy != 0).
        if(dispHierarchy)
            std::cout << "Physical group '" << name << "' (" << physicalGroupsDim << "D) with tag '" << physicalGroupsTag << "' has a '" << bcType << "' type and the value " << bcValue << "\n";

        std::vector<int> entitiesTags;
        gmsh::model::getEntitiesForPhysicalGroup(physicalGroupsDim, physicalGroupsTag, entitiesTags);

        for (int k = 0; k < entitiesTags.size(); k++)
        {
            // Print the tag of the entities of each physical group (if dispHierarchy != 0).
            if(dispHierarchy)
                std::cout << "\t contains entity '" << entitiesTags[k] << "' ";

            std::vector<std::size_t> elementTags;
            std::vector<std::size_t> nodeTags;
            gmsh::model::mesh::getElementsByType(1, elementTags, nodeTags, entitiesTags[k]);
            
            // Print the number of element of the given entity (if dispHierarchy != 0).
            if(dispHierarchy)
                std::cout << "which has " << elementTags.size() << " elements:\n";

            for(int el = 0; el < elementTags.size(); el++)
            {
                // Stores the boundary type and the index of the element.
                int elTag = elementTags[el];
                
                elementStruct element;
                element.bcType = bcType;
                element.tag = elTag;
                element.nodeTags[0] = nodeTags[2*el];
                element.nodeTags[1] = nodeTags[2*el + 1];
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
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*(element.nodeTags[0]-1)] << "; " << allCoord[3*(element.nodeTags[0]-1) + 1] << ")\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*(element.nodeTags[1]-1)] << "; " << allCoord[3*(element.nodeTags[1]-1) + 1] << ")\n";
                }
            }
            if(dispHierarchy)
                std::cout << "\n";
        }
    }

    // Declaration of matrices A, B and vector c;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(elementVector.size(), 1);
    
    //  Double loop over each element in order to compute the element h and g that will be given to the matrices A and B.
    for(int i = 0; i < elementVector.size(); i++)
    {
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
                double jMidX = (j1X + j2X) / 2;
                double jMidY = (j1Y + j2Y) / 2;
                double jDeltaX = j2X - j1X;
                double jDeltaY = j2Y - j1Y;
                double jLength = sqrt(pow(jDeltaX, 2.0) + pow(jDeltaY, 2.0));

                // Compute the angles between mid point of element 'i' and nodes 1/2 of element 'j' (in radian).
                double a1 = atan2((j1Y - iMidY),(j1X - iMidX));
                double a2 = atan2((j2Y - iMidY),(j2X - iMidX));
                while(a2 < a1)  // This is not general and should be reviewed !
                    a2 += 2*M_PI;

                // Angle verification
                if(dispAngles)
                {
                    double angleConst = 180 / M_PI;
                    double a1_deg = angleConst * a1;
                    double a2_deg = angleConst * a2;
                    std::cout << "Angles from element with tag '" << elementVector[i].tag << "' to element with tag '" << elementVector[j].tag << "':\n";
                    std::cout << "\tNode: tag '" << j1Tag << "': " << a1_deg << "\n";
                    std::cout << "\tNode: tag '" << j2Tag << "': " << a2_deg << "\n";
                }
                
                h = (a2 - a1) / (2*M_PI);
                
                // Computation of the integration points (parametric) and weights.
                std::string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.
                std::vector<double> localCoord;
                std::vector<double> weights;
                gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

                // Computation of the sum (of integrand evaluated at integration point * weight)
                double gSum = 0;
                for(int k = 0; k < weights.size(); k++)
                {
                    // Parametric coordinate
                    double xi = localCoord[3*k];
                    // Global coordinates
                    double x = jMidX + jDeltaX*xi/2;
                    double y = jMidY + jDeltaY*xi/2;

                    // Distance between mid point of element 'i' and integration point.
                    double r = sqrt(pow(iMidX - x, 2.0) + pow(iMidY - y, 2.0));
                    gSum += log(r) * weights[k];
                }
                
                g = gSum * jLength / (4*M_PI);
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
    Eigen::MatrixXd b = B*c;
    //Eigen::MatrixXd x = A.householderQr().solve(b);
    Eigen::MatrixXd x = A.fullPivLu().solve(b);
    
    // Matrices verification
    if(dispMatrices)
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

    // Assign each value of x to the right element.
    for(int i = 0; i < elementVector.size(); i++)
    {
        if(i < dirichletCount)
            elementVector[i].bcValue[1] = x(i);
        else
            elementVector[i].bcValue[0] = x(i);
        
        if(dispValues)
            std::cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
    }

    /*
     *  Visualisation
     */

    std::vector<std::string> names;
    gmsh::model::list(names);
    std::cout << "the file contains " << names.size() << " model names\n";

    int viewNodeTag = gmsh::view::add("NodeData");
    std::vector<std::vector<double>> data(allNodeTags.size());

    for(int i = 0; i < allNodeTags.size(); i++)
    {
        // Coordinates of the node
        double x = allCoord[3*i];
        double y = allCoord[3*i + 1];

        double phi = 0;
        for(int j = 0; j < elementVector.size(); j++)
        {
            double h;
            double g;
            
            // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element i here).
            int j1Tag = elementVector[j].nodeTags[0];
            int j2Tag = elementVector[j].nodeTags[1];
            if(j1Tag == i+1 || j2Tag == i+1)
            {
                phi = elementVector[j].bcValue[0];
                break;
            }

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

            double a1 = atan2((j1Y - y),(j1X - x));
            double a2 = atan2((j2Y - y),(j2X - x));
            while(a2 < a1)
                a2 += 2*M_PI;
            h = (a2 - a1) / (2*M_PI);
            
            // Computation of the integration points (parametric) and weights.
            std::string integrationType = "Gauss4"; // Could vary depending on the distance between the elements.
            std::vector<double> localCoord;
            std::vector<double> weights;
            gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

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
        data[i].resize(1);
        data[i][0] = phi;
    }
    gmsh::view::addModelData(viewNodeTag, 0, names[0], "NodeData", allNodeTags, data, 0.0, 1);
    gmsh::view::write(viewNodeTag, "results.msh");

    gmsh::option::setNumber("View[0].IntervalsType", 3);
    gmsh::option::setNumber("View[0].AdaptVisualizationGrid", 1);
    gmsh::option::setNumber("View[0].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[0].TargetError", -0.0001);

    if(openGUI)
        gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}