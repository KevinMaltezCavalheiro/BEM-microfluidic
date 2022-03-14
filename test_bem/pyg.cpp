#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <math.h>
#include <Eigen/Dense>

// test for commit
// Each element has all its information gathered as in this structure.
struct elementStruct{
    std::string bcType; // "dirichlet" or "neumann"
    double bcValue [2]; // value of u/u_n
    int globalIndex;    // each element has a tag, but also an index that represent an order (from 0 to n-1)
    int typeIndex;  // index of the "dirichlet"/"neumann" element
    int nodeTags [2];   // tags of the first/second nodes of the element
};

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // Parameters for displaying variables (for debugging).
    int dispNodeCoord = 1;  // Coordinates of the nodes.
    int dispHierarchy = 1;  // Hierarchy (physical groups -> entities -> elements -> nodes).
    int dispAngles = 1;     // Angles used for the computation of h. /!\ Displays ~n^2 angles.
    int dispMatrices = 1;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
    int dispValues = 1;     // Elements and their values of u and u_n.
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

    // Initialises a map: elementTag -> elementStruct
    std::map<int, elementStruct> elementMap;
    
    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    int dirichletIndex = 0;
    int neumannIndex = 0;

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
            std::cout << "Physical group '" << name << "' (" << physicalGroupsDim << "D) with tag '" << physicalGroupsTag << "' is a '" << bcType << "' element and has the value " << bcValue << "\n";

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
                element.nodeTags[0] = nodeTags[2*el];
                element.nodeTags[1] = nodeTags[2*el + 1];
                element.globalIndex = dirichletIndex + neumannIndex;
                if(bcType == "dirichlet")
                {
                    element.bcValue[0] = bcValue;
                    element.typeIndex = dirichletIndex;
                    dirichletIndex++;
                }
                else if(bcType == "neumann")
                {
                    element.bcValue[1] = bcValue;
                    element.typeIndex = neumannIndex;
                    neumannIndex++;
                }
                elementMap[elTag] = element;

                
                // Print the tag of a given element and the tags of the nodes that composes it as well as their coordinates (if dispHierarchy != 0).
                if(dispHierarchy)
                {
                    std::cout << "\t\t-Element '" << elTag << "' (" << elementMap[elTag].bcType << ") has index '" << elementMap[elTag].globalIndex << "':\n";
                    std::cout << "\t\t\t-> Node '" << elementMap[elTag].nodeTags[0] << "' (" << allCoord[3*(elementMap[elTag].nodeTags[0]-1)] << "; " << allCoord[3*(elementMap[elTag].nodeTags[0]-1) + 1] << ")\n";
                    std::cout << "\t\t\t-> Node '" << elementMap[elTag].nodeTags[1] << "' (" << allCoord[3*(elementMap[elTag].nodeTags[1]-1)] << "; " << allCoord[3*(elementMap[elTag].nodeTags[1]-1) + 1] << ")\n";
                }
            }
            if(dispHierarchy)
                std::cout << "\n";
        }
    }

    // Declaration of matrices A, B and vector c;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(elementMap.size(), elementMap.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(elementMap.size(), elementMap.size());
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(elementMap.size(), 1);
    
    //  Double loop over each element (accessed by their tag i/j) in order to compute the element h and g that will be given to the matrices A and B.
    //  /!\ The tags start at 1, so the loops are from i = 1 to i = elementMap.size()
    for(int i = 1; i <= elementMap.size(); i++)
    {
            // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element i here).
            int i1Tag = elementMap[i].nodeTags[0];
            int i2Tag = elementMap[i].nodeTags[1];

            // Get the coordinates of the nodes from their tag.
            double i1X = allCoord[3*(i1Tag-1)];
            double i1Y = allCoord[3*(i1Tag-1) + 1];
            double i2X = allCoord[3*(i2Tag-1)];
            double i2Y = allCoord[3*(i2Tag-1) + 1];

            // Compute the coordinates of the mid point of element i and its length.
            double iMidX = (i1X + i2X) / 2;
            double iMidY = (i1Y + i2Y) / 2;
            double iLength = sqrt(pow(i1X - i2X, 2.0) + pow(i1Y - i2Y, 2.0));

        for(int j = 1; j <= elementMap.size(); j++)
        {
            // Will be computed and stored in A and B.
            double h;
            double g;

            if(i == j)  // Formula for diagonal elements.
            {
                h = -0.5;
                g = iLength * (log(iLength/2)-1) / (2*M_PI);
            }
            else    // General case
            {
                // Indices '1'/'2' stand for 'Node 1'/'Node 2' (of element j here).
                int j1Tag = elementMap[j].nodeTags[0];
                int j2Tag = elementMap[j].nodeTags[1];

                // Get the coordinates of the nodes from their tag.
                double j1X = allCoord[3*(j1Tag-1)];
                double j1Y = allCoord[3*(j1Tag-1) + 1];
                double j2X = allCoord[3*(j2Tag-1)];
                double j2Y = allCoord[3*(j2Tag-1) + 1];

                // Compute the coordinates of the mid point of element j, deltaX and deltaY and its length.
                double jMidX = (j1X + j2X) / 2;
                double jMidY = (j1Y + j2Y) / 2;
                double jDeltaX = j2X - j1X;
                double jDeltaY = j2Y - j1Y;
                double jLength = sqrt(pow(jDeltaX, 2.0) + pow(jDeltaY, 2.0));

                // Compute the angles between mid point of element i and nodes 1/2 of element j (in radian).
                double a1 = atan2((j1Y - iMidY),(j1X - iMidX));
                double a2 = atan2((j2Y - iMidY),(j2X - iMidX));

                // Angle verification  
                if(dispAngles)
                {
                    double angleConst = 180 / M_PI;
                    double a1_deg = angleConst * a1;
                    double a2_deg = angleConst * a2;
                    std::cout << "Angles from element with tag '" << i << "' to element with tag '" << j << "':\n";
                    std::cout << "\tNode: tag '" << j1Tag << "': " << a1_deg << "\n";
                    std::cout << "\tNode: tag '" << j2Tag << "': " << a2_deg << "\n";
                }
                
                h = (a2 - a1) / (2*M_PI);

                double sum = 0;
                
                // Computation of the integration points (parametric) and weights.
                std::string integrationType = "Gauss2"; // Could vary depending on the distance between the elements.
                std::vector<double> localCoord;
                std::vector<double> weights;
                gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

                // Computation of the sum (of integrand evaluated at integration point * weight)
                for(int k = 0; k < weights.size(); k++)
                {
                    double xi = localCoord[3*k];    // Parametric coordinate
                    // Global coordinates
                    double x = jMidX + jDeltaX*xi/2;    
                    double y = jMidY + jDeltaY*xi/2;

                    // Distance between mid point of element i and integration point.
                    double r = sqrt(pow(iMidX - x, 2.0) + pow(iMidY - y, 2.0));
                    sum += log(r) * weights[k];
                }
                
                g = sum * jLength / (4*M_PI);
            }
            
            // Storage of h and g at the right place (in matrices A and B). The method is described below.
            
            // Store h at the right place.
            // the i-component is the globalIndex of element i.
            // if element j is dirichlet --> in matrix B (with a minus sign).
                // the j-component is the dirichletIndex of element j.
            // if element j is neumann --> in matrix A
                // the j-component is the neumannIndex of element j.

            // Store g at the right place.
            // the i-component is the globalIndex of element i.
            // if element j is dirichlet --> in matrix A (with a minus sign).
                // the j-component is the dirichletIndex of element j + the number of neumann elements.
            // if element j is neumann --> in matrix B.
                // the j-component is the neumannIndex of element j + the number of dirichlet elements.

            if(elementMap[j].bcType == "dirichlet")
            {
                B(elementMap[i].globalIndex, elementMap[j].typeIndex) = -h;
                A(elementMap[i].globalIndex, elementMap[j].typeIndex + neumannIndex) = -g;
                c(elementMap[j].typeIndex) = elementMap[j].bcValue[0];
            }
            else if(elementMap[j].bcType == "neumann")
            {
                A(elementMap[i].globalIndex, elementMap[j].typeIndex) = h;
                B(elementMap[i].globalIndex, elementMap[j].typeIndex + dirichletIndex) = g;
                c(elementMap[j].typeIndex + dirichletIndex) = elementMap[j].bcValue[1];
            }
        }
    }

    // Compute b then x.
    Eigen::MatrixXd b = B*c;
    Eigen::MatrixXd x = A.householderQr().solve(b);
    
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
    // This is done by doing a loop over the element tags (from 1 to elementMap.size() !) and seeking which component of x is interresting for element i.
    // Then, thanks to the type of element i, one can determine is the value of x represent u or u_n.
    for(int i = 1; i <= elementMap.size(); i++)
    {
        if(elementMap[i].bcType == "dirichlet")
            elementMap[i].bcValue[1] = x(elementMap[i].typeIndex + neumannIndex);
        else if(elementMap[i].bcType == "neumann")
            elementMap[i].bcValue[0] = x(elementMap[i].typeIndex);
        
        if(dispValues)
            std::cout << "Element '" << i << "': from node " << elementMap[i].nodeTags[0] << " to node " << elementMap[i].nodeTags[1] << " (" << elementMap[i].bcType << "):\tDirichlet: " << elementMap[i].bcValue[0] << "\tNeumann: " << elementMap[i].bcValue[1] << "\n";
    }

    if(openGUI)
        gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}
