#include <gmsh.h>
#include <iostream>
#include <map>
#include <sstream>
#include <math.h>
#include <Eigen/Dense>

// Each element has all its information gathered as in this structure.
struct elementStruct{
    int tag;            // tag of the element
    std::string bcType; // "dirichlet" or "neumann" (at the end, should replace this by an int or boolean)
    double bcValue [2]; // values of [u, u_n]
    int nodeTags [2];   // tags of the [first, second] nodes of the element

    double length;
    double x1, y1, x2, y2;
    double midX, midY;
    double deltaX, deltaY;
};

void dispNodeCoordFun(const std::vector<size_t> &allNodeTags, const std::vector<double> &allCoord);
void dispMatricesFun(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const Eigen::MatrixXd &c, const Eigen::MatrixXd &b, const Eigen::MatrixXd &x);
void dispValuesFun(const std::vector<elementStruct> &elementVector);

void getBC(const std::vector<double> &allCoord, std::vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy);
void computeGeometricParam(elementStruct &element, const std::vector<double> &allCoord);
double computeH(const double x, const double y, const elementStruct &element);
double computeG(const double x, const double y, const elementStruct &element, const std::string integrationType);

void computeData(const std::vector<elementStruct> &elementVector, const std::vector<int> &elementTypes, const std::vector<std::vector<size_t>> &elementTags, std::vector<size_t> &tags, std::vector<std::vector<double>> &data);
void displayData(const std::string viewName, const std::string modelName, const std::string fileName, const std::vector<size_t> &tags, const std::vector<std::vector<double>> &data);

void dispNodeCoordFun(const std::vector<size_t> &allNodeTags, const std::vector<double> &allCoord)
{
    for(int i = 0; i < allNodeTags.size(); i++)
        std::cout << "Node '" << allNodeTags[i] << "': (" << allCoord[3*i] << "; " << allCoord[3*i + 1] << "; " << allCoord[3*i + 2] << ")\n";
    std::cout << "\n";
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
    for(int i = 0; i < elementVector.size(); i++)
    {
        std::cout << "Element '" << elementVector[i].tag << "': from node " << elementVector[i].nodeTags[0] << " to node " << elementVector[i].nodeTags[1] << " (" << elementVector[i].bcType << "):\tDirichlet: " << elementVector[i].bcValue[0] << "\tNeumann: " << elementVector[i].bcValue[1] << "\n";
    }
}


void getBC(const std::vector<double> &allCoord, std::vector<elementStruct> &elementVector, int &dirichletCount, const bool dispHierarchy)
{
    // Gets the dim and the tags of the physical groups
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

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

    // Initially, these are 0 and are incremented each time a dirichlet/neumann-element is encountered.
    for (int i = 0; i < physicalGroupsDimTags.size(); i++)
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
                // Stores the tag, the boundary type, the boundary value and the node tags of the element.
                int elTag = elementTags[el];
                
                elementStruct element;
                element.tag = elTag;
                element.bcType = bcType;
                element.nodeTags[0] = nodeTags[2*el];
                element.nodeTags[1] = nodeTags[2*el + 1];
                computeGeometricParam(element, allCoord);
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
                    std::cout << "\t\t-Element '" << elTag << "' (" << element.bcType << "):\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[0] << "' (" << allCoord[3*(element.nodeTags[0]-1)] << "; " << allCoord[3*(element.nodeTags[0]-1) + 1] << ")\n";
                    std::cout << "\t\t\t-> Node '" << element.nodeTags[1] << "' (" << allCoord[3*(element.nodeTags[1]-1)] << "; " << allCoord[3*(element.nodeTags[1]-1) + 1] << ")\n";
                }
            }
            if(dispHierarchy)
                std::cout << "\n";
        }
    }
}


void computeGeometricParam(elementStruct &element, const std::vector<double> &allCoord)
{
    int i1Tag = element.nodeTags[0];
    int i2Tag = element.nodeTags[1];

    // Get the coordinates of the nodes from their tag.
    element.x1 = allCoord[3*(i1Tag-1)];
    element.y1 = allCoord[3*(i1Tag-1) + 1];
    element.x2 = allCoord[3*(i2Tag-1)];
    element.y2 = allCoord[3*(i2Tag-1) + 1];

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
    
    double cosAlpha = (deltaX1*deltaX2 + deltaY1*deltaY2) / sqrt((pow(deltaX1,2.0) + pow(deltaY1,2.0))*(pow(deltaX2,2.0) + pow(deltaY2,2.0)));
    if(cosAlpha < -1.0)
        cosAlpha = -1.0;
    else if(cosAlpha > 1.0)
        cosAlpha = 1.0;

    double alpha = acos(cosAlpha);

    double crossProduct = deltaX1*deltaY2 - deltaY1*deltaX2;
    if(crossProduct < 0)
        alpha = -alpha;
    
    return alpha / (2*M_PI);
}


double computeG(const double x, const double y, const elementStruct &element, const std::string integrationType)
{
    // Computation of the integration points (parametric) and weights.
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
        double gX = element.midX + element.deltaX*xi/2;
        double gY = element.midY + element.deltaY*xi/2;

        // Distance between mid point of element 'i' and integration point.
        double r = sqrt(pow(x - gX, 2.0) + pow(y - gY, 2.0));
        gSum += log(r) * weights[k];
    }
    
    return gSum * element.length / (4*M_PI);
}


void computeData(const std::vector<elementStruct> &elementVector, const std::vector<int> &elementTypes, const std::vector<std::vector<size_t>> &elementTags, std::vector<size_t> &tags, std::vector<std::vector<double>> &data)
{
    int elementCount = 0;
    for(int i = 0; i < elementTypes.size(); i++)
    {
        int elType = elementTypes[i];
        std::vector<double> barycenters;
        gmsh::model::mesh::getBarycenters(elType, -1, false, false, barycenters);
        
        for(int j = 0; j < elementTags[i].size(); j++)
        {
            double x = barycenters[3*j];
            double y = barycenters[3*j + 1];

            double phi = 0;
            for(int k = 0; k < elementVector.size(); k++)
            {
                double h = computeH(x, y, elementVector[k]);
                double g = computeG(x, y, elementVector[k], "Gauss4");
                phi += h*elementVector[k].bcValue[0] - g*elementVector[k].bcValue[1];
            }
            // Potential at the node is phi.
            tags[elementCount] = elementTags[i][j];
            data[elementCount] = std::vector<double> (1);
            data[elementCount][0] = phi;
            elementCount++;
        }
    }
}

void displayData(const std::string viewName, const std::string modelName, const std::string fileName, const std::vector<size_t> &tags, const std::vector<std::vector<double>> &data)
{
    int viewTag = gmsh::view::add(viewName);
    gmsh::view::addModelData(viewTag, 0, modelName, "ElementData", tags, data);
    gmsh::view::write(viewTag, fileName);

    gmsh::option::setNumber("View[0].IntervalsType", 3);
    // gmsh::option::setNumber("View[0].AdaptVisualizationGrid", 1); // Problem with this command !
    gmsh::option::setNumber("View[0].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[0].TargetError", -0.0001);
}


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    // Parameters for displaying variables (for debugging).
    bool dispNodeCoord = false;  // Coordinates of the nodes.
    bool dispHierarchy = false;  // Hierarchy (physical groups -> entities -> elements -> nodes).
    bool dispMatrices = false;   // Matrices A, B, c, b, x. /!\ May lead to visualisation problems for large matrices.
    bool dispValues = false;     // Elements and their values of u and u_n.
    bool openGUI = true;        // Graphical interface of GMSH.

    gmsh::initialize();

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    std::vector<elementStruct> elementVector;

    // Gets the tags, the coordinates and the parametric coordinates of the nodes.
    std::vector<std::size_t> allNodeTags;
    std::vector<double> allCoord;
    std::vector<double> allParametricCoord;
    gmsh::model::mesh::getNodes(allNodeTags, allCoord, allParametricCoord);
    
    // Print the coordinates of all nodes if dispNodeCoord != 0
    if(dispNodeCoord)
        dispNodeCoordFun(allNodeTags, allCoord);

    int dirichletCount = 0;
    getBC(allCoord, elementVector, dirichletCount, dispHierarchy);

    // Declaration of matrices A, B and vector c (the matrix type is perhaps not the most efficient).
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(elementVector.size(), 1);

    //  Double loop over each element in order to compute the element h and g that will be given to the matrices A and B.
    for(int i = 0; i < elementVector.size(); i++)
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
                g = computeG(elementVector[i].midX, elementVector[i].midY, elementVector[j], "Gauss4");
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

    Eigen::MatrixXd b = B*c;
    Eigen::MatrixXd x = A.fullPivLu().solve(b);
    
    // Matrices verification
    if(dispMatrices)
        dispMatricesFun(A, B, c, b, x);

    // Assign each value of x to the right element.
    for(int i = 0; i < elementVector.size(); i++)
    {
        if(i < dirichletCount)
            elementVector[i].bcValue[1] = x(i);
        else
            elementVector[i].bcValue[0] = x(i);
    }
    if(dispValues)
        dispValuesFun(elementVector);

    /*
     *  Visualisation
     */
    std::vector<std::string> names;
    gmsh::model::list(names);
    
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);
    
    int nb2dElements = 0;
    for(int i = 0; i < elementTypes.size(); i++)
        nb2dElements += elementTags[i].size();

    std::vector<std::vector<double>> phi(nb2dElements);
    std::vector<std::size_t> tags(nb2dElements);

    computeData(elementVector, elementTypes, elementTags, tags, phi);
    displayData("phi", names[0], "results.msh", tags, phi);

    if(openGUI)
        gmsh::fltk::run();
    
    gmsh::finalize();
    return 0;
}