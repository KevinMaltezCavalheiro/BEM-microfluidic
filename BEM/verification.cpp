#include "functionsBEM.hpp"

// Displays the time taken between the values of "start" and "end" (in seconds) if "showTime" is set to 'true'. "functionName" is
// the name of the function that has been timed.
void displayTime(const double start,
                 const double end,
                 std::string functionName, const bool showTime)
{
    if(showTime)
    {
        double duration = end - start;
        std::cout << functionName << "\n";
        std::cout << "Elapsed time: " << 1000000*duration << " [µs]\tor " << duration << " [s]\n\n";
    }
}

// Displays the different errors induced by the boundary element method compared with a simple analytical solution on a certain
// model geometry ("BEM_convergence.geo"). "domainEntityTags" is a vector containing the tag of each entity of the domain.
// "domainElementTypes" stores the types of the elements of the the domain. The first dimension of this vector represents the
// different entities and the second, the different element types. "domainElementTags" stores the tags of the elements of the
// domain. The first dimension of this vector represents the different entities, the second, the different element types and the
// third, the different element tags. "elementalPhi" and "elementalElectricField" are the used data for comparing numerical and
// analytical solutions. The convergence function can only be used with the "ElementData" data type.
void convergenceFun(const std::vector<int> &domainEntityTags,
                    const std::vector<std::vector<int>> &domainElementTypes,
                    const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> &elementalPhi,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> &elementalElectricField)
{
    double l = 1.0;
    double phiLeft = 100.0;
    double phiRight = 200.0;
    double exactE = (phiLeft - phiRight) / l;

    double xLeft = 0.25;
    double xRight = 0.75;
    double yBottom = 0.25;
    double yTop = 0.75;

    double errorPhi = 0.0;
    double errorE = 0.0;
    double errorEWindow = 0.0;
    int count = 0;
    int countWindow = 0;

    std::vector<double> barycenters;

    for(std::size_t i = 0; i < domainEntityTags.size(); i++)
    {
        for(std::size_t j = 0; j < domainElementTypes[i].size(); j++)
        {
            gmsh::model::mesh::getBarycenters(domainElementTypes[i][j], domainEntityTags[i], false, false, barycenters);
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                double numericalPhi = elementalPhi[i][j][k][0];
                double numericalE = elementalElectricField[i][j][k][0];

                double x = barycenters[3*k];
                double y = barycenters[3*k + 1];

                double exactPhi = phiLeft + (phiRight - phiLeft) * x / l;

                double relErrorPhi = abs( (exactPhi - numericalPhi) / exactPhi );
                double relErrorE = abs( (exactE - numericalE) / exactE);

                if(x > xLeft && x < xRight && y > yBottom && y < yTop)
                {
                    errorEWindow += relErrorE*relErrorE;
                    countWindow++;
                }

                errorPhi += relErrorPhi*relErrorPhi;
                errorE += relErrorE*relErrorE;
                count++;
            }
        }
    }
    errorPhi = sqrt(errorPhi / count);
    errorE = sqrt(errorE / count);
    errorEWindow = sqrt(errorEWindow / countWindow);

    std::cout << "\tRelative error on the electric potential: " << errorPhi << ".\n";
    std::cout << "\tRelative error on the electric field: " << errorE << ".\n";
    std::cout << "\tRelative error on the electric field (without the edges of the domain): " << errorEWindow << ".\n\n";
}