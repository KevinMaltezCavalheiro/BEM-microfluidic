#include "functionsBEM.hpp"

// Computes the 'h' element analytically. "x" and "y" are the coordinates at which 'h' should be computed. "element" is the
// element 'j' while will helps to compute 'h'.
double computeAnalyticalH(const double x, const double y,
                           const elementStruct &element, const bool ccw)
{
    double deltaX1 = element.x1 - x;
    double deltaX2 = element.x2 - x;
    double deltaY1 = element.y1 - y;
    double deltaY2 = element.y2 - y;

    double cosAlpha = (deltaX1*deltaX2 + deltaY1*deltaY2) /
                      (sqrt((pow(deltaX1,2.0) + pow(deltaY1,2.0)) *
                             (pow(deltaX2,2.0) + pow(deltaY2,2.0))));
    if(cosAlpha < -1.0) cosAlpha = -1.0;
    else if(cosAlpha > 1.0) cosAlpha = 1.0;

    double alpha = acos(cosAlpha);
    double crossProduct = deltaX1*deltaY2 - deltaY1*deltaX2;
    if(crossProduct < 0) alpha = -alpha;

    if(!ccw) alpha = -alpha;

    return alpha / (2*M_PI);
}

// Computes the 'g' element analytically. "x" and "y" are the coordinates at which 'g' should be computed. "element" is the
// element 'j' while will helps to compute 'g'.
double computeAnalyticalG(const double x,
                          const double y,
                          const elementStruct &element)
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
// Computes the 'hx' component of the gradient analytically. "x" and "y" are the coordinates at which 'gradHx' should be
// computed. "element" is the element 'j' while will helps to compute 'gradHx'.
double computeAnalyticalGradHx(const double x,
                               const double y,
                               const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;
    double nx = element.nx;
    double ny = element.ny;


    double gradHx_analytique = (lj*(nx*(pow(xi,2) + xj*xjp1 - xi*(xj + xjp1) - (yi - yj)*(yi - yjp1)) -
                                    ny*(xj*yi + xjp1*yi - xjp1*yj - xj*yjp1 + xi*(-2*yi + yj + yjp1))))/
                                (2.*Pi*(pow(xi - xj,2) + pow(yi - yj,2))*
                                  (pow(xi - xjp1,2) + pow(yi - yjp1,2)));

    return gradHx_analytique;
}
// Computes the 'hy' component of the gradient analytically. "x" and "y" are the coordinates at which 'gradHy' should be
// computed. "element" is the element 'j' while will helps to compute 'gradHy'.
double computeAnalyticalGradHy(const double x,
                               const double y,
                               const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;
    double nx = element.nx;
    double ny = element.ny;

    double gradHy_analytique = -(lj*(ny*(pow(xi,2) + xj*xjp1 - xi*(xj + xjp1) - (yi - yj)*(yi - yjp1)) +
                                     nx*(xj*yi + xjp1*yi - xjp1*yj - xj*yjp1 + xi*(-2*yi + yj + yjp1))))/
                                (2.*Pi*(pow(xi - xj,2) + pow(yi - yj,2))*
                                 (pow(xi - xjp1,2) + pow(yi - yjp1,2)));

    return gradHy_analytique;
}
// Computes the 'gx' component of the gradient analytically. "x" and "y" are the coordinates at which 'gradGx' should be
// computed. "element" is the element 'j' while will helps to compute 'gradGx'.
double computeAnalyticalGradGx(const double x,
                               const double y,
                               const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double gradGx_analytique = (lj*(-2*(yj - yjp1)*(atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1))) +
              atan(((xi - xjp1)*(-xj + xjp1) + (yi - yjp1)*(-yj + yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1)))) -
         (xj - xjp1)*(log(pow(xi - xj,2) + pow(yi - yj,2)) - log(pow(xi - xjp1,2) + pow(yi - yjp1,2)))))/(4.*Pi*(pow(xj - xjp1,2) + pow(yj - yjp1,2)));

    return gradGx_analytique;
}
// Computes the 'gy' component of the gradient analytically. "x" and "y" are the coordinates at which 'gradGy' should be
// computed. "element" is the element 'j' while will helps to compute 'gradGy'.
double computeAnalyticalGradGy(const double x,
                               const double y,
                               const elementStruct &element)
{
    double lj = element.length;
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;
    double Pi = M_PI;

    double gradGy_analytique = (lj*(2*(xj - xjp1)*(atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1))) +
              atan(((xi - xjp1)*(-xj + xjp1) + (yi - yjp1)*(-yj + yjp1))/(xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1)))) -
         (yj - yjp1)*(log(pow(xi - xj,2) + pow(yi - yj,2)) - log(pow(xi - xjp1,2) + pow(yi - yjp1,2)))))/(4.*Pi*(pow(xj - xjp1,2) + pow(yj - yjp1,2)));
    return gradGy_analytique;
}

// Computes the 'g' element numerically. "x" and "y" are the coordinates at which 'g' should be computed. "element" is the
// element 'j' while will helps to compute 'g'. "integrationType" is the type of numerical integration used (in case of numerical
// computation).
double computeNumericalG(const double x,
                         const double y,
                         const elementStruct &element,
                         const std::string integrationType)
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
// Computes the 'hx' component of the gradient numerically. "x" and "y" are the coordinates at which 'gradHx' should be computed.
// "element" is the element 'j' while will helps to compute 'gradHx'. "integrationType" is the type of numerical integration used
// (in case of numerical computation).
double computeNumericalGradHx(const double x,
                              const double y,
                              const elementStruct &element,
                              const std::string integrationType)
{
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;

    double nx = element.nx;
    double ny = element.ny;

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
// Computes the 'hy' component of the gradient numerically. "x" and "y" are the coordinates at which 'gradHy' should be computed.
// "element" is the element 'j' while will helps to compute 'gradHy'. "integrationType" is the type of numerical integration used
// (in case of numerical computation).
double computeNumericalGradHy(const double x,
                              const double y,
                              const elementStruct &element,
                              const std::string integrationType)
{
    double xi = x;
    double yi = y;
    double xj = element.x1;
    double xjp1 = element.x2;
    double yj = element.y1;
    double yjp1 = element.y2;

    double nx = element.nx;
    double ny = element.ny;

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
// Computes the 'gx' component of the gradient numerically. "x" and "y" are the coordinates at which 'gradGx' should be computed.
// "element" is the element 'j' while will helps to compute 'gradGx'. "integrationType" is the type of numerical integration used
// (in case of numerical computation).
double computeNumericalGradGx(const double x,
                              const double y,
                              const elementStruct &element,
                              const std::string integrationType)
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
// Computes the 'gy' component of the gradient numerically. "x" and "y" are the coordinates at which 'gradGy' should be computed.
// "element" is the element 'j' while will helps to compute 'gradGy'. "integrationType" is the type of numerical integration used
// (in case of numerical computation).
double computeNumericalGradGy(const double x,
                              const double y,
                              const elementStruct &element,
                              const std::string integrationType)
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
