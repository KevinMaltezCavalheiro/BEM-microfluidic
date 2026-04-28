#include "functionsBEM.hpp"

// Computes the 'h' element analytically. "x" and "y" are the coordinates at which 'h' should be computed. "element" is the
// element 'j' while will helps to compute 'h'.
double computeAnalyticalH(const double x, const double y,
                          const elementStruct &element, const bool ccw)
{
    const double dX1 = element.x1 - x;
    const double dX2 = element.x2 - x;
    const double dY1 = element.y1 - y;
    const double dY2 = element.y2 - y;

    // Produit scalaire et distances au carré (plus rapide que pow)
    const double dot = dX1 * dX2 + dY1 * dY2;
    const double distSq1 = dX1 * dX1 + dY1 * dY1;
    const double distSq2 = dX2 * dX2 + dY2 * dY2;

    // Utilisation de M_1_PI (1/pi) prédéfini pour éviter la division
    // sqrt(a*b) est identique à sqrt(a)*sqrt(b) mais plus rapide ici
    double cosAlpha = dot / std::sqrt(distSq1 * distSq2);

    // Clamp rapide
    if(cosAlpha < -1.0) cosAlpha = -1.0;
    else if(cosAlpha > 1.0) cosAlpha = 1.0;

    double alpha = std::acos(cosAlpha);

    // Produit vectoriel 2D pour le signe
    if((dX1 * dY2 - dY1 * dX2) < 0) alpha = -alpha;
    if(!ccw) alpha = -alpha;

    // alpha / (2 * PI) -> alpha * (0.5 * 1/PI)
    return alpha * (0.5 * M_1_PI);
}

// Computes the 'g' element analytically. "x" and "y" are the coordinates at which 'g' should be computed. "element" is the
// element 'j' while will helps to compute 'g'.
double computeAnalyticalG(const double x, const double y, const elementStruct &element)
{
    // On garde tes noms de variables exacts
    const double lj = element.length;
    const double xi = x;
    const double yi = y;
    const double xj = element.x1;
    const double xjp1 = element.x2;
    const double yj = element.y1;
    const double yjp1 = element.y2;
    // On remplace Pi par une constante pour éviter la division
    const double inv_4Pi = 0.25 * M_1_PI;

    // Étape 1 : Remplacer pow(..., 2) par (a*a) - C'est CA le gain de temps
    const double dx_seg = xj - xjp1;
    const double dy_seg = yj - yjp1;
    const double lj_sq = dx_seg * dx_seg + dy_seg * dy_seg;

    // Étape 2 : Pré-calculer les termes qui apparaissent plusieurs fois (sans changer leur place)
    const double common_atan_den = -(xj*yi) + xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*yjp1;
    const double rj_sq = (xi - xj)*(xi - xj) + (yi - yj)*(yi - yj);
    const double rjp1_sq = (xi - xjp1)*(xi - xjp1) + (yi - yjp1)*(yi - yjp1);

    // Étape 3 : Ta formule EXACTE avec juste les variables remplacées
    // J'ai recopié ton bloc sans changer l'ordre ni les signes
    double g_analytique = (lj * (
        -2.0 * lj_sq
        + 2.0 * (common_atan_den) * std::atan(((xi - xjp1)*(xj - xjp1) + (yi - yjp1)*(yj - yjp1)) / common_atan_den)
        + 2.0 * (-(xjp1*yi) - xi*yj + xjp1*yj + xj*(yi - yjp1) + xi*yjp1) * std::atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1)) / (xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1)))
        + (-((xi - xj)*(xj - xjp1)) - (yi - yj)*(yj - yjp1)) * std::log(rj_sq)
        + ((xi - xjp1)*(xj - xjp1) + (yi - yjp1)*(yj - yjp1)) * std::log(rjp1_sq)
    )) / (4.0 * M_PI * lj_sq);

    return g_analytique;
}

double computeAnalyticalGradHx(const double x, const double y, const elementStruct &element)
{
    const double xi = x, yi = y;
    const double xj = element.x1, yj = element.y1;
    const double xjp1 = element.x2, yjp1 = element.y2;

    const double dxj = xi - xj;
    const double dyj = yi - yj;
    const double dxjp1 = xi - xjp1;
    const double dyjp1 = yi - yjp1;

    // Remplacement de pow(..., 2) par multiplication
    const double dist_j_sq = dxj*dxj + dyj*dyj;
    const double dist_jp1_sq = dxjp1*dxjp1 + dyjp1*dyjp1;

    const double common_denom = 2.0 * M_PI * dist_j_sq * dist_jp1_sq;

    // Numérateur sans pow
    const double term1 = xi*xi + xj*xjp1 - xi*(xj + xjp1) - (yi - yj)*(yi - yjp1);
    const double term2 = xj*yi + xjp1*yi - xjp1*yj - xj*yjp1 + xi*(-2.0*yi + yj + yjp1);

    double gradHx = (element.length * (element.nx * term1 - element.ny * term2)) / common_denom;

    if(!element.normalOutward) gradHx = -gradHx;
    return gradHx;
}

double computeAnalyticalGradHy(const double x, const double y, const elementStruct &element)
{
    const double xi = x, yi = y;
    const double xj = element.x1, yj = element.y1;
    const double xjp1 = element.x2, yjp1 = element.y2;

    const double dist_j_sq = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);
    const double dist_jp1_sq = (xi-xjp1)*(xi-xjp1) + (yi-yjp1)*(yi-yjp1);

    const double common_denom = 2.0 * M_PI * dist_j_sq * dist_jp1_sq;

    const double term1 = xi*xi + xj*xjp1 - xi*(xj + xjp1) - (yi - yj)*(yi - yjp1);
    const double term2 = xj*yi + xjp1*yi - xjp1*yj - xj*yjp1 + xi*(-2.0*yi + yj + yjp1);

    double gradHy = -(element.length * (element.ny * term1 + element.nx * term2)) / common_denom;

    if(!element.normalOutward) gradHy = -gradHy;
    return gradHy;
}

double computeAnalyticalGradGx(const double x, const double y, const elementStruct &element)
{
    const double xi = x, yi = y;
    const double xj = element.x1, yj = element.y1;
    const double xjp1 = element.x2, yjp1 = element.y2;

    const double lj = element.length;
    const double lj_sq = lj * lj; // Utilisation du membre length au carré

    // Dénominateur atan (le même que dans ta formule de G)
    const double den_atan = xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1);

    double atan_sum = 0;
    if(std::abs(den_atan) > 1e-15) {
        atan_sum = std::atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1)) / den_atan) +
                   std::atan(((xi - xjp1)*(xjp1 - xj) + (yi - yjp1)*(yjp1 - yj)) / den_atan);
    }

    const double rj_sq = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);
    const double rjp1_sq = (xi-xjp1)*(xi-xjp1) + (yi-yjp1)*(yi-yjp1);

    const double log_diff = std::log(rj_sq) - std::log(rjp1_sq);

    double gradGx = (lj * (-2.0 * (yj - yjp1) * atan_sum - (xj - xjp1) * log_diff)) / (4.0 * M_PI * lj_sq);

    return gradGx;
}

double computeAnalyticalGradGy(const double x, const double y, const elementStruct &element)
{
    const double xi = x, yi = y;
    const double xj = element.x1, yj = element.y1;
    const double xjp1 = element.x2, yjp1 = element.y2;
    const double lj = element.length;
    const double lj_sq = lj * lj;

    const double den_atan = xjp1*yi + xi*yj - xjp1*yj - xi*yjp1 + xj*(-yi + yjp1);

    double atan_sum = 0;
    if(std::abs(den_atan) > 1e-15) {
        atan_sum = std::atan(((xi - xj)*(xj - xjp1) + (yi - yj)*(yj - yjp1)) / den_atan) +
                   std::atan(((xi - xjp1)*(xjp1 - xj) + (yi - yjp1)*(yjp1 - yj)) / den_atan);
    }

    const double rj_sq = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);
    const double rjp1_sq = (xi-xjp1)*(xi-xjp1) + (yi-yjp1)*(yi-yjp1);

    const double log_diff = std::log(rj_sq) - std::log(rjp1_sq);

    double gradGy = (lj * (2.0 * (xj - xjp1) * atan_sum - (yj - yjp1) * log_diff)) / (4.0 * M_PI * lj_sq);

    return gradGy;
}

// Computes the 'g' element numerically. "x" and "y" are the coordinates at which 'g' should be computed. "element" is the
// element 'j' while will helps to compute 'g'. "integrationType" is the type of numerical integration used (in case of numerical
// computation).
// --- G NUMÉRIQUE ---
double computeNumericalG(const double x, const double y,
                         const elementStruct &element,
                         const std::string integrationType,
                         const std::vector<double> &localCoord,
                         const std::vector<double> &weights)
{
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        double xi_param = localCoord[3*k];
        double gX = element.midX + element.deltaX * xi_param * 0.5;
        double gY = element.midY + element.deltaY * xi_param * 0.5;

        double r_sq = (x - gX)*(x - gX) + (y - gY)*(y - gY);
        if(r_sq > 1e-20) gSum += std::log(r_sq) * weights[k];
    }
    // log(r) = 0.5 * log(r_sq), d'où le (gSum * 0.5)
    return (gSum * 0.5) * element.length * (0.25 * M_1_PI);
}

// --- GRAD HX NUMÉRIQUE ---
double computeNumericalGradHx(const double x, const double y,
                              const elementStruct &element,
                              const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights)
{
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        double xi_param = localCoord[3*k];
        double hX = element.midX + element.deltaX * xi_param * 0.5;
        double hY = element.midY + element.deltaY * xi_param * 0.5;

        double dx = x - hX;
        double dy = y - hY;
        double r_sq = dx*dx + dy*dy;
        double r_4 = r_sq * r_sq;

        if(r_sq > 1e-20) {
            double term1 = -(0.5 * r_sq - dx*dx) / r_4;
            double term2 = (dx * dy) / r_4;
            gSum += (term1 * element.nx + term2 * element.ny) * weights[k];
        }
    }
    double gradHx = gSum * element.length * (0.5 * M_1_PI);
    return element.normalOutward ? gradHx : -gradHx;
}

// --- GRAD HY NUMÉRIQUE ---
double computeNumericalGradHy(const double x, const double y,
                              const elementStruct &element,
                              const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights)
{
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        double xi_param = localCoord[3*k];
        double hX = element.midX + element.deltaX * xi_param * 0.5;
        double hY = element.midY + element.deltaY * xi_param * 0.5;

        double dx = x - hX;
        double dy = y - hY;
        double r_sq = dx*dx + dy*dy;
        double r_4 = r_sq * r_sq;

        if(r_sq > 1e-20) {
            double term1 = -(0.5 * r_sq - dy*dy) / r_4;
            double term2 = (dx * dy) / r_4;
            gSum += (term1 * element.ny + term2 * element.nx) * weights[k];
        }
    }
    double gradHy = gSum * element.length * (0.5 * M_1_PI);
    return element.normalOutward ? gradHy : -gradHy;
}

// --- GRAD GX NUMÉRIQUE ---
double computeNumericalGradGx(const double x, const double y,
                              const elementStruct &element,
                              const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights)
{
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        double xi_param = localCoord[3*k];
        double gX = element.midX + element.deltaX * xi_param * 0.5;
        double gY = element.midY + element.deltaY * xi_param * 0.5;

        double dx = x - gX;
        double r_sq = dx*dx + (y - gY)*(y - gY);

        if(r_sq > 1e-20) gSum += (dx / r_sq) * weights[k];
    }
    return gSum * element.length * (0.25 * M_1_PI);
}

// --- GRAD GY NUMÉRIQUE ---
double computeNumericalGradGy(const double x, const double y,
                              const elementStruct &element,
                              const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights)
{
    double gSum = 0;
    for(size_t k = 0; k < weights.size(); k++)
    {
        double xi_param = localCoord[3*k];
        double gX = element.midX + element.deltaX * xi_param * 0.5;
        double gY = element.midY + element.deltaY * xi_param * 0.5;

        double dy = y - gY;
        double r_sq = (x - gX)*(x - gX) + dy*dy;

        if(r_sq > 1e-20) gSum += (dy / r_sq) * weights[k];
    }
    return gSum * element.length * (0.25 * M_1_PI);
}
