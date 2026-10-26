#include "functionsBEM.hpp"
#include <cmath>
#include <cassert>

// a*ln(rsq) avec limite 0 quand le point est sur un nœud (a ~ r -> 0)
inline double xlogx(const double a, const double rsq) { return (rsq > 0.0) ? a*std::log(rsq) : 0.0; }

// Signe d'orientation commun à toute la famille H (formules établies pour n = (ty,-tx)/L)
inline double orient(const elementStruct &e) { return e.normalOutward ? 1.0 : -1.0; }

// ------------------------------ famille G (indépendante de l'orientation) ------------------------------
double computeAnalyticalG(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    const double a1 =   g.tx*g.d2x + g.ty*g.d2y;     // L² - tx*u - ty*v
    const double a2 = -(g.tx*g.d1x + g.ty*g.d1y);    // tx*u + ty*v
    return (-2.0*g.L2 + 2.0*g.c*g.theta + xlogx(a1, g.r2sq) + xlogx(a2, g.r1sq))
           / (4.0*M_PI*g.L);
}

double computeAnalyticalGradGx(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    assert(g.r1sq > 0.0 && g.r2sq > 0.0);            // log-divergent aux nœuds
    const double th  = g.onElement ? 0.0 : g.theta;  // valeur principale
    const double ell = std::log(g.r1sq / g.r2sq);
    return (g.tx*ell - 2.0*g.ty*th) / (4.0*M_PI*g.L);
}

double computeAnalyticalGradGy(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    assert(g.r1sq > 0.0 && g.r2sq > 0.0);
    const double th  = g.onElement ? 0.0 : g.theta;
    const double ell = std::log(g.r1sq / g.r2sq);
    return (g.ty*ell + 2.0*g.tx*th) / (4.0*M_PI*g.L);
}

// ------------------------------ famille H (impaire : signe s) ------------------------------
double computeAnalyticalH(const double x, const double y, const elementStruct &e, const bool ccw)
{
    const SegGeom g = makeSegGeom(x, y, e);
    const double s  = ccw ? 1.0 : -1.0;
    const double th = g.onElement ? M_PI : g.theta;
    return s * th * (0.5*M_1_PI);
}

double computeAnalyticalGradHx(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    assert(g.r1sq > 0.0 && g.r2sq > 0.0);
    return orient(e) * (g.d2y/g.r2sq - g.d1y/g.r1sq) * (0.5*M_1_PI);
}

double computeAnalyticalGradHy(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    assert(g.r1sq > 0.0 && g.r2sq > 0.0);
    return orient(e) * (g.d1x/g.r1sq - g.d2x/g.r2sq) * (0.5*M_1_PI);
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

SecondKernels computeNumericalSecondKernels(
    const double x,
    const double y,
    const elementStruct &element,
    const std::string integrationType,
    const std::vector<double> &localCoord,
    const std::vector<double> &weights)
{
    (void)integrationType;

    assert(localCoord.size() >= 3 * weights.size());

    double GxxSum = 0.0;
    double GxySum = 0.0;
    double HxxSum = 0.0;
    double HxySum = 0.0;

    for(std::size_t k = 0; k < weights.size(); ++k)
    {
        const double xi = localCoord[3 * k];

        // Point d'intégration sur le segment
        const double X = element.midX + element.deltaX * xi * 0.5;
        const double Y = element.midY + element.deltaY * xi * 0.5;

        const double dx = x - X;
        const double dy = y - Y;
        const double r2 = dx * dx + dy * dy;

        // Singularité : cette quadrature ne convient pas
        // lorsque le point d'observation appartient au segment.
        if(r2 <= 1e-20)
            continue;

        const double invR2 = 1.0 / r2;
        const double invR4 = invR2 * invR2;
        const double invR6 = invR4 * invR2;

        const double dx2 = dx * dx;
        const double dy2 = dy * dy;

        // Famille G
        GxxSum += (dy2 - dx2) * invR4 * weights[k];
        GxySum -= 2.0 * dx * dy * invR4 * weights[k];

        // Famille H
        const double Ax = dx * (-dx2 + 3.0 * dy2) * invR6;
        const double Ay = dy * (-3.0 * dx2 + dy2) * invR6;

        const double Bx = Ay;
        const double By = -Ax;

        HxxSum += (element.nx * Ax + element.ny * Bx) * weights[k];
        HxySum += (element.nx * Ay + element.ny * By) * weights[k];
    }

    const double scaleG = element.length * (0.25 * M_1_PI);
    const double scaleH = element.length * (0.5 * M_1_PI);

    const double orientation = element.normalOutward ? 1.0 : -1.0;

    SecondKernels result;

    result.Gxx = scaleG * GxxSum;
    result.Gxy = scaleG * GxySum;
    result.Gyy = -result.Gxx;

    result.Hxx = orientation * scaleH * HxxSum;
    result.Hxy = orientation * scaleH * HxySum;
    result.Hyy = -result.Hxx;

    return result;
}


Kernels computeAnalyticalKernels(
    const double x,
    const double y,
    const elementStruct &e,
    const bool ccw)
{
    const SegGeom g = makeSegGeom(x, y, e);

    const double s  = ccw ? 1.0 : -1.0;
    const double th = g.onElement ? M_PI : g.theta;

    const double a1 = g.tx * g.d2x + g.ty * g.d2y;
    const double a2 = -(g.tx * g.d1x + g.ty * g.d1y);

    Kernels k;

    k.G = (-2.0 * g.L2
           + 2.0 * g.c * g.theta
           + xlogx(a1, g.r2sq)
           + xlogx(a2, g.r1sq))
          / (4.0 * M_PI * g.L);

    k.H = s * th * (0.5 * M_1_PI);

    return k;
}


// Noyaux d'ordre 1 ET 2 en une seule passe : une seule géométrie (makeSegGeom)
// par couple (point, élément). Mêmes formules que les fonctions séparées.
GradKernels computeAnalyticalGradKernels(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    assert(g.r1sq > 0.0 && g.r2sq > 0.0);

    const double s     = orient(e);
    const double c2pi  = 0.5 * M_1_PI;
    const double th    = g.onElement ? 0.0 : g.theta;   // valeur principale
    const double ell   = std::log(g.r1sq / g.r2sq);
    const double invG  = 1.0 / (g.r1sq * g.r2sq);
    const double i1    = 1.0 / (g.r1sq * g.r1sq), i2 = 1.0 / (g.r2sq * g.r2sq);

    GradKernels k;

    // ordre 1
    k.Gx = (g.tx*ell - 2.0*g.ty*th) / (4.0*M_PI*g.L);
    k.Gy = (g.ty*ell + 2.0*g.tx*th) / (4.0*M_PI*g.L);
    k.Hx = s * (g.d2y/g.r2sq - g.d1y/g.r1sq) * c2pi;
    k.Hy = s * (g.d1x/g.r1sq - g.d2x/g.r2sq) * c2pi;

    // ordre 2
    k.sk.Gxx =  g.L * (g.d1y*g.d2y - g.d1x*g.d2x) * invG * c2pi;
    k.sk.Gxy = -g.L * (g.d1x*g.d2y + g.d1y*g.d2x) * invG * c2pi;
    k.sk.Gyy = -k.sk.Gxx;
    k.sk.Hxx = s * (2.0*g.d2x*g.d2y*i2 - 2.0*g.d1x*g.d1y*i1) * c2pi;
    k.sk.Hxy = s * ((g.d1x*g.d1x - g.d1y*g.d1y)*i1 - (g.d2x*g.d2x - g.d2y*g.d2y)*i2) * c2pi;
    k.sk.Hyy = -k.sk.Hxx;

    return k;
}


// Noyaux d'ordre 1 seulement (E nodal) : une seule géométrie par couple (point, élément).
FirstKernels computeAnalyticalFirstKernels(const double x, const double y, const elementStruct &e)
{
    const SegGeom g = makeSegGeom(x, y, e);
    assert(g.r1sq > 0.0 && g.r2sq > 0.0);

    const double s    = orient(e);
    const double c2pi = 0.5 * M_1_PI;
    const double th   = g.onElement ? 0.0 : g.theta;   // valeur principale
    const double ell  = std::log(g.r1sq / g.r2sq);

    FirstKernels k;
    k.Gx = (g.tx*ell - 2.0*g.ty*th) / (4.0*M_PI*g.L);
    k.Gy = (g.ty*ell + 2.0*g.tx*th) / (4.0*M_PI*g.L);
    k.Hx = s * (g.d2y/g.r2sq - g.d1y/g.r1sq) * c2pi;
    k.Hy = s * (g.d1x/g.r1sq - g.d2x/g.r2sq) * c2pi;
    return k;
}
