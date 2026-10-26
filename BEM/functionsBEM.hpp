#pragma once

#include <gmsh.h>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <omp.h>
#include <algorithm>
#include <cstdlib>
#include <array>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include "parameters.hpp"

// ============================================================
//  ELEMENT STRUCT
// ============================================================

struct elementStruct
{
    std::size_t tag = 0;
    double bcValue[2] = {HUGE_VAL, HUGE_VAL}; // [phi, phi_n]
    std::size_t nodeTags[2] = {0, 0};
    std::size_t midNodeTag = 0;

    double length = 0.0;
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    double midX = 0.0, midY = 0.0;
    double deltaX = 0.0, deltaY = 0.0;
    double nx = 0.0, ny = 0.0;
    bool normalOutward = true;

    // Domaine auquel appartient cette représentation BIE
    int domainPhysicalGroupTag = -1;

    // E = external, A = interface primary side,
    // B = interface secondary side, H = hole
    char interfaceSide = 'E';
    std::string interfaceName = "";
    int interfaceIndex = -1;
    int globalIndex = -1;

    double En = HUGE_VAL;
    double Et = HUGE_VAL;
    double Ex = HUGE_VAL;
    double Ey = HUGE_VAL;
};

// ============================================================
//  INTERFACE STRUCT
// ============================================================
struct interfaceStruct
{
    std::string name;
    std::string domainA;
    std::string domainB;
    double      kRatio = 1.0;           // epsilon_A / epsilon_B

    std::vector<std::size_t> elemIndicesA;
    std::vector<std::size_t> elemIndicesB;

    std::size_t nbElements() const { return elemIndicesA.size(); }
};

// ============================================================
//  PER-DOMAIN DATA
// ============================================================
struct domainData
{
    std::string name;
    int    physicalGroupDim = -1;
    int    physicalGroupTag = -1;
    double epsilon          = 8.8541878128e-12;

    std::vector<elementStruct> elementVector;

    std::size_t nbDirichletElements = 0;
    std::size_t nbExternalElements  = 0;
    std::size_t rowOffset           = 0;

    std::size_t nbRows() const { return elementVector.size(); }
};

struct DomainMesh {
    std::vector<std::size_t> nodeTags;
    std::map<std::size_t, std::size_t> nodeIndices;
    std::vector<double> nodeCoords;
    std::vector<int> entityTags;
    std::vector<std::vector<int>> elemTypes;
    std::vector<std::vector<std::vector<std::size_t>>> elemTags;
    std::map<std::size_t, std::vector<std::size_t>> elemNodeTags;
    std::vector<std::vector<std::vector<double>>> barycenters; // [entité][type] -> xyz à plat
    std::vector<std::vector<int>> nbNodes; // [entité][type] -> nb de nœuds par élément
};

DomainMesh loadDomainMesh(const domainData &dom);

constexpr double SEG_TOL = 1e-9;   // relatif à L²

struct SegGeom {
    double tx, ty, L, L2;          // t = P_{j+1} - P_j
    double d1x, d1y, d2x, d2y;     // P_j - P_i , P_{j+1} - P_i   (d1 = (-u,-v))
    double r1sq, r2sq;
    double c, dot;                 // c = d1 x d2 = tx*v - ty*u ; dot = d1.d2
    double theta;                  // angle signé d1 -> d2, dans (-pi, pi]
    bool   onElement;              // point à l'intérieur du segment
};

inline SegGeom makeSegGeom(const double x, const double y, const elementStruct &e)
{
    SegGeom g;
    g.tx = e.x2 - e.x1;  g.ty = e.y2 - e.y1;
    g.L2 = g.tx*g.tx + g.ty*g.ty;
    g.L  = std::sqrt(g.L2);
    g.d1x = e.x1 - x;  g.d1y = e.y1 - y;
    g.d2x = e.x2 - x;  g.d2y = e.y2 - y;
    g.r1sq = g.d1x*g.d1x + g.d1y*g.d1y;
    g.r2sq = g.d2x*g.d2x + g.d2y*g.d2y;
    g.c   = g.d1x*g.d2y - g.d1y*g.d2x;
    g.dot = g.d1x*g.d2x + g.d1y*g.d2y;
    g.theta = std::atan2(g.c, g.dot);
    g.onElement = (std::abs(g.c) <= SEG_TOL * g.L2) && (g.dot < 0.0);
    return g;
}

// ---- Noyaux d'ordre 2 (un seul makeSegGeom par paire point/élément) ----
struct SecondKernels { double Gxx, Gxy, Gyy, Hxx, Hxy, Hyy; };

SecondKernels computeNumericalSecondKernels(const double x,
                                            const double y,
                                            const elementStruct &element,
                                            const std::string integrationType,
                                            const std::vector<double> &localCoord,
                                            const std::vector<double> &weights);

struct Kernels { double G; double H; };
Kernels computeAnalyticalKernels(const double x,
                                 const double y,
                                 const elementStruct &e,
                                 const bool ccw);

// Noyaux d'ordre 1 et 2 calculés ensemble (une seule géométrie par couple point/élément)
struct GradKernels { double Hx, Hy, Gx, Gy; SecondKernels sk; };
GradKernels computeAnalyticalGradKernels(const double x, const double y,
                                         const elementStruct &element);

// Noyaux d'ordre 1 seulement (champ E nodal)
struct FirstKernels { double Hx, Hy, Gx, Gy; };
FirstKernels computeAnalyticalFirstKernels(const double x, const double y,
                                           const elementStruct &element);

struct Quadrature {
    std::string type;
    std::vector<double> localCoord, weights;   // localCoord : 3 valeurs (u,v,w) par point
};

inline Quadrature loadQuadrature(const std::string &type)
{
    Quadrature q;
    q.type = type;
    gmsh::model::mesh::getIntegrationPoints(1, type, q.localCoord, q.weights);
    return q;
}

// ---- Vues supplémentaires : tenseur de Maxwell + grad|E|² ----
using NodalScalarData = std::vector<std::vector<std::vector<std::vector<double>>>>;
using NodalVectorData = std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>;

struct extraViewsStruct
{
    // Champs à publier : recopiés depuis BEMParams::fields (parameters.hpp)
    std::set<Field> fields;

    bool wants(const Field f) const { return fields.count(f) > 0; }

    // Tags Gmsh
    int phiViewTag = -1;
    int electricFieldViewTag = -1;
    int electricFieldNodeViewTag = -1;

    int txxViewTag = -1;
    int tyyViewTag = -1;
    int txyViewTag = -1;
    int tzzViewTag = -1;
    int gradE2ViewTag = -1;
    int stressJumpViewTag = -1;

    void createViews()
    {
        // Dépendance : le tenseur de Maxwell se calcule à partir du champ E nodal
        if(wants(Field::Maxwell))
            fields.insert(Field::ENodal);

        // Réinitialisation des tags
        phiViewTag = -1;
        electricFieldViewTag = -1;
        electricFieldNodeViewTag = -1;
        txxViewTag = -1;
        tyyViewTag = -1;
        txyViewTag = -1;
        tzzViewTag = -1;
        gradE2ViewTag = -1;
        stressJumpViewTag = -1;

        if(wants(Field::Phi))
            phiViewTag = gmsh::view::add("Phi");

        if(wants(Field::ENodal))
            electricFieldNodeViewTag =
                gmsh::view::add("Electric Field (Nodal)");

        if(wants(Field::EBary))
            electricFieldViewTag =
                gmsh::view::add("Electric Field (Barycenters)");

        if(wants(Field::Maxwell))
        {
            txxViewTag = gmsh::view::add("Maxwell Txx [Pa]");
            tyyViewTag = gmsh::view::add("Maxwell Tyy [Pa]");
            txyViewTag = gmsh::view::add("Maxwell Txy [Pa]");
            tzzViewTag = gmsh::view::add("Maxwell Tzz [Pa]");
        }

        if(wants(Field::GradE2))
            gradE2ViewTag = gmsh::view::add("grad |E|^2 [V^2/m^3]");

        if(wants(Field::StressJump))
            stressJumpViewTag =
                gmsh::view::add("Maxwell Normal Stress Jump [Pa]");

        std::cout << "[Gmsh views] Phi=" << phiViewTag
                  << ", E=" << electricFieldViewTag
                  << ", E nodal=" << electricFieldNodeViewTag
                  << ", Txx=" << txxViewTag
                  << ", gradE2=" << gradE2ViewTag
                  << ", stress jump=" << stressJumpViewTag
                  << '\n';
    }
};

void computeMaxwellStressNodal(NodalScalarData &Txx, NodalScalarData &Tyy, NodalScalarData &Txy, NodalScalarData &Tzz, const NodalVectorData &nodalE, const double epsilon);

void computeElementGradE2(std::vector<std::vector<std::vector<std::vector<double>>>> &gradE2,
                          const DomainMesh &mesh,
                          const std::vector<elementStruct> &elementVector,
                          const Quadrature &quad,
                          const bool &solutionType,
                          const double prefactor = 1.0);

// À appeler HORS de tout bloc "omp critical(gmsh_api)" : il pose ses propres sections critiques.
void computeAndDisplayExtraViews(const extraViewsStruct &extra,
                                 const std::string &modelName,
                                 const DomainMesh &mesh,
                                 const NodalVectorData &nodalE,
                                 const std::vector<elementStruct> &elementVector,
                                 const Quadrature &quad,
                                 const double epsilon,
                                 const bool &solutionType);

void postProcessDomain(const domainData &dom,
                       const DomainMesh &mesh,
                       const Quadrature &quad,
                       const extraViewsStruct &extra,
                       const std::string &modelName,
                       const bool solutionType);

// ============================================================
//  SOLVER ENTRY POINTS
// ============================================================

void solverBEM(std::map<int, double> &electrostaticPressure,
               std::map<int, std::pair<double,double>> &nodalDisplacements,
               const bool postProcessing,
               const int iteration,
               const bool untangleMesh,
               const BEMParams &params = BEMParams());

void multiDomainSolverBEM(std::map<int,double> &electrostaticPressure,
                          std::map<int,std::pair<double,double>> &nodalDisplacements,
                          const bool postProcessing,
                          const int iteration,
                          std::vector<domainData> &domains,
                          std::vector<interfaceStruct> &interfaces,
                          const extraViewsStruct &extra,
                          const BEMParams &params);

// ============================================================
//  INTERFACE DETECTION
// ============================================================

void detectInterfaces(std::vector<interfaceStruct> &interfaces,
                      const std::vector<std::string> &domainNames);

// ============================================================
//  BOUNDARY ELEMENT CONSTRUCTION
// ============================================================

void getBC(std::map<std::string, std::pair<bool, double>> &boundaryConditions,
           const std::string domainName);

double getEpsilon(const std::string domainName);

void setViewDefaults(const int tag);

std::size_t fillElementVector(std::vector<elementStruct> &elementVector,
                              int &domainPhysicalGroupDim,
                              int &domainPhysicalGroupTag,
                              std::map<int, std::pair<double,double>> &nodalDisplacements,
                              const std::string domainName);

void addInterfaceElements(domainData &domainA,
                          domainData &domainB,
                          interfaceStruct &iface,
                          std::map<int, std::pair<double,double>> &nodalDisplacements);

void computeGeometricParam(elementStruct &element,
                           std::map<std::size_t, std::size_t> &nodeIndices,
                           const std::vector<double> &nodeCoords,
                           const std::pair<double,double> &nodalDisp1,
                           const std::pair<double,double> &nodalDisp2);

void assembleGlobalSystem(Eigen::MatrixXd &A_global,
                          Eigen::MatrixXd &b_global,
                          const std::vector<domainData> &domains,
                          const std::vector<interfaceStruct> &interfaces,
                          const std::string integrationType,
                          const bool solutionType,
                          std::size_t &ifPotentialOffset,
                          std::size_t &ifFluxOffset);

void distributeGlobalSolution(std::vector<domainData> &domains,
                              const std::vector<interfaceStruct> &interfaces,
                              const Eigen::MatrixXd &X_global,
                              const std::size_t ifPotentialOffset,
                              const std::size_t ifFluxOffset);

// ============================================================
//  POST-PROCESSING
// ============================================================

void getModel(std::vector<std::size_t> &domainNodeTags,
              std::map<std::size_t, std::size_t> &domainNodeIndices,
              std::vector<double> &domainNodeCoords,
              std::vector<int> &domainEntityTags,
              std::vector<std::vector<int>> &domainElementTypes,
              std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
              std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
              const int domainPhysicalGroupDim,
              const int domainPhysicalGroupTag);

void computeNodeToElements(std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
                           const std::vector<elementStruct> &elementVector);

// NOTE: all post-processing functions take the FULL elementVector.
// The BIE filtering ('E'+'A' only, no 'B') is handled internally.
void computeInternalPhi(std::map<int, double> &internalNodalPhi,
                        const std::vector<double> &domainNodeCoords,
                        const std::vector<std::size_t> &domainNodeTags,
                        const std::map<std::size_t, std::size_t> &domainNodeIndices,
                        const std::vector<elementStruct> &elementVector,
                        const Quadrature &quad,
                        std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
                        const bool &solutionType);

void elementNodeData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                     const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                     const std::vector<std::vector<int>> &domainElementTypes,
                     const std::vector<std::vector<int>> &nbNodesPerType,
                     const std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                     const std::vector<elementStruct> &elementVector,
                     std::map<int, double> &internalNodalPhi,
                     const int domainPhysicalGroupTag,
                     std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements);

void phiInitialization(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                       std::map<std::size_t, std::vector<std::size_t>> &position,
                       const std::vector<std::vector<int>> &domainElementTypes,
                       const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                       const std::vector<std::vector<int>> &nbNodesPerType);

void elementNodalPhiAssociation(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                                const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                                const std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                                std::map<int, double> &internalNodalPhi);

void computeElementData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                        std::vector<std::vector<std::vector<std::vector<double>>>> &electricField,
                        const DomainMesh &mesh,
                        const std::vector<elementStruct> &elementVector,
                        const bool solutionType,
                        const Quadrature &quad,
                        const bool computePhi,
                        const bool computeElectricField);

void displayData(const int viewTag,
                 const std::string modelName,
                 const std::string fileName,
                 const std::string dataType,
                 const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                 const std::vector<std::vector<std::vector<std::vector<double>>>> &data);

// ============================================================
//  UTILITY
// ============================================================

void displayTime(const double start,
                 const double end,
                 std::string functionName,
                 const bool showTime);

// ============================================================
//  KERNEL COMPUTATIONS
// ============================================================

double computeAnalyticalH(const double x, const double y, const elementStruct &element, const bool ccw = true);
double computeAnalyticalG(const double x, const double y, const elementStruct &element);
double computeAnalyticalGradHx(const double x, const double y, const elementStruct &element);
double computeAnalyticalGradHy(const double x, const double y, const elementStruct &element);
double computeAnalyticalGradGx(const double x, const double y, const elementStruct &element);
double computeAnalyticalGradGy(const double x, const double y, const elementStruct &element);

double computeNumericalG(const double x, const double y, const elementStruct &element,
                         const std::string integrationType,
                         const std::vector<double> &localCoord,
                         const std::vector<double> &weights);
double computeNumericalGradHx(const double x, const double y, const elementStruct &element,
                               const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights);
double computeNumericalGradHy(const double x, const double y, const elementStruct &element,
                               const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights);
double computeNumericalGradGx(const double x, const double y, const elementStruct &element,
                               const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights);
double computeNumericalGradGy(const double x, const double y, const elementStruct &element,
                               const std::string integrationType,
                              const std::vector<double> &localCoord,
                              const std::vector<double> &weights);
void computeNormalOrientation(std::vector<elementStruct> &elementVector,
                               const int domainPhysicalGroupTag);

  void computeElectrostaticPressureMultiDomain(
          std::map<int, double>            &electrostaticPressure,
          const std::vector<domainData>    &domains);

  void computeElectrostaticPressureJump(
          std::map<int, double>         &electrostaticPressure,
          const std::vector<domainData> &domains,
          const interfaceStruct         &iface);

void computeInternalElectricField(
    std::map<std::size_t, std::vector<double>> &internalNodalE,
    const std::vector<double> &domainNodeCoords,
    const std::vector<std::size_t> &domainNodeTags,
    const std::map<std::size_t, std::size_t> &domainNodeIndices,
    const std::vector<elementStruct> &elementVector,
    const Quadrature &quad,
    std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
    const bool &solutionType
);

void computeBoundaryElectricField(
    std::vector<elementStruct> &elementVector
);

void computeCartesianElectricField(
    std::vector<elementStruct> &elementVector
);

void elementNodeElectricField(
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
        &electricFieldNode,
    const std::vector<std::vector<std::vector<std::size_t>>>
        &domainElementTags,
    const std::map<std::size_t, std::vector<std::size_t>>
        &domainElementNodeTags,
    const std::vector<elementStruct> &elementVector,
    const std::map<std::size_t, std::vector<double>> &internalNodalE
);

void displayNodeElectricField(
    const int viewTag,
    const std::string modelName,
    const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
    const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &data);

void printInterfaceFieldDiagnostics(
    const interfaceStruct& iface,
    const domainData& domainA,
    const domainData& domainB,
    const double epsA,
    const double epsB);

struct maxwellLoadStruct
{
    double Fx = 0.0;  // N si les dimensions sont en m
    double Fy = 0.0;
};

maxwellLoadStruct computeInclusionMaxwellForce(
    const domainData& domainA,
    const domainData& domainB,
    const interfaceStruct& iface,
    const double thickness);

double maxwellNormalStress(
    const elementStruct& el,
    const double epsilon,
    const double nx,
    const double ny);

void displayInterfaceStressJump(
    const int viewTag,
    const std::string& modelName,
    const interfaceStruct& iface,
    const domainData& domainA,
    const std::map<int, double>& electrostaticPressure);


