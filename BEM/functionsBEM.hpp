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

// ============================================================
//  ELEMENT STRUCT
// ============================================================
struct elementStruct
{
    std::size_t tag;
    double      bcValue[2];     // [phi, phi_n]
    std::size_t nodeTags[2];
    std::size_t midNodeTag;
    double      length;
    double      x1, y1, x2, y2;
    double      midX, midY;
    double      deltaX, deltaY;
    double      nx, ny;
    bool normalOutward;

    // 'E' = external boundary, 'A' = interface primary side, 'B' = interface secondary side
    char        interfaceSide  = 'E';
    std::string interfaceName  = "";
    int         interfaceIndex = -1;
    int         globalIndex    = -1;
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

// ============================================================
//  SOLVER ENTRY POINTS
// ============================================================

void solverBEM(std::map<int, double> &electrostaticPressure,
               int &nbViews,
               std::map<int, std::pair<double,double>> &nodalDisplacements,
               const bool postProcessing,
               const int iteration,
               const bool untangleMesh);

void singleDomainSolverBEM(std::map<int, double> &electrostaticPressure,
                           const int nbViews,
                           std::map<int, std::pair<double,double>> &nodalDisplacements,
                           const bool postProcessing,
                           const int iteration,
                           const std::string domainName,
                           const int phiViewTag,
                           const int electricFieldViewTag,
                           const bool computePhi,
                           const bool computeElectricField);

void multiDomainSolverBEM(std::map<int, double> &electrostaticPressure,
                          int &nbViews,
                          std::map<int, std::pair<double,double>> &nodalDisplacements,
                          const bool postProcessing,
                          const int iteration,
                          std::vector<domainData> &domains,
                          std::vector<interfaceStruct> &interfaces,
                          const int phiViewTag,
                          const int electricFieldViewTag,
                          const bool computePhi,
                          const bool computeElectricField);

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

// ============================================================
//  SYSTEM ASSEMBLY
// ============================================================

void fillSystem(Eigen::MatrixXd &A,
                Eigen::MatrixXd &B,
                Eigen::MatrixXd &c,
                const std::vector<elementStruct> &elementVector,
                const std::size_t nbDirichletElements,
                const std::string integrationType,
                const bool solutionType);

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
//  ELECTROSTATIC PRESSURE
// ============================================================

void computeElectrostaticPressure(std::map<int, double> &electrostaticPressure,
                                  std::vector<elementStruct> &elementVector,
                                  const Eigen::MatrixXd &x,
                                  const std::size_t nbDirichletElements,
                                  const std::string domainName);

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
                        std::map<std::size_t, std::size_t> &domainNodeIndices,
                        const std::vector<elementStruct> &elementVector,
                        const std::string integrationType,
                        std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
                        const bool &solutionType);

void elementNodeData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                     const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                     const std::vector<std::vector<int>> &domainElementTypes,
                     std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                     const std::vector<elementStruct> &elementVector,
                     std::map<int, double> &internalNodalPhi,
                     const int domainPhysicalGroupTag,
                     std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements);

void phiInitialization(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                       std::map<std::size_t, std::vector<std::size_t>> &position,
                       const std::vector<std::vector<int>> &domainElementTypes,
                       const std::vector<std::vector<std::vector<size_t>>> &domainElementTags);

void boundaryElementTags(std::vector<std::size_t> &domainBoundaryElementTags,
                         const std::vector<std::vector<int>> &domainElementTypes,
                         const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                         const std::vector<elementStruct> &elementVector,
                         const int domainPhysicalGroupTag);

void elementNodalPhiAssociation(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                                const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                                std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                                std::map<int, double> &internalNodalPhi);

void computeElementData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                        std::vector<std::vector<std::vector<std::vector<double>>>> &electricField,
                        const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                        const std::vector<std::vector<int>> &domainElementTypes,
                        const std::vector<int> &domainEntityTags,
                        const std::vector<elementStruct> &elementVector,
                        const bool solutionType,
                        const std::string integrationType,
                        const bool computePhi,
                        const bool computeElectricField);

void displayData(const int viewTag,
                 const std::string modelName,
                 const std::string fileName,
                 const std::string dataType,
                 const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                 const std::vector<std::vector<std::vector<std::vector<double>>>> &data,
                 const int viewIndex);

// ============================================================
//  UTILITY
// ============================================================

void displayTime(const double start,
                 const double end,
                 std::string functionName,
                 const bool showTime);

void convergenceFun(const std::vector<int> &domainEntityTags,
                    const std::vector<std::vector<int>> &domainElementTypes,
                    const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> &elementalPhi,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> &elementalElectricField);

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

void boundaryElementTags(std::vector<size_t> &domainBoundaryElTags,
                         const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                         const std::vector<std::vector<int>> &domainElementTypes,
                         const std::map<size_t, std::vector<size_t>> &domainElementNodeTags,
                         const std::map<size_t, std::vector<size_t>> &bieNodeMap);

struct Force2D {
    double fx = 0.0;
    double fy = 0.0;
};
Force2D computeDEPForceMST(const std::vector<elementStruct> &elementVector,
                           const std::map<int, double> &internalNodalPhi, // Utilise ta map résultante
                           const std::string &interfaceName,
                           double epsilon_oil);
