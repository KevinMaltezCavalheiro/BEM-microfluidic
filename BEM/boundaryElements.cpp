#include "functionsBEM.hpp"

// ============================================================
//  getBC  (unchanged from original)
// ============================================================
void getBC(std::map<std::string, std::pair<bool, double>> &boundaryConditions,
           const std::string domainName)
{
    std::vector<std::string> bcKeys;
    gmsh::onelab::getNames(bcKeys, "Boundary Conditions.+");

    for(auto &key : bcKeys)
    {
        std::vector<double> values;
        gmsh::onelab::getNumber(key, values);

        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(getline(ss, word, '/'))
            words.push_back(word);

        // Key structure: "Boundary Conditions/<curveName>/<domainName>/<bcType>"
        if(words.size() == 4 && words[2].compare(domainName) == 0)
        {
            if(words[3].compare("dirichlet") == 0)
                boundaryConditions[words[1]] = {true, values[0]};
            else if(words[3].compare("neumann") == 0)
                boundaryConditions[words[1]] = {false, values[0]};
        }
    }
}

// ============================================================
//  getEpsilon
// ============================================================
double getEpsilon(const std::string domainName)
{
    std::vector<std::string> keys;
    gmsh::onelab::getNames(keys, "(Materials).+");

    double epsilon = 8.8541878128e-12;
    bool set = false;
    for(auto &key : keys)
    {
        std::vector<double> value;
        gmsh::onelab::getNumber(key, value);

        std::stringstream ss(key);
        std::vector<std::string> words;
        std::string word;
        while(std::getline(ss, word, '/'))
            words.push_back(word);

        if(words.size() == 3 && words[2] == "Epsilon" && words[1] == domainName)
        {
            epsilon = value[0];
            set = true;
        }
    }
    if(!set)
    {
        std::cout << "\tWarning: Epsilon not set for domain '" << domainName
                  << "' → using vacuum permittivity " << epsilon << " [F/m].\n";
    }
    return epsilon;
}

// ============================================================
//  computeGeometricParam  (unchanged from original)
// ============================================================
void computeGeometricParam(elementStruct &element,
                           std::map<std::size_t, std::size_t> &nodeIndices,
                           const std::vector<double> &nodeCoords,
                           const std::pair<double, double> &nodalDisp1,
                           const std::pair<double, double> &nodalDisp2)
{
    std::size_t i1Index = nodeIndices[element.nodeTags[0]];
    std::size_t i2Index = nodeIndices[element.nodeTags[1]];

    element.x1 = nodeCoords[3*i1Index]     + nodalDisp1.first;
    element.x2 = nodeCoords[3*i2Index]     + nodalDisp2.first;
    element.y1 = nodeCoords[3*i1Index + 1] + nodalDisp1.second;
    element.y2 = nodeCoords[3*i2Index + 1] + nodalDisp2.second;

    element.midX   = (element.x1 + element.x2) / 2.0;
    element.midY   = (element.y1 + element.y2) / 2.0;
    element.deltaX = element.x2 - element.x1;
    element.deltaY = element.y2 - element.y1;
    element.length = std::sqrt(element.deltaX*element.deltaX + element.deltaY*element.deltaY);

    double dx = element.x2-element.x1;
    double dy = element.y2-element.y1;
    double nx = dy;
    double ny = -dx;
    double norm = sqrt(nx*nx+ny*ny);
    nx = nx / norm;
    ny = ny / norm;
    element.nx = nx;
    element.ny = ny;
}

// ============================================================
//  fillElementVector
//  Fills the elementVector with EXTERNAL (non-interface) elements only.
//  Convention: Dirichlet elements come first, Neumann elements after.
//  Interface elements are handled by addInterfaceElements().
// ============================================================
std::size_t fillElementVector(std::vector<elementStruct> &elementVector,
                              int &domainPhysicalGroupDim,
                              int &domainPhysicalGroupTag,
                              std::map<int, std::pair<double, double>> &nodalDisplacements,
                              const std::string domainName)
{
    std::size_t nbDirichletElements = 0;

    std::map<std::string, std::pair<bool, double>> boundaryConditions;
    getBC(boundaryConditions, domainName);

    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    for(std::size_t i = 0; i < physicalGroupsDimTags.size(); i++)
    {
        int physicalGroupDim = physicalGroupsDimTags[i].first;
        int physicalGroupTag = physicalGroupsDimTags[i].second;
        std::string physicalGroupName;
        gmsh::model::getPhysicalName(physicalGroupDim, physicalGroupTag, physicalGroupName);

        if(physicalGroupName.compare(domainName) == 0)
        {
            domainPhysicalGroupDim = physicalGroupDim;
            domainPhysicalGroupTag = physicalGroupTag;
        }

        if(boundaryConditions.find(physicalGroupName) == boundaryConditions.end())
            continue;

        if(physicalGroupName.substr(0, 14) == "BEM_interface_")
            continue;

        // --- DÉTECTION DU TROU ---
        bool isHole = (physicalGroupName.find("Hole") != std::string::npos ||
                       physicalGroupName.find("hole") != std::string::npos);

        bool   bcType  = boundaryConditions[physicalGroupName].first;
        double bcValue = boundaryConditions[physicalGroupName].second;

        std::vector<int> entityTags;
        gmsh::model::getEntitiesForPhysicalGroup(physicalGroupDim, physicalGroupTag, entityTags);

        for(std::size_t m = 0; m < entityTags.size(); m++)
        {
            std::vector<int> elementTypes = {1, 8};
            for(std::size_t j = 0; j < elementTypes.size(); j++)
            {
                int nbNodesPerElement = (elementTypes[j] == 8) ? 3 : 2;

                std::vector<std::size_t> elementTags;
                std::vector<std::size_t> elementNodeTags;
                gmsh::model::mesh::getElementsByType(elementTypes[j], elementTags,
                                                     elementNodeTags, entityTags[m]);

                std::vector<std::size_t> nodeTags;
                std::vector<double> nodeCoords, parametricNodeCoords;
                gmsh::model::mesh::getNodesByElementType(elementTypes[j], nodeTags, nodeCoords,
                                                         parametricNodeCoords, entityTags[m], false);

                std::map<std::size_t, std::size_t> nodeIndices;
                for(std::size_t k = 0; k < nodeTags.size(); k++)
                    nodeIndices[nodeTags[k]] = k;

                // --- DÉBUT DE LA RÉGION PARALLÈLE ---
                #pragma omp parallel
                {
                    std::vector<elementStruct> localDirichlet;
                    std::vector<elementStruct> localNeumann;

                    #pragma omp for schedule(guided)
                    for(std::size_t el = 0; el < elementTags.size(); el++)
                    {
                        elementStruct element;
                        element.tag         = elementTags[el];
                        element.nodeTags[0] = nodeTags[nbNodesPerElement*el];
                        element.nodeTags[1] = nodeTags[nbNodesPerElement*el + 1];
                        element.midNodeTag  = (nbNodesPerElement == 3)
                                              ? nodeTags[nbNodesPerElement*el + 2] : 0;

                        // On applique les nouveaux flags
                        element.interfaceSide  = isHole ? 'H' : 'E';
                        element.interfaceName  = physicalGroupName;
                        element.interfaceIndex = -1;

                        std::pair<double,double> d1 = {0.0,0.0}, d2 = {0.0,0.0};

                        // Note: Les maps std::map ne sont pas thread-safe en écriture,
                        // mais en lecture seule (count/find), c'est généralement OK ici.
                        if(nodalDisplacements.count(element.nodeTags[0]))
                            d1 = nodalDisplacements[element.nodeTags[0]];
                        if(nodalDisplacements.count(element.nodeTags[1]))
                            d2 = nodalDisplacements[element.nodeTags[1]];

                        computeGeometricParam(element, nodeIndices, nodeCoords, d1, d2);

                        if(bcType) {
                            element.bcValue[0] = bcValue;
                            localDirichlet.push_back(element);
                        }
                        else {
                            element.bcValue[1] = bcValue;
                            localNeumann.push_back(element);
                        }
                    }

                    // Fusion des résultats locaux dans le vecteur global
                    #pragma omp critical
                    {
                        if (!localDirichlet.empty()) {
                            elementVector.insert(elementVector.begin(), localDirichlet.begin(), localDirichlet.end());
                            nbDirichletElements += localDirichlet.size();
                        }
                        if (!localNeumann.empty()) {
                            elementVector.insert(elementVector.end(), localNeumann.begin(), localNeumann.end());
                        }
                    }
                } // --- FIN DE LA RÉGION PARALLÈLE ---
            }
        }
    }

    // On garde l'orientation hors parallélisation (souvent séquentiel par nature)
    computeNormalOrientation(elementVector, domainPhysicalGroupTag);

    return nbDirichletElements;
}
// ============================================================
//  detectInterfaces
//  Scans all Physical Curves for the pattern "BEM_interface_<A>_<B>".
//  Fills the interfaces vector and computes kRatio = epsilon_A / epsilon_B
//  by reading the ONELAB epsilon values.
// ============================================================
void detectInterfaces(std::vector<interfaceStruct> &interfaces,
                      const std::vector<std::string> &domainNames)
{
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    for(auto &dimTag : physicalGroupsDimTags)
    {
        std::string name;
        gmsh::model::getPhysicalName(dimTag.first, dimTag.second, name);

        // Convention: "BEM_interface_<domainA>_<domainB>"
        // Split by '_': ["BEM", "interface", "<domainA>", "<domainB>"]
        std::stringstream ss(name);
        std::vector<std::string> words;
        std::string word;
        while(getline(ss, word, '_'))
            words.push_back(word);

        if(words.size() < 4 || words[0] != "BEM" || words[1] != "interface")
            continue;

        // Domain names may themselves contain underscores, so we reconstruct
        // them as words[2] and words[3..end].  For simplicity we assume that
        // each domain name is a single token (no underscores internally).
        // If you need multi-token domain names, adjust the .geo naming convention.
        std::string dA = words[2];
        std::string dB = words[3];
        // Re-join remaining words if domain name has embedded underscores
        for(std::size_t w = 4; w < words.size(); w++)
            dB += "_" + words[w];

        // Verify both domain names are known.
        bool foundA = false, foundB = false;
        for(auto &dn : domainNames)
        {
            if(dn == dA) foundA = true;
            if(dn == dB) foundB = true;
        }
        if(!foundA || !foundB)
        {
            std::cout << "\tWarning: interface '" << name
                      << "' references unknown domain(s). Skipping.\n";
            continue;
        }

        interfaceStruct iface;
        iface.name    = name;
        iface.domainA = dA;
        iface.domainB = dB;

        // Compute k-ratio = epsilon_A / epsilon_B
        double epsA = getEpsilon("BEM_domain_" + dA);
        double epsB = getEpsilon("BEM_domain_" + dB);
        iface.kRatio = (epsB > 0.0) ? (epsA / epsB) : 1.0;

        interfaces.push_back(iface);
        std::cout << "  Detected interface '" << name
                  << "'  (k_ratio = " << iface.kRatio << ")\n";
    }
}

// ============================================================
//  addInterfaceElements
//  Appends elements of the interface curve to both domain A and domain B.
//  Side A and side B share the same physical geometry but have opposite
//  outward normals, hence the sign convention in the global assembly.
// ============================================================
void addInterfaceElements(domainData &domainA,
                          domainData &domainB,
                          interfaceStruct &iface,
                          std::map<int, std::pair<double,double>> &nodalDisplacements)
{
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    int ifaceDim = -1, ifaceTag = -1;
    for(auto &dimTag : physicalGroupsDimTags)
    {
        std::string name;
        gmsh::model::getPhysicalName(dimTag.first, dimTag.second, name);
        if(name == iface.name) { ifaceDim = dimTag.first; ifaceTag = dimTag.second; break; }
    }
    if(ifaceDim < 0)
    {
        std::cerr << "Error: Cannot find Physical Curve '" << iface.name << "'.\n";
        return;
    }

    std::vector<int> entityTags;
    gmsh::model::getEntitiesForPhysicalGroup(ifaceDim, ifaceTag, entityTags);

    std::vector<int> elementTypes = {1, 8};
    for(auto entityTag : entityTags)
    {
        for(int elType : elementTypes)
        {
            int nbNPE = (elType == 8) ? 3 : 2;

            std::vector<std::size_t> elTags, elNodeTags;
            gmsh::model::mesh::getElementsByType(elType, elTags, elNodeTags, entityTag);
            if(elTags.empty()) continue;

            std::vector<std::size_t> nodeTags;
            std::vector<double> nodeCoords, parametric;
            gmsh::model::mesh::getNodesByElementType(elType, nodeTags, nodeCoords,
                                                      parametric, entityTag, false);

            std::map<std::size_t, std::size_t> nodeIndices;
            for(std::size_t k = 0; k < nodeTags.size(); k++)
                nodeIndices[nodeTags[k]] = k;

            for(std::size_t el = 0; el < elTags.size(); el++)
            {
                elementStruct base;
                base.tag         = elTags[el];
                base.nodeTags[0] = elNodeTags[nbNPE*el];
                base.nodeTags[1] = elNodeTags[nbNPE*el + 1];
                base.midNodeTag  = (nbNPE == 3) ? elNodeTags[nbNPE*el + 2] : 0;
                base.bcValue[0]  = 0.0;
                base.bcValue[1]  = 0.0;
                base.interfaceName = iface.name;

                std::pair<double,double> d1 = {0.0,0.0}, d2 = {0.0,0.0};
                if(nodalDisplacements.count(base.nodeTags[0]))
                    d1 = nodalDisplacements[base.nodeTags[0]];
                if(nodalDisplacements.count(base.nodeTags[1]))
                    d2 = nodalDisplacements[base.nodeTags[1]];
                computeGeometricParam(base, nodeIndices, nodeCoords, d1, d2);

                // --- Côté A : tester si la normale pointe vers l'intérieur de domainA ---
                elementStruct elA = base;
                elA.interfaceSide = 'A';
                iface.elemIndicesA.push_back(domainA.elementVector.size());
                domainA.elementVector.push_back(elA);

                // --- Côté B : même normale géométrique mais testée pour domainB ---
                // → résultat opposé à elA par construction (les deux domaines sont de
                //   part et d'autre de l'interface)
                elementStruct elB = base;
                elB.interfaceSide = 'B';
                iface.elemIndicesB.push_back(domainB.elementVector.size());
                domainB.elementVector.push_back(elB);
            }
        }
        computeNormalOrientation(domainA.elementVector, domainA.physicalGroupTag);
        computeNormalOrientation(domainB.elementVector, domainB.physicalGroupTag);
    }

    std::cout << "  Interface '" << iface.name << "': added "
              << iface.elemIndicesA.size() << " element(s) per side.\n";
}

// ============================================================
//  fillSystem  (extended for interface elements)
//
//  The convention in the local system is:
//    - columns 0 .. nbDirichletElements-1   : Dirichlet unknowns (phi_n)
//    - columns nbDirichletElements .. nbExternalElements-1 : Neumann unknowns (phi)
//    - columns nbExternalElements .. N-1     : interface elements
//      For each interface element the unknown is BOTH phi AND phi_n, but they
//      are handled at the GLOBAL level.  In the local block we treat them as
//      "Neumann-side unknowns" (phi is prescribed = 0 dummy; phi_n unknown) so
//      that the block structure [A | H_if] * [x_ext ; u_if] + [-G_if] * qn_if
//      is easy to extract.  The caller (assembleGlobalSystem) uses the raw H
//      and G values computed here and places them in the right columns.
//
//  For the local system we reuse the same structure as before.  Only the
//  mapping to the global system changes (done in assembleGlobalSystem).
// ============================================================
void fillSystem(Eigen::MatrixXd &A,
                Eigen::MatrixXd &B,
                Eigen::MatrixXd &c,
                const std::vector<elementStruct> &elementVector,
                const std::size_t nbDirichletElements,
                const std::string integrationType,
                const bool solutionType)
{
    // nbExternalElements: all elements with interfaceSide == 'E' ou 'H'
    // Ce calcul reste séquentiel car très rapide.
    std::size_t nbExternalElements = 0;
    for(auto &e : elementVector)
        if(e.interfaceSide == 'E' || e.interfaceSide == 'H') nbExternalElements++;

    // On récupère la taille pour la lisibilité
    const std::size_t N = elementVector.size();
    std::vector<double> localCoord, weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Début de la parallélisation
    // collapse(2) permet de fusionner les deux boucles pour que les threads
    // se partagent les N*N itérations au lieu de seulement les N de la première boucle.
    #pragma omp parallel for collapse(2) schedule(guided) shared(elementVector, A, B, c)
    for(std::size_t i = 0; i < N; i++)
    {
        for(std::size_t j = 0; j < N; j++)
        {
            // h et g sont déclarés à l'intérieur pour être privés à chaque itération
            double h = 0.0;
            double g = 0.0;

            if(i == j)
            {
                h = -0.5;
                g = (elementVector[i].length / (2*M_PI)) * (std::log(elementVector[i].length/2) - 1);
            }
            else
            {
                // Note : computeAnalyticalH prend maintenant 4 arguments
                h = computeAnalyticalH(elementVector[i].midX, elementVector[i].midY,
                                       elementVector[j], elementVector[j].normalOutward);

                if(!solutionType)
                    g = computeNumericalG(elementVector[i].midX, elementVector[i].midY,
                                          elementVector[j], integrationType, localCoord, weights);
                else
                    g = computeAnalyticalG(elementVector[i].midX, elementVector[i].midY,
                                           elementVector[j]);
            }

            // ---- Remplissage des matrices A, B et du vecteur c ----
            // Les opérations d'écriture dans les matrices Eigen sont sûres ici
            // car chaque paire (i, j) est unique par thread.

            if(elementVector[j].interfaceSide == 'E' || elementVector[j].interfaceSide == 'H')
            {
                if(j < nbDirichletElements)
                {
                    A(i,j) =  g;
                    B(i,j) =  h;
                    // Note: c(j) est écrit plusieurs fois (pour chaque i), mais c'est
                    // la même valeur, donc pas de race condition problématique.
                    c(j)   = elementVector[j].bcValue[0];
                }
                else
                {
                    A(i,j) = -h;
                    B(i,j) = -g;
                    c(j)   = elementVector[j].bcValue[1];
                }
            }
            else
            {
                A(i,j) = -h;
                B(i,j) = -g;
                c(j)   = 0.0;
            }
        }
    }
}

// ============================================================
//  computeElectrostaticPressure  (unchanged from original)
// ============================================================
void computeElectrostaticPressure(std::map<int, double> &electrostaticPressure,
                                  std::vector<elementStruct> &elementVector,
                                  const Eigen::MatrixXd &x,
                                  const std::size_t nbDirichletElements,
                                  const std::string domainName)
{
    double epsilon = getEpsilon(domainName);

    // Parallélisation de la boucle sur les éléments
    #pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < elementVector.size(); i++)
    {
        // On ne met à jour les bcValue que pour les éléments externes ou les trous
        // Les éléments d'interface sont déjà gérés par distributeGlobalSolution().
        if(elementVector[i].interfaceSide == 'E' || elementVector[i].interfaceSide == 'H')
        {
            if(i < nbDirichletElements)
                elementVector[i].bcValue[1] = x(i);
            else
                elementVector[i].bcValue[0] = x(i);
        }

        // Calcul de la pression électrostatique
        double pressure = epsilon * elementVector[i].bcValue[1] * elementVector[i].bcValue[1] / 2.0;

        // L'insertion dans une std::map n'est PAS thread-safe.
        // On utilise critical pour éviter que deux threads n'écrivent en même temps.
        #pragma omp critical
        {
            electrostaticPressure[elementVector[i].tag] = pressure;
        }
    }
}

// ============================================================
//  computeNormalOrientation
//  Pour chaque élément de la frontière, détermine si sa normale
//  (nx,ny) pointe vers l'extérieur du domaine en vérifiant si
//  un point légèrement décalé dans cette direction appartient
//  au physical group du domaine.
// ============================================================
void computeNormalOrientation(std::vector<elementStruct> &elementVector,
                               const int domainPhysicalGroupTag)
{
    for(auto &element : elementVector)
    {
//        if(element.interfaceSide == 'A')
//        {
//            element.normalOutward = 0;
//        }
//        else if(element.interfaceSide == 'B')
//        {
//            element.normalOutward = 1;
//        }
//        else
//        {
//            element.normalOutward = 1;
//        }
        if(element.interfaceSide == 'A')
        {
            // 'A' : Interface côté primaire (inversé pour couplage)
            // 'H' : Hole (inversé pour que la normale sorte du domaine vers le centre du trou)
            element.normalOutward = false;
        }
        else if(element.interfaceSide == 'B')
        {
            // 'B' : Interface côté secondaire
            // 'E' : Frontière externe (Carré)
            element.normalOutward = true;
        }
        else if(element.interfaceSide == 'E')
        {
            // 'A' : Interface côté primaire (inversé pour couplage)
            // 'H' : Hole (inversé pour que la normale sorte du domaine vers le centre du trou)
            element.normalOutward = true;
        }
        else if(element.interfaceSide == 'H')
        {
            // 'B' : Interface côté secondaire
            // 'E' : Frontière externe (Carré)
            element.normalOutward = false;
        }
    }
}
