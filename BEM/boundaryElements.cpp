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
std::size_t fillElementVector(
    std::vector<elementStruct>& elementVector,
    int& domainPhysicalGroupDim,
    int& domainPhysicalGroupTag,
    std::map<int, std::pair<double, double>>& nodalDisplacements,
    const std::string domainName)
{
    std::size_t nbDirichletElements = 0;

    std::map<std::string, std::pair<bool, double>> boundaryConditions;
    getBC(boundaryConditions, domainName);

    // ---------------------------------------------------------
    // 1. Identifier le domaine courant et ses surfaces 2D
    // ---------------------------------------------------------
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    bool domainFound = false;
    std::vector<int> domainSurfaceTags;

    for (const auto& [dim, tag] : physicalGroupsDimTags)
    {
        if (dim != 2)
            continue;

        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);

        if (name != domainName)
            continue;

        domainPhysicalGroupDim = dim;
        domainPhysicalGroupTag = tag;

        gmsh::model::getEntitiesForPhysicalGroup(
            dim, tag, domainSurfaceTags);

        domainFound = true;
        break;
    }

    if (!domainFound)
    {
        std::cerr << "[ERROR] Domain not found: "
                  << domainName << '\n';
        std::exit(EXIT_FAILURE);
    }

    // ---------------------------------------------------------
    // 2. Parcourir les groupes physiques portant une BC
    // ---------------------------------------------------------
    for (const auto& [physicalGroupDim, physicalGroupTag]
         : physicalGroupsDimTags)
    {
        std::string physicalGroupName;
        gmsh::model::getPhysicalName(
            physicalGroupDim, physicalGroupTag, physicalGroupName);

        auto bcIt = boundaryConditions.find(physicalGroupName);
        if (bcIt == boundaryConditions.end())
            continue;

        // Les interfaces sont ajoutées séparément par
        // addInterfaceElements().
        if (physicalGroupName.rfind("BEM_interface_", 0) == 0)
            continue;

        // Les frontières externes doivent être des courbes 1D.
        if (physicalGroupDim != 1)
            continue;

        const bool isHole =
            physicalGroupName.find("Hole") != std::string::npos ||
            physicalGroupName.find("hole") != std::string::npos;

        const bool bcType = bcIt->second.first;
        const double bcValue = bcIt->second.second;

        std::vector<int> entityTags;
        gmsh::model::getEntitiesForPhysicalGroup(
            physicalGroupDim, physicalGroupTag, entityTags);

        // -----------------------------------------------------
        // 3. Traiter chaque courbe du groupe physique
        // -----------------------------------------------------
        for (const int curveTag : entityTags)
        {
            // Vérifier que la courbe est adjacente à une surface
            // du domaine courant.
            std::vector<int> upward, downward;
            gmsh::model::getAdjacencies(
                1, curveTag, upward, downward);

            const bool belongsToDomain =
                std::any_of(
                    upward.begin(), upward.end(),
                    [&](const int surfaceTag)
                    {
                        return std::find(
                            domainSurfaceTags.begin(),
                            domainSurfaceTags.end(),
                            surfaceTag
                        ) != domainSurfaceTags.end();
                    });

            if (!belongsToDomain)
                continue;

            // Éléments linéaires et quadratiques (ordre 2).
            const std::vector<int> elementTypes = {1, 8};

            for (const int elementType : elementTypes)
            {
                const int nbNodesPerElement =
                    (elementType == 8) ? 3 : 2;

                std::vector<std::size_t> elementTags;
                std::vector<std::size_t> elementNodeTags;

                gmsh::model::mesh::getElementsByType(
                    elementType, elementTags,
                    elementNodeTags, curveTag);

                if (elementTags.empty())
                    continue;

                std::vector<std::size_t> nodeTags;
                std::vector<double> nodeCoords, parametricNodeCoords;

                gmsh::model::mesh::getNodesByElementType(
                    elementType, nodeTags, nodeCoords,
                    parametricNodeCoords, curveTag, false);

                std::map<std::size_t, std::size_t> nodeIndices;
                for (std::size_t k = 0; k < nodeTags.size(); ++k)
                    nodeIndices[nodeTags[k]] = k;

                // -------------------------------------------------
                // 4. Construction parallèle des éléments BIE
                // -------------------------------------------------
                #pragma omp parallel
                {
                    std::vector<elementStruct> localDirichlet;
                    std::vector<elementStruct> localNeumann;

                    #pragma omp for schedule(guided)
                    for (std::size_t el = 0;
                         el < elementTags.size(); ++el)
                    {
                        elementStruct element;

                        element.tag = elementTags[el];
                        element.domainPhysicalGroupTag =
                            domainPhysicalGroupTag;

                        element.nodeTags[0] =
                            elementNodeTags[nbNodesPerElement * el];

                        element.nodeTags[1] =
                            elementNodeTags[nbNodesPerElement * el + 1];

                        element.midNodeTag =
                            (nbNodesPerElement == 3)
                            ? elementNodeTags[nbNodesPerElement * el + 2]
                            : 0;

                        element.interfaceSide = isHole ? 'H' : 'E';
                        element.interfaceName = physicalGroupName;
                        element.interfaceIndex = -1;

                        std::pair<double, double> d1 = {0.0, 0.0};
                        std::pair<double, double> d2 = {0.0, 0.0};

                        auto it1 = nodalDisplacements.find(
                            static_cast<int>(element.nodeTags[0]));
                        if (it1 != nodalDisplacements.end())
                            d1 = it1->second;

                        auto it2 = nodalDisplacements.find(
                            static_cast<int>(element.nodeTags[1]));
                        if (it2 != nodalDisplacements.end())
                            d2 = it2->second;

                        computeGeometricParam(
                            element, nodeIndices, nodeCoords, d1, d2);

                        if (bcType)
                        {
                            element.bcValue[0] = bcValue;
                            localDirichlet.push_back(element);
                        }
                        else
                        {
                            element.bcValue[1] = bcValue;
                            localNeumann.push_back(element);
                        }
                    }

                    #pragma omp critical
                    {
                        if (!localDirichlet.empty())
                        {
                            elementVector.insert(
                                elementVector.begin(),
                                localDirichlet.begin(),
                                localDirichlet.end());

                            nbDirichletElements += localDirichlet.size();
                        }

                        if (!localNeumann.empty())
                        {
                            elementVector.insert(
                                elementVector.end(),
                                localNeumann.begin(),
                                localNeumann.end());
                        }
                    }
                }
            }
        }
    }

    // Calcul de l'orientation des normales après construction.
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

        std::stringstream ss(name);
        std::vector<std::string> words;
        std::string word;
        while(getline(ss, word, '_'))
            words.push_back(word);

        // On attend au moins BEM_interface_A_B (4 mots)
        if(words.size() < 4 || words[0] != "BEM" || words[1] != "interface")
            continue;

        interfaceStruct iface;
        iface.name = name; // On garde le nom complet (ex: BEM_interface_1_2_jux)

        // On récupère dA et dB
        iface.domainA = words[2];
        iface.domainB = words[3];

        // On vérifie si les noms de domaines sont valides
        auto itA = std::find(domainNames.begin(), domainNames.end(), iface.domainA);
        auto itB = std::find(domainNames.begin(), domainNames.end(), iface.domainB);

        if(itA == domainNames.end() || itB == domainNames.end()) {
            std::cout << "\tWarning: interface '" << name << "' inconnue. Skipping.\n";
            continue;
        }

        // Calcul du k-ratio
        double epsA = getEpsilon("BEM_domain_" + iface.domainA);
        double epsB = getEpsilon("BEM_domain_" + iface.domainB);
        iface.kRatio = (epsB > 0.0) ? (epsA / epsB) : 1.0;

        interfaces.push_back(iface);
        std::cout << "  Detected interface '" << name << "' (type: "
                  << (name.find("_inc") != std::string::npos ? "Inclusion" : "Juxtaposition")
                  << ", k_ratio = " << iface.kRatio << ")\n";
    }
}


void addInterfaceElements(
    domainData& domainA,
    domainData& domainB,
    interfaceStruct& iface,
    std::map<int, std::pair<double, double>>& nodalDisplacements)
{
    iface.elemIndicesA.clear();
    iface.elemIndicesB.clear();
    // ---------------------------------------------------------
    // 1. Trouver le groupe physique de l'interface
    // ---------------------------------------------------------
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    int ifaceDim = -1;
    int ifaceTag = -1;

    for (const auto& [dim, tag] : physicalGroupsDimTags)
    {
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);

        if (name == iface.name)
        {
            ifaceDim = dim;
            ifaceTag = tag;
            break;
        }
    }

    if (ifaceDim != 1 || ifaceTag < 0)
    {
        std::cerr
            << "[ERROR] Cannot find Physical Curve '"
            << iface.name << "' (dimension 1).\n";
        return;
    }

    // ---------------------------------------------------------
    // 2. Récupérer les courbes de l'interface
    // ---------------------------------------------------------
    std::vector<int> entityTags;
    gmsh::model::getEntitiesForPhysicalGroup(
        ifaceDim, ifaceTag, entityTags);

    const std::vector<int> elementTypes = {1, 8};

    // ---------------------------------------------------------
    // 3. Construire les deux représentations BIE
    // ---------------------------------------------------------
    for (const int entityTag : entityTags)
    {
        for (const int elType : elementTypes)
        {
            const int nbNPE = (elType == 8) ? 3 : 2;

            std::vector<std::size_t> elTags;
            std::vector<std::size_t> elNodeTags;

            gmsh::model::mesh::getElementsByType(
                elType, elTags, elNodeTags, entityTag);

            if (elTags.empty())
                continue;

            std::vector<std::size_t> nodeTags;
            std::vector<double> nodeCoords, parametric;

            gmsh::model::mesh::getNodesByElementType(
                elType, nodeTags, nodeCoords,
                parametric, entityTag, false);

            std::map<std::size_t, std::size_t> nodeIndices;

            for (std::size_t k = 0; k < nodeTags.size(); ++k)
                nodeIndices[nodeTags[k]] = k;

            for (std::size_t el = 0; el < elTags.size(); ++el)
            {
                elementStruct base;

                base.tag = elTags[el];

                base.nodeTags[0] =
                    elNodeTags[nbNPE * el];

                base.nodeTags[1] =
                    elNodeTags[nbNPE * el + 1];

                base.midNodeTag =
                    (nbNPE == 3)
                    ? elNodeTags[nbNPE * el + 2]
                    : 0;

                base.bcValue[0] = 0.0;
                base.bcValue[1] = 0.0;
                base.interfaceName = iface.name;
                base.interfaceIndex = -1;

                // Les déplacements sont consultés en lecture seule.
                auto it1 = nodalDisplacements.find(
                    static_cast<int>(base.nodeTags[0]));

                auto it2 = nodalDisplacements.find(
                    static_cast<int>(base.nodeTags[1]));

                std::pair<double, double> d1 = {0.0, 0.0};
                std::pair<double, double> d2 = {0.0, 0.0};

                if (it1 != nodalDisplacements.end())
                    d1 = it1->second;

                if (it2 != nodalDisplacements.end())
                    d2 = it2->second;

                computeGeometricParam(
                    base, nodeIndices, nodeCoords, d1, d2);

                // ---------------------------------------------
                // Côté A : représentation dans le domaine A
                // ---------------------------------------------
                elementStruct elA = base;

                elA.domainPhysicalGroupTag =
                    domainA.physicalGroupTag;

                elA.interfaceSide = 'A';

                iface.elemIndicesA.push_back(
                    domainA.elementVector.size());

                domainA.elementVector.push_back(elA);

                // ---------------------------------------------
                // Côté B : représentation dans le domaine B
                // ---------------------------------------------
                elementStruct elB = base;

                elB.domainPhysicalGroupTag =
                    domainB.physicalGroupTag;

                elB.interfaceSide = 'B';

                iface.elemIndicesB.push_back(
                    domainB.elementVector.size());

                domainB.elementVector.push_back(elB);
            }
        }
    }

    // ---------------------------------------------------------
    // 4. Réorienter les normales après ajout de l'interface
    // ---------------------------------------------------------
    computeNormalOrientation(
        domainA.elementVector,
        domainA.physicalGroupTag);

    computeNormalOrientation(
        domainB.elementVector,
        domainB.physicalGroupTag);

    std::cout
        << "  Interface '" << iface.name << "': added "
        << iface.elemIndicesA.size()
        << " element(s) per side.\n";
}

void computeNormalOrientation(std::vector<elementStruct> &elementVector,
                               const int domainPhysicalGroupTag)
{
    for(auto &element : elementVector)
    {
        // On cherche le flag dans le nom de l'interface
        bool isInclusion = (element.interfaceName.find("_inc") != std::string::npos);

        // --- Cas de la frontière externe 'E' ou des trous 'H' ---
        if(element.interfaceSide == 'E')
        {
            element.normalOutward = true;
            continue;
        }
        if(element.interfaceSide == 'H')
        {
            element.normalOutward = false;
            continue;
        }

        // --- Cas des interfaces 'A' et 'B' ---
        if(isInclusion)
        {
            // LOGIQUE GOUTTE (Domaine A entoure Domaine B)
            if(element.interfaceSide == 'A')      element.normalOutward = false;
            else if(element.interfaceSide == 'B') element.normalOutward = true;
        }
        else
        {
            // LOGIQUE JUXTAPOSITION (Par défaut ou si _jux est présent)
            // Ici, pour le Domaine A, la normale Gmsh est déjà sortante
            if(element.interfaceSide == 'A')      element.normalOutward = true;
            else if(element.interfaceSide == 'B') element.normalOutward = false;
        }
    }
}

// ============================================================
//  MULTI-DOMAIN  — call this AFTER distributeGlobalSolution()
//
//  Iterates over all domains and stores the pressure for every
//  element.  For external / Dirichlet / Neumann elements this is
//  simply (eps/2)*phi_n^2.  The physically meaningful quantity
//  for the FEM coupling is the *jump* on shared interfaces; call
//  computeElectrostaticPressureJump() instead of (or in addition
//  to) this function when you need that.
// ============================================================
void computeElectrostaticPressureMultiDomain(
        std::map<int, double>            &electrostaticPressure,
        const std::vector<domainData>    &domains)
{
    for(const auto &dom : domains)
    {
        const double eps = dom.epsilon;

        #pragma omp parallel for schedule(guided)
        for(std::size_t i = 0; i < dom.elementVector.size(); i++)
        {
            const double phi_n    = dom.elementVector[i].bcValue[1];
            const double pressure = eps * phi_n * phi_n * 0.5;

            #pragma omp critical
            {
                electrostaticPressure[dom.elementVector[i].tag] = pressure;
            }
        }
    }
}

void computeElectrostaticPressureJump(
        std::map<int, double>         &electrostaticPressure,
        const std::vector<domainData> &domains,
        const interfaceStruct         &iface)
{
    // Find domain A and domain B
    const domainData *pA = nullptr, *pB = nullptr;
    for(const auto &d : domains)
    {
        if(d.name == iface.domainA) pA = &d;
        if(d.name == iface.domainB) pB = &d;
    }
    if(!pA || !pB)
    {
        std::cerr << "[computeElectrostaticPressureJump] domain not found for interface "
                  << iface.name << "\n";
        return;
    }

    const double epsA = pA->epsilon;
    const double epsB = pB->epsilon;
    const std::size_t nIf = iface.nbElements();

    #pragma omp parallel for schedule(guided)
    for(std::size_t m = 0; m < nIf; m++)
    {
        const std::size_t jA = iface.elemIndicesA[m];
        const std::size_t jB = iface.elemIndicesB[m];

        // Net Maxwell traction pointing from B into A (outward of B = droplet interior)
        const auto& elA = pA->elementVector.at(jA);
        const auto& elB = pB->elementVector.at(jB);

        // Normale commune dirigée de B vers A
        const double nx = elB.normalOutward ? elB.nx : -elB.nx;
        const double ny = elB.normalOutward ? elB.ny : -elB.ny;

        // Contrainte normale de Maxwell des deux côtés
        const double stressA = maxwellNormalStress(
            elA, epsA, nx, ny);

        const double stressB = maxwellNormalStress(
            elB, epsB, nx, ny);

        // Saut de contrainte normale [Pa]
        const double p_jump = stressA - stressB;

        #pragma omp critical
        {
            electrostaticPressure[elA.tag] = p_jump;
        }
    }
}
