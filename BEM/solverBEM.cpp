#include "functionsBEM.hpp"

// ============================================================
//  solverBEM  – top-level dispatcher
// ============================================================
void solverBEM(std::map<int, double> &electrostaticPressure,
               int &nbViews,
               std::map<int, std::pair<double,double>> &nodalDisplacements,
               const bool postProcessing,
               const int iteration,
               const bool untangleMesh)
{
    gmsh::vectorpair physicalGroupsDimTags;
    gmsh::model::getPhysicalGroups(physicalGroupsDimTags);

    std::vector<std::string> domainNames;
    std::vector<std::pair<int,int>> domainDimTags;

    for(auto &dimTag : physicalGroupsDimTags)
    {
        std::string name;
        gmsh::model::getPhysicalName(dimTag.first, dimTag.second, name);

        std::stringstream ss(name);
        std::vector<std::string> words;
        std::string word;
        while(getline(ss, word, '_'))
            words.push_back(word);

        if(words.size() == 3 && words[0] == "BEM" && words[1] == "domain")
        {
            domainNames.push_back(words[2]);
            domainDimTags.push_back(dimTag);
        }
    }

    if(domainNames.empty())
    {
        std::cerr << "Error: No BEM domain found!\n"
                  << "  Convention: Physical Surface(\"BEM_domain_<n>\").\n";
        exit(1);
    }

    if(untangleMesh)
    {
        for(std::size_t i = 0; i < domainNames.size(); i++)
        {
            int dim = domainDimTags[i].first;
            int tag = domainDimTags[i].second;
            std::vector<int> tags;
            gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
            gmsh::vectorpair dimTagsBis(tags.size());
            for(std::size_t k = 0; k < tags.size(); k++)
                dimTagsBis[k] = {dim, tags[k]};
            gmsh::model::mesh::optimize("UntangleMeshGeometry", false, 1, dimTagsBis);
        }
    }

    if(iteration == 1)
    {
        std::cout << "\n--------------------[BEM SOLVER]--------------------\n";
        std::cout << "Found " << domainNames.size() << " BEM domain(s).\n";
    }

    bool computePhi           = true;
    bool computeElectricField = true;
    int phiViewTag            = 0;
    int electricFieldViewTag  = 0;
    if(postProcessing)
    {
        if(computePhi)           phiViewTag          = gmsh::view::add("Phi");
        if(computeElectricField) electricFieldViewTag = gmsh::view::add("Electric Field");
    }

    if(domainNames.size() == 1)
    {
        singleDomainSolverBEM(electrostaticPressure, nbViews, nodalDisplacements,
                              postProcessing, iteration,
                              "BEM_domain_" + domainNames[0],
                              phiViewTag, electricFieldViewTag,
                              computePhi, computeElectricField);
    }
    else
    {
        std::vector<domainData> domains(domainNames.size());
        for(std::size_t i = 0; i < domainNames.size(); i++)
        {
            domains[i].name             = domainNames[i];
            domains[i].physicalGroupDim = domainDimTags[i].first;
            domains[i].physicalGroupTag = domainDimTags[i].second;
            domains[i].epsilon          = getEpsilon("BEM_domain_" + domainNames[i]);
        }

        std::vector<interfaceStruct> interfaces;
        detectInterfaces(interfaces, domainNames);

        multiDomainSolverBEM(electrostaticPressure, nbViews, nodalDisplacements,
                             postProcessing, iteration,
                             domains, interfaces,
                             phiViewTag, electricFieldViewTag,
                             computePhi, computeElectricField);
    }

    if(iteration == 1)
        std::cout << "--------------------[BEM SOLVER]--------------------\n\n";

    if(postProcessing)
    {
        if(computePhi)           nbViews++;
        if(computeElectricField) nbViews++;
    }
}

// ============================================================
//  singleDomainSolverBEM  (identical to original)
// ============================================================
void singleDomainSolverBEM(std::map<int, double> &electrostaticPressure,
                            const int nbViews,
                            std::map<int, std::pair<double,double>> &nodalDisplacements,
                            const bool postProcessing,
                            const int iteration,
                            const std::string domainName,
                            const int phiViewTag,
                            const int electricFieldViewTag,
                            const bool computePhi,
                            const bool computeElectricField)
{
    const bool showTime           = true;
    const bool convergence        = false;
    const bool visualisation      = true;
    const bool solutionType       = true;
    const std::string integrationType = "Gauss30";

    double startTotal = omp_get_wtime();
    double start, end;

    std::vector<elementStruct> elementVector;
    int domainPhysicalGroupDim, domainPhysicalGroupTag;

    start = omp_get_wtime();
    std::size_t nbDirichletElements =
        fillElementVector(elementVector, domainPhysicalGroupDim, domainPhysicalGroupTag,
                          nodalDisplacements, domainName);
    end = omp_get_wtime();
    displayTime(start, end, "'fillElementVector'", showTime);

    if(iteration == 1)
        std::cout << "\n- BEM domain '" << domainName << "': "
                  << elementVector.size() << " boundary elements.\n";

    start = omp_get_wtime();
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(elementVector.size(), elementVector.size());
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(elementVector.size(), 1);
    fillSystem(A, B, c, elementVector, nbDirichletElements, integrationType, solutionType);
    end = omp_get_wtime();
    displayTime(start, end, "'fillSystem'", showTime);

    start = omp_get_wtime();
    Eigen::MatrixXd b = B * c;
    Eigen::MatrixXd x = A.partialPivLu().solve(b);
    end = omp_get_wtime();
    displayTime(start, end, "Solve linear system", showTime);

    start = omp_get_wtime();
    computeElectrostaticPressure(electrostaticPressure, elementVector, x,
                                 nbDirichletElements, domainName);
    end = omp_get_wtime();
    displayTime(start, end, "'computeElectrostaticPressure'", showTime);

    if(postProcessing)
    {
        int viewIndex = nbViews;

        std::vector<std::size_t> domainNodeTags;
        std::map<std::size_t, std::size_t> domainNodeIndices;
        std::vector<double> domainNodeCoords;
        std::vector<int> domainEntityTags;
        std::vector<std::vector<int>> domainElementTypes;
        std::vector<std::vector<std::vector<std::size_t>>> domainElementTags;
        std::map<std::size_t, std::vector<std::size_t>> domainElementNodeTags;

        start = omp_get_wtime();
        getModel(domainNodeTags, domainNodeIndices, domainNodeCoords, domainEntityTags,
                 domainElementTypes, domainElementTags, domainElementNodeTags,
                 domainPhysicalGroupDim, domainPhysicalGroupTag);
        end = omp_get_wtime();
        displayTime(start, end, "'getModel'", showTime);

        std::vector<std::vector<std::vector<std::vector<double>>>> elementalPhi;
        std::vector<std::vector<std::vector<std::vector<double>>>> elementNodalPhi;
        std::vector<std::vector<std::vector<std::vector<double>>>> elementalElectricField;

        std::vector<std::string> names;
        gmsh::model::list(names);

        if(visualisation)
        {
            std::map<std::size_t, std::vector<std::size_t>> boundaryNodeToElements;
            computeNodeToElements(boundaryNodeToElements, elementVector);

            start = omp_get_wtime();
            std::map<int, double> internalNodalPhi;
            computeInternalPhi(internalNodalPhi, domainNodeCoords, domainNodeTags,
                               domainNodeIndices, elementVector, integrationType,
                               boundaryNodeToElements, solutionType);
            end = omp_get_wtime();
            displayTime(start, end, "'computeInternalPhi'", showTime);

            start = omp_get_wtime();
            elementNodeData(elementNodalPhi, domainElementTags, domainElementTypes,
                            domainElementNodeTags, elementVector, internalNodalPhi,
                            domainPhysicalGroupTag, boundaryNodeToElements);
            end = omp_get_wtime();
            displayTime(start, end, "'elementNodeData'", showTime);

            start = omp_get_wtime();
            computeElementData(elementalPhi, elementalElectricField, domainElementTags,
                               domainElementTypes, domainEntityTags, elementVector,
                               solutionType, integrationType, false, computeElectricField);
            end = omp_get_wtime();
            displayTime(start, end, "'computeElementData'", showTime);

            if(computePhi)
            {
                displayData(phiViewTag, names[0], "results.msh", "ElementNodeData",
                            domainElementTags, elementNodalPhi, viewIndex);
                viewIndex++;
            }
            if(computeElectricField)
                displayData(electricFieldViewTag, names[0], "results.msh", "ElementData",
                            domainElementTags, elementalElectricField, viewIndex);
        }
        else
        {
            start = omp_get_wtime();
            computeElementData(elementalPhi, elementalElectricField, domainElementTags,
                               domainElementTypes, domainEntityTags, elementVector,
                               solutionType, integrationType, computePhi, computeElectricField);
            end = omp_get_wtime();
            displayTime(start, end, "'computeElementData'", showTime);

            if(computePhi)
            {
                displayData(phiViewTag, names[0], "results.msh", "ElementData",
                            domainElementTags, elementalPhi, viewIndex);
                viewIndex++;
            }
            if(computeElectricField)
                displayData(electricFieldViewTag, names[0], "results.msh", "ElementData",
                            domainElementTags, elementalElectricField, viewIndex);

            if(convergence)
            {
                std::string convModel = "BEM_convergence";
                if(names[0] == convModel)
                    convergenceFun(domainEntityTags, domainElementTypes, domainElementTags,
                                   elementalPhi, elementalElectricField);
                else
                    std::cout << "\tWarning: convergence check requires '"
                              << convModel << ".geo'.\n";
            }
        }
    }

    double endTotal = omp_get_wtime();
    displayTime(startTotal, endTotal, "total BEM time", showTime);
}

// ============================================================
//  assembleGlobalSystem
// ============================================================
void assembleGlobalSystem(Eigen::MatrixXd &A_global,
                          Eigen::MatrixXd &b_global,
                          const std::vector<domainData> &domains,
                          const std::vector<interfaceStruct> &interfaces,
                          const std::string integrationType,
                          const bool solutionType,
                          std::size_t &ifPotentialOffset,
                          std::size_t &ifFluxOffset)
{
    const std::size_t nDomains    = domains.size();
    const std::size_t nInterfaces = interfaces.size();

    std::vector<double> localCoord, weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // 1. CALCUL DES OFFSETS (Séquentiel car rapide et interdépendant)
    std::vector<std::size_t> rowOffsets(nDomains + 1, 0);
    for(std::size_t i = 0; i < nDomains; i++)
        rowOffsets[i+1] = rowOffsets[i] + domains[i].nbRows();

    const std::size_t totalRows = rowOffsets[nDomains];

    std::vector<std::size_t> extColOffsets(nDomains + 1, 0);
    for(std::size_t i = 0; i < nDomains; i++)
        extColOffsets[i+1] = extColOffsets[i] + domains[i].nbExternalElements;

    const std::size_t totalExtCols = extColOffsets[nDomains];

    std::vector<std::size_t> ifColOffsets(nInterfaces + 1, 0);
    for(std::size_t k = 0; k < nInterfaces; k++)
        ifColOffsets[k+1] = ifColOffsets[k] + interfaces[k].nbElements();

    const std::size_t totalIfElements = ifColOffsets[nInterfaces];

    ifPotentialOffset = totalExtCols;
    ifFluxOffset      = totalExtCols + totalIfElements;
    const std::size_t totalCols = totalExtCols + 2 * totalIfElements;

    // Vérification de la matrice carrée
    if(totalRows != totalCols) {
        std::cerr << "BEM ERROR: system is NOT square! (" << totalRows << "x" << totalCols << ")\n";
        exit(1);
    }

    A_global = Eigen::MatrixXd::Zero(totalRows, totalCols);
    b_global = Eigen::MatrixXd::Zero(totalRows, 1);
    std::cout << "  Global system size: " << totalRows << " × " << totalCols << "\n";

    // Helper thread-safe pour trouver l'index du domaine
    auto getDomIdx = [&](const std::string &name) -> std::size_t {
        for(std::size_t i = 0; i < nDomains; i++)
            if(domains[i].name == name) return i;
        return size_t(-1);
    };

    // 2. ASSEMBLAGE PARALLÈLE (Boucle sur les domaines)
    // Chaque thread s'occupe d'un domaine complet (lignes indépendantes)
    #pragma omp parallel for schedule(dynamic)
    for(std::size_t di = 0; di < nDomains; di++)
    {
        const domainData &dom = domains[di];
        const std::size_t N   = dom.elementVector.size();
        const std::size_t nD  = dom.nbDirichletElements;
        const std::size_t nE  = dom.nbExternalElements;
        const std::size_t rOff = rowOffsets[di];
        const std::size_t cOff = extColOffsets[di];

        // Matrices locales pour H et G pour ce domaine
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N, N);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N, N);

        // Calcul des intégrales élémentaires (H et G)
        for(std::size_t i = 0; i < N; i++) {
            for(std::size_t j = 0; j < N; j++) {
                if(i == j) {
                    H(i,j) = -0.5;
                    G(i,j) = (dom.elementVector[i].length / (2.0*M_PI)) *
                             (std::log(dom.elementVector[i].length / 2.0) - 1.0);
                } else {
                    H(i,j) = computeAnalyticalH(dom.elementVector[i].midX,
                                                dom.elementVector[i].midY,
                                                dom.elementVector[j],
                                                dom.elementVector[j].normalOutward);
                    if(solutionType)
                        G(i,j) = computeAnalyticalG(dom.elementVector[i].midX,
                                                    dom.elementVector[i].midY,
                                                    dom.elementVector[j]);
                    else
                        G(i,j) = computeNumericalG(dom.elementVector[i].midX,
                                                   dom.elementVector[i].midY,
                                                   dom.elementVector[j], integrationType, localCoord, weights);
                }
            }
        }

        // Remplissage des colonnes externes (Dirichlet / Neumann)
        for(std::size_t i = 0; i < N; i++) {
            for(std::size_t j = 0; j < nE; j++) {
                if(j < nD) {
                    A_global(rOff + i, cOff + j) += G(i, j);
                    b_global(rOff + i)           += H(i, j) * dom.elementVector[j].bcValue[0];
                } else {
                    A_global(rOff + i, cOff + j) += -H(i, j);
                    b_global(rOff + i)           += -G(i, j) * dom.elementVector[j].bcValue[1];
                }
            }
        }

        // Remplissage des blocs d'interface
        for(std::size_t k = 0; k < nInterfaces; k++) {
            const interfaceStruct &iface = interfaces[k];
            const std::size_t idxA = getDomIdx(iface.domainA);
            const std::size_t idxB = getDomIdx(iface.domainB);

            const bool isDomainA = (di == idxA);
            const bool isDomainB = (di == idxB);
            if(!isDomainA && !isDomainB) continue;

            const std::vector<std::size_t> &ifElemIdx = isDomainA ? iface.elemIndicesA : iface.elemIndicesB;
            const std::size_t nIf         = ifElemIdx.size();
            const std::size_t potColBase  = ifPotentialOffset + ifColOffsets[k];
            const std::size_t fluxColBase = ifFluxOffset      + ifColOffsets[k];

            for(std::size_t i = 0; i < N; i++) {
                for(std::size_t m = 0; m < nIf; m++) {
                    const std::size_t j = ifElemIdx[m];
                    A_global(rOff + i, potColBase + m) += -H(i, j);

                    if(isDomainA)
                        A_global(rOff + i, fluxColBase + m) += G(i, j);
                    else
                        A_global(rOff + i, fluxColBase + m) += -G(i, j) * iface.kRatio;
                }
            }
        }
    } // Fin du pragma omp parallel for

    // 3. DEBUG & ANALYSE (Séquentiel)
    std::cout << "--- Matrix Layout Debug ---\n";
    std::cout << "Total Rows: " << totalRows << " | Total Cols: " << totalCols << "\n";

    // Test de singularité (Peut être lent sur de très grandes matrices)
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A_global);
    if(!lu.isInvertible()) {
        std::cerr << "WARNING: Global Matrix is nearly singular! (Det: " << A_global.determinant() << ")\n";
    }
}

// ============================================================
//  distributeGlobalSolution
// ============================================================
void distributeGlobalSolution(std::vector<domainData> &domains,
                              const std::vector<interfaceStruct> &interfaces,
                              const Eigen::MatrixXd &X_global,
                              const std::size_t ifPotentialOffset,
                              const std::size_t ifFluxOffset)
{
    const std::size_t nDomains    = domains.size();
    const std::size_t nInterfaces = interfaces.size();

    // 1. Recalcul des offsets (Séquentiel)
    std::vector<std::size_t> extColOffsets(nDomains + 1, 0);
    for(std::size_t i = 0; i < nDomains; i++)
        extColOffsets[i+1] = extColOffsets[i] + domains[i].nbExternalElements;

    std::vector<std::size_t> ifColOffsets(nInterfaces + 1, 0);
    for(std::size_t k = 0; k < nInterfaces; k++)
        ifColOffsets[k+1] = ifColOffsets[k] + interfaces[k].nbElements();

    // Helper thread-safe
    auto getDomIdx = [&](const std::string &name) -> std::size_t {
        for(std::size_t i = 0; i < nDomains; i++)
            if(domains[i].name == name) return i;
        return 0;
    };

    // 2. Éléments externes (Parallélisé par domaine)
    #pragma omp parallel for schedule(static)
    for(std::size_t di = 0; di < nDomains; di++)
    {
        domainData &dom = domains[di];
        const std::size_t cOff = extColOffsets[di];
        const std::size_t nD   = dom.nbDirichletElements;
        const std::size_t nE   = dom.nbExternalElements;

        for(std::size_t j = 0; j < nE; j++)
        {
            double val = X_global(cOff + j);
            if(j < nD)
                dom.elementVector[j].bcValue[1] = val; // Dirichlet -> qn inconnu
            else
                dom.elementVector[j].bcValue[0] = val; // Neumann -> phi inconnu
        }
    }

    // 3. Éléments d'interface (Parallélisé par interface)
    #pragma omp parallel for schedule(static)
    for(std::size_t k = 0; k < nInterfaces; k++)
    {
        const interfaceStruct &iface    = interfaces[k];
        const std::size_t      nIf      = iface.nbElements();
        const std::size_t      potBase  = ifPotentialOffset + ifColOffsets[k];
        const std::size_t      fluxBase = ifFluxOffset      + ifColOffsets[k];

        const std::size_t idxA = getDomIdx(iface.domainA);
        const std::size_t idxB = getDomIdx(iface.domainB);

        for(std::size_t m = 0; m < nIf; m++)
        {
            const double u_if  = X_global(potBase  + m);
            const double qn_if = X_global(fluxBase + m);

            const std::size_t jA = iface.elemIndicesA[m];
            const std::size_t jB = iface.elemIndicesB[m];

            // Pas de collision ici car jA et jB sont uniques à chaque interface/domaine
            domains[idxA].elementVector[jA].bcValue[0] = u_if;
            domains[idxA].elementVector[jA].bcValue[1] = qn_if;

            domains[idxB].elementVector[jB].bcValue[0] = u_if;
            domains[idxB].elementVector[jB].bcValue[1] = -qn_if; // Flux opposé sur l'autre face
        }

        // On garde le debug pour le premier élément (exécuté par un thread à la fois)
        #pragma omp critical(debug_output)
        {
            std::cout << "[DEBUG DISTRIBUTE] Interface: " << iface.name << " (Thread " << omp_get_thread_num() << ")\n"
                      << "  -> X_global (u, qn) m=0: " << X_global(potBase) << ", " << X_global(fluxBase) << std::endl;
        }
    }
}

// ============================================================
//  multiDomainSolverBEM
// ============================================================
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
                          const bool computeElectricField)
{
    const bool showTime           = true;
    const bool solutionType       = true;
    const std::string integrationType = "Gauss30";

    double startTotal = omp_get_wtime();
    double start, end;

    // ---- Étape 1 : remplissage des éléments externes (Parallèle par domaine) ----
    start = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic)
    for(std::size_t d = 0; d < domains.size(); d++)
    {
        auto &dom = domains[d];
        dom.nbDirichletElements =
            fillElementVector(dom.elementVector,
                              dom.physicalGroupDim,
                              dom.physicalGroupTag,
                              nodalDisplacements,
                              "BEM_domain_" + dom.name);
        dom.nbExternalElements = dom.elementVector.size();
    }
    end = omp_get_wtime();
    displayTime(start, end, "fillElementVector (all domains)", showTime);

    // ---- Étape 2 : éléments d'interface (Séquentiel car nécessite domainByName) ----
    start = omp_get_wtime();
    auto domainByName = [&](const std::string &n) -> domainData& {
        for(auto &d : domains) if(d.name == n) return d;
        std::cerr << "Domain '" << n << "' not found.\n"; exit(1);
    };

    for(auto &iface : interfaces)
    {
        domainData &dA = domainByName(iface.domainA);
        domainData &dB = domainByName(iface.domainB);
        addInterfaceElements(dA, dB, iface, nodalDisplacements);
    }
    end = omp_get_wtime();
    displayTime(start, end, "addInterfaceElements (all interfaces)", showTime);

    // ---- Étape 3 : assemblage global (Déjà parallélisé en interne) ----
    start = omp_get_wtime();
    Eigen::MatrixXd A_global, b_global;
    std::size_t ifPotentialOffset, ifFluxOffset;
    assembleGlobalSystem(A_global, b_global, domains, interfaces,
                         integrationType, solutionType,
                         ifPotentialOffset, ifFluxOffset);
    end = omp_get_wtime();
    displayTime(start, end, "assembleGlobalSystem", showTime);

    // ---- Étape 4 : résolution (Eigen utilise OpenMP si lié avec BLAS) ----
    start = omp_get_wtime();
    // partialPivLu est robuste, mais pour de très gros systèmes, préférez un solveur itératif
    Eigen::MatrixXd X_global = A_global.partialPivLu().solve(b_global);
    end = omp_get_wtime();
    displayTime(start, end, "Solve global system", showTime);

    // ---- Étape 5 : distribution de la solution (Déjà parallélisé en interne) ----
    start = omp_get_wtime();
    distributeGlobalSolution(domains, interfaces, X_global, ifPotentialOffset, ifFluxOffset);
    end = omp_get_wtime();
    displayTime(start, end, "distributeGlobalSolution", showTime);

    // ---- Étape 6 : pression électrostatique (Parallèle) ----
    // Step 6a: store per-domain pressure (useful for visualisation / debugging)
    computeElectrostaticPressureMultiDomain(electrostaticPressure, domains);

    // Step 6b: overwrite with the *jump* pressure on every interface that
    //          is also a BEM_FEM_boundary.  The FEM solver's applyElecPressure()
    //          looks up element tags that belong to BEM_FEM_boundary, which in
    //          the .geo is the same set of curves as BEM_interface_2_3
    //          (i.e. the droplet surface).  We therefore only call the jump
    //          function for interfaces whose name contains "interface_2_3"
    //          (or whichever interface is shared with the FEM domain).
    for(const auto &iface : interfaces)
        {
            // On vérifie si le nom de l'interface contient le nom du groupe de couplage
            // ou si l'interface lie un domaine FEM.
            if(iface.domainA == "FEM_domain" || iface.domainB == "FEM_domain")
            {
                computeElectrostaticPressureJump(electrostaticPressure, domains, iface);
            }
        }

    // ---- Étape 7 : post-processing (Parallèle par domaine) ----
    if(postProcessing)
    {
        start = omp_get_wtime();
        int viewIndexBase = nbViews;
        std::vector<std::string> modelNames;
        gmsh::model::list(modelNames);
        std::string modelName = modelNames[0];

        // On parallélise le calcul pour chaque domaine
        #pragma omp parallel for schedule(dynamic)
        for(std::size_t d = 0; d < domains.size(); d++)
        {
            auto &dom = domains[d];

            // Variables locales au thread (très important !)
            std::vector<std::size_t> dNodeTags;
            std::map<std::size_t, std::size_t> dNodeIndices;
            std::vector<double> dNodeCoords;
            std::vector<int> dEntityTags;
            std::vector<std::vector<int>> dElemTypes;
            std::vector<std::vector<std::vector<std::size_t>>> dElemTags;
            std::map<std::size_t, std::vector<std::size_t>> dElemNodeTags;

            // Chargement de la géométrie locale
            getModel(dNodeTags, dNodeIndices, dNodeCoords, dEntityTags,
                     dElemTypes, dElemTags, dElemNodeTags,
                     dom.physicalGroupDim, dom.physicalGroupTag);

            // Ces fonctions utilisent déjà l'optimisation par unordered_map
            std::map<std::size_t, std::vector<std::size_t>> bNodeToElems;
            computeNodeToElements(bNodeToElems, dom.elementVector);

            std::map<int, double> internalPhi;
            computeInternalPhi(internalPhi, dNodeCoords, dNodeTags,
                               dNodeIndices, dom.elementVector, integrationType,
                               bNodeToElems, solutionType);

//            for(const auto &iface : interfaces) {
//                // Si ce domaine (dom) est le porteur (Side A) de cette interface
//                if(iface.domainA == dom.name && iface.name.find("_inc") != std::string::npos) {
//
//                    Force2D F = computeDEPForceMST(dom.elementVector, internalPhi, iface.name, epsilon_oil);
//
//                    #pragma omp critical(print_force)
//                    {
//                        std::cout << "[Force DEP] Interface: " << iface.name
//                                  << " | Fx: " << F.fx << " N/m"
//                                  << " | Fy: " << F.fy << " N/m" << std::endl;
//                    }
//                }
//            }

            // Calcul des données de visualisation
            std::vector<std::vector<std::vector<std::vector<double>>>> nPhi, ePhi, eField;

            elementNodeData(nPhi, dElemTags, dElemTypes, dElemNodeTags,
                            dom.elementVector, internalPhi, dom.physicalGroupTag, bNodeToElems);

            computeElementData(ePhi, eField, dElemTags, dElemTypes, dEntityTags,
                               dom.elementVector, solutionType, integrationType, false, computeElectricField);

            // Section critique pour les appels à l'API GMSH (pas thread-safe)
            #pragma omp critical(gmsh_api)
            {
                if(computePhi)
                    displayData(phiViewTag, modelName, "results.msh", "ElementNodeData",
                                dElemTags, nPhi, viewIndexBase);

                if(computeElectricField)
                    displayData(electricFieldViewTag, modelName, "results.msh", "ElementData",
                                dElemTags, eField, viewIndexBase + (computePhi ? 1 : 0));
            }
        }
        end = omp_get_wtime();
        displayTime(start, end, "Post-processing (parallel domains)", showTime);
    }

    double endTotal = omp_get_wtime();
    displayTime(startTotal, endTotal, "total multi-domain BEM time", showTime);
}
