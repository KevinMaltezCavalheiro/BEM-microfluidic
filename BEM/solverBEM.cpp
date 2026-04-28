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
    const std::string integrationType = "Gauss2";

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

    // ---- Offsets des lignes (une ligne par point de collocation = un élément) ----
    // Chaque domaine contribue nbRows() = nbExternalElements + nbInterfaceElements lignes.
    std::vector<std::size_t> rowOffsets(nDomains + 1, 0);
    for(std::size_t i = 0; i < nDomains; i++)
        rowOffsets[i+1] = rowOffsets[i] + domains[i].nbRows();
    const std::size_t totalRows = rowOffsets[nDomains];

    // ---- Offsets des colonnes pour les inconnues externes ----
    // Domaine i : d'abord ses nbDirichletElements colonnes (qn_j inconnu),
    //             puis ses (nbExternalElements - nbDirichletElements) colonnes (phi_j inconnu).
    std::vector<std::size_t> extColOffsets(nDomains + 1, 0);
    for(std::size_t i = 0; i < nDomains; i++)
        extColOffsets[i+1] = extColOffsets[i] + domains[i].nbExternalElements;
    const std::size_t totalExtCols = extColOffsets[nDomains];

    // ---- Offsets des colonnes pour les inconnues d'interface ----
    std::vector<std::size_t> ifColOffsets(nInterfaces + 1, 0);
    for(std::size_t k = 0; k < nInterfaces; k++)
        ifColOffsets[k+1] = ifColOffsets[k] + interfaces[k].nbElements();
    const std::size_t totalIfElements = ifColOffsets[nInterfaces];

    ifPotentialOffset = totalExtCols;
    ifFluxOffset      = totalExtCols + totalIfElements;
    const std::size_t totalCols = totalExtCols + 2 * totalIfElements;

    if(totalRows != totalCols)
    {
        std::cerr << "\nBEM assembler ERROR: system is NOT square!\n"
                  << "  totalRows = " << totalRows << "\n"
                  << "  totalCols = " << totalCols << "\n"
                  << "  totalExtCols     = " << totalExtCols << "\n"
                  << "  totalIfElements  = " << totalIfElements << "\n"
                  << "  (each interface contributes N_if potentials + N_if fluxes = 2*N_if columns)\n"
                  << "  Hint: check that nbExternalElements is set correctly for each domain,\n"
                  << "  and that the interface mesh is conformal (same nodes on both sides).\n";
        exit(1);
    }

    A_global = Eigen::MatrixXd::Zero(totalRows, totalCols);
    b_global = Eigen::MatrixXd::Zero(totalRows, 1);

    std::cout << "  Global system size: " << totalRows << " × " << totalCols << "\n";

    // ---- Helper : trouve l'index du domaine par son nom ----
    auto domIdx = [&](const std::string &name) -> std::size_t {
        for(std::size_t i = 0; i < nDomains; i++)
            if(domains[i].name == name) return i;
        std::cerr << "assembleGlobalSystem: domain '" << name << "' not found.\n"; exit(1);
    };

    // ---- Boucle sur les domaines ----
    for(std::size_t di = 0; di < nDomains; di++)
    {
        const domainData   &dom  = domains[di];
        const std::size_t   N    = dom.elementVector.size();   // nb éléments totaux du domaine
        const std::size_t   nD   = dom.nbDirichletElements;    // nb éléments Dirichlet externes
        const std::size_t   nE   = dom.nbExternalElements;     // nb éléments externes totaux
        const std::size_t   rOff = rowOffsets[di];             // première ligne dans A_global
        const std::size_t   cOff = extColOffsets[di];          // première colonne externe dans A_global

        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N, N);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N, N);

        //#pragma omp parallel for collapse(2) schedule(guided)
        for(std::size_t i = 0; i < N; i++)
        {
            for(std::size_t j = 0; j < N; j++)
            {
                if(i == j)
                {
                    H(i,j) = -0.5;
                    G(i,j) = (dom.elementVector[i].length / (2.0*M_PI)) *
                              (std::log(dom.elementVector[i].length / 2.0) - 1.0);
                }
                else
                {
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
                                                   dom.elementVector[j], integrationType);
                }
            }
        }

        // =====================================================================
        // COLONNE EXTERNE : éléments j = 0 .. nE-1  (Dirichlet puis Neumann)
        // =====================================================================
        for(std::size_t i = 0; i < N; i++)
        {
            for(std::size_t j = 0; j < nE; j++)
            {
                if(j < nD)
                {
                    // Dirichlet : qn_j inconnu, phi_j connu
                    A_global(rOff + i, cOff + j) += G(i, j);
                    b_global(rOff + i)            += H(i, j) * dom.elementVector[j].bcValue[0];
                }
                else
                {
                    // Neumann : phi_j inconnu, qn_j connu
                    A_global(rOff + i, cOff + j) += -H(i, j);
                    b_global(rOff + i)            += -G(i, j) * dom.elementVector[j].bcValue[1];
                }
            }
        }

        // =====================================================================
        // COLONNES D'INTERFACE : u_if et qn_if
        // =====================================================================
        for(std::size_t k = 0; k < nInterfaces; k++)
        {
            const interfaceStruct &iface = interfaces[k];
            const std::size_t idxA = domIdx(iface.domainA);
            const std::size_t idxB = domIdx(iface.domainB);

            const bool isDomainA = (di == idxA);
            const bool isDomainB = (di == idxB);
            if(!isDomainA && !isDomainB) continue;

            const std::vector<std::size_t> &ifElemIdx = isDomainA ?
                                                         iface.elemIndicesA :
                                                         iface.elemIndicesB;
            const std::size_t nIf       = ifElemIdx.size();
            const std::size_t potColBase  = ifPotentialOffset + ifColOffsets[k];
            const std::size_t fluxColBase = ifFluxOffset      + ifColOffsets[k];

            for(std::size_t i = 0; i < N; i++)
            {
                for(std::size_t m = 0; m < nIf; m++)
                {
                    const std::size_t j = ifElemIdx[m];

                    A_global(rOff + i, potColBase + m) += -H(i, j);

                    if(isDomainA)
                        A_global(rOff + i, fluxColBase + m) += G(i, j);
                    else
                        A_global(rOff + i, fluxColBase + m) += -G(i, j) * iface.kRatio;
                }
            }
        }
    } // fin boucle domaines
    std::cout << "--- Matrix Layout Debug ---" << std::endl;
    std::cout << "Total Rows: " << totalRows << " | Total Cols: " << totalCols << std::endl;
    std::cout << "External Columns (Total): " << totalExtCols << std::endl;
    for(size_t i=0; i < nDomains; ++i)
        std::cout << "  Dom " << i << " (" << domains[i].name << ") starts at row " << rowOffsets[i] << " and extCol " << extColOffsets[i] << std::endl;

    std::cout << "Interface Potentials start at: " << ifPotentialOffset << std::endl;
    std::cout << "Interface Fluxes start at: " << ifFluxOffset << std::endl;
    for(size_t k=0; k < nInterfaces; ++k)
        std::cout << "  Interface " << k << " (" << interfaces[k].name << ") offset: " << ifColOffsets[k] << " (size: " << interfaces[k].nbElements() << ")" << std::endl;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A_global);
    if(!lu.isInvertible()) {
        std::cerr << "WARNING: Global Matrix is nearly singular!" << std::endl;
        std::cerr << "Determinant: " << A_global.determinant() << std::endl;
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

    // Recalcul des offsets (miroir de assembleGlobalSystem)
    std::vector<std::size_t> extColOffsets(nDomains + 1, 0);
    for(std::size_t i = 0; i < nDomains; i++)
        extColOffsets[i+1] = extColOffsets[i] + domains[i].nbExternalElements;

    std::vector<std::size_t> ifColOffsets(nInterfaces + 1, 0);
    for(std::size_t k = 0; k < nInterfaces; k++)
        ifColOffsets[k+1] = ifColOffsets[k] + interfaces[k].nbElements();

    auto domIdx = [&](const std::string &name) -> std::size_t {
        for(std::size_t i = 0; i < nDomains; i++)
            if(domains[i].name == name) return i;
        return 0;
    };

    // ---- Éléments externes ----
    for(std::size_t di = 0; di < nDomains; di++)
    {
        domainData   &dom  = domains[di];
        const std::size_t cOff = extColOffsets[di];
        const std::size_t nD   = dom.nbDirichletElements;
        const std::size_t nE   = dom.nbExternalElements;

        for(std::size_t j = 0; j < nE; j++)
        {
            double val = X_global(cOff + j);
            if(j < nD)
                // Dirichlet j : l'inconnu était qn → bcValue[1]
                dom.elementVector[j].bcValue[1] = val;
            else
                // Neumann j : l'inconnu était phi → bcValue[0]
                dom.elementVector[j].bcValue[0] = val;
        }
    }

    // ---- Éléments d'interface ----
    for(std::size_t k = 0; k < nInterfaces; k++)
    {
        const interfaceStruct &iface   = interfaces[k];
        const std::size_t      nIf     = iface.nbElements();
        const std::size_t      potBase  = ifPotentialOffset + ifColOffsets[k];
        const std::size_t      fluxBase = ifFluxOffset      + ifColOffsets[k];

        const std::size_t idxA = domIdx(iface.domainA);
        const std::size_t idxB = domIdx(iface.domainB);

        for(std::size_t m = 0; m < nIf; m++)
                {
                    const double u_if  = X_global(potBase  + m);
                    const double qn_if = X_global(fluxBase + m);

                    const std::size_t jA = iface.elemIndicesA[m];
                    const std::size_t jB = iface.elemIndicesB[m];

                    domains[idxA].elementVector[jA].bcValue[0] = u_if;
                    domains[idxA].elementVector[jA].bcValue[1] = qn_if;

                    domains[idxB].elementVector[jB].bcValue[0] = u_if;
                    domains[idxB].elementVector[jB].bcValue[1] = -qn_if;

                    // --- AJOUT DEBUG ---
                    if(m == 0) { // On n'affiche que le premier élément pour ne pas flood
                        std::cout << "[DEBUG DISTRIBUTE] Interface: " << iface.name << " Elem index local interface: " << m << "\n"
                                  << "  -> X_global (u, qn): " << u_if << ", " << qn_if << "\n"
                                  << "  -> Dom A [" << domains[idxA].name << "] El " << jA << " reçoit bcValue: "
                                  << domains[idxA].elementVector[jA].bcValue[0] << ", " << domains[idxA].elementVector[jA].bcValue[1] << "\n"
                                  << "  -> Dom B [" << domains[idxB].name << "] El " << jB << " reçoit bcValue: "
                                  << domains[idxB].elementVector[jB].bcValue[0] << ", " << domains[idxB].elementVector[jB].bcValue[1] << std::endl;
                    }
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
    const std::string integrationType = "Gauss2";

    double startTotal = omp_get_wtime();
    double start, end;

    // ---- Étape 1 : éléments externes ----
    start = omp_get_wtime();
    auto domainByName = [&](const std::string &n) -> domainData& {
        for(auto &d : domains) if(d.name == n) return d;
        std::cerr << "Domain '" << n << "' not found.\n"; exit(1);
    };

    for(auto &dom : domains)
    {
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

    // ---- Étape 2 : éléments d'interface ----
    start = omp_get_wtime();
    for(auto &iface : interfaces)
    {
        domainData &dA = domainByName(iface.domainA);
        domainData &dB = domainByName(iface.domainB);
        addInterfaceElements(dA, dB, iface, nodalDisplacements);
    }
    end = omp_get_wtime();
    displayTime(start, end, "addInterfaceElements (all interfaces)", showTime);

    if(iteration == 1)
    {
        for(auto &dom : domains)
            std::cout << "  Domain '" << dom.name << "': "
                      << dom.nbExternalElements << " external + "
                      << (dom.elementVector.size() - dom.nbExternalElements)
                      << " interface = "
                      << dom.elementVector.size() << " total elements.\n";
    }

    // ---- Étape 3 : assemblage global ----
    start = omp_get_wtime();
    Eigen::MatrixXd A_global, b_global;
    std::size_t ifPotentialOffset, ifFluxOffset;
    assembleGlobalSystem(A_global, b_global, domains, interfaces,
                         integrationType, solutionType,
                         ifPotentialOffset, ifFluxOffset);
    end = omp_get_wtime();
    displayTime(start, end, "assembleGlobalSystem", showTime);

    // ---- Étape 4 : résolution ----
    start = omp_get_wtime();
    Eigen::MatrixXd X_global = A_global.partialPivLu().solve(b_global);
    std::cout << std::endl;
    std::cout << "X_global" << std::endl;
    std::cout << X_global << std::endl;
    std::cout << std::endl;
    end = omp_get_wtime();
    displayTime(start, end, "Solve global system", showTime);

    // ---- Étape 5 : distribution de la solution ----
    start = omp_get_wtime();
    distributeGlobalSolution(domains, interfaces, X_global, ifPotentialOffset, ifFluxOffset);
    end = omp_get_wtime();
    displayTime(start, end, "distributeGlobalSolution", showTime);

    // ---- DEBUG: vérification après distribution ----
    std::cout << "\n--- DEBUG distributeGlobalSolution ---\n";
    for(auto &dom : domains)
    {
        std::cout << "Domaine '" << dom.name << "' (" << dom.elementVector.size() << " elements):\n";
        double phiMin =  1e300, phiMax = -1e300;
        double qnMin  =  1e300, qnMax  = -1e300;
        int nbIface = 0;
        for(auto &el : dom.elementVector)
        {
            if(std::abs(el.bcValue[0]) < 1e15) { phiMin = std::min(phiMin, el.bcValue[0]); phiMax = std::max(phiMax, el.bcValue[0]); }
            qnMin = std::min(qnMin, el.bcValue[1]); qnMax = std::max(qnMax, el.bcValue[1]);
            if(!el.interfaceName.empty()) nbIface++;
        }
        std::cout << "  phi range : [" << phiMin << " , " << phiMax << "]\n";
        std::cout << "  qn  range : [" << qnMin  << " , " << qnMax  << "]\n";
        std::cout << "  interface elements : " << nbIface << "\n";
    }
    // Continuité à l'interface
    {
        auto domIdxF = [&](const std::string &n) -> int {
            for(int ii = 0; ii < (int)domains.size(); ii++)
                if(domains[ii].name == n) return ii;
            return -1;
        };
        for(auto &iface : interfaces)
        {
            int idxA = domIdxF(iface.domainA), idxB = domIdxF(iface.domainB);
            if(idxA < 0 || idxB < 0) continue;
            double maxDeltaPhi = 0.0, maxFluxErr = 0.0;
            for(std::size_t m = 0; m < iface.nbElements(); m++)
            {
                double phiA = domains[idxA].elementVector[iface.elemIndicesA[m]].bcValue[0];
                double phiB = domains[idxB].elementVector[iface.elemIndicesB[m]].bcValue[0];
                double qnA  = domains[idxA].elementVector[iface.elemIndicesA[m]].bcValue[1];
                double qnB  = domains[idxB].elementVector[iface.elemIndicesB[m]].bcValue[1];
                maxDeltaPhi = std::max(maxDeltaPhi, std::abs(phiA - phiB));
                maxFluxErr  = std::max(maxFluxErr, std::abs(qnA + qnB));
            }
            std::cout << "Interface '" << iface.name << "':\n";
            std::cout << "  max |phi_A - phi_B|               = " << maxDeltaPhi << "  (expected ~0)\n";
            std::cout << "  max |eps_A*qn_A + eps_B*qn_B|     = " << maxFluxErr  << "  (expected ~0)\n";
        }
    }
    std::cout << "--- END DEBUG ---\n\n";

    // ---- Étape 6 : pression électrostatique ----
    for(auto &dom : domains)
    {
        const double eps = dom.epsilon;
        for(auto &el : dom.elementVector)
        {
            const double pressure = eps * el.bcValue[1] * el.bcValue[1] / 2.0;
            //#pragma omp critical
            electrostaticPressure[el.tag] = pressure;
        }
    }

    // ---- Étape 7 : post-processing ----
    if(postProcessing)
    {
        int viewIndex = nbViews;
        std::vector<std::string> names;
        gmsh::model::list(names);

        for(auto &dom : domains)
        {
            std::vector<std::size_t> domainNodeTags;
            std::map<std::size_t, std::size_t> domainNodeIndices;
            std::vector<double> domainNodeCoords;
            std::vector<int> domainEntityTags;
            std::vector<std::vector<int>> domainElementTypes;
            std::vector<std::vector<std::vector<std::size_t>>> domainElementTags;
            std::map<std::size_t, std::vector<std::size_t>> domainElementNodeTags;

            getModel(domainNodeTags, domainNodeIndices, domainNodeCoords, domainEntityTags,
                     domainElementTypes, domainElementTags, domainElementNodeTags,
                     dom.physicalGroupDim, dom.physicalGroupTag);

            std::map<std::size_t, std::vector<std::size_t>> boundaryNodeToElements;
            computeNodeToElements(boundaryNodeToElements, dom.elementVector);

            std::map<int, double> internalNodalPhi;
            computeInternalPhi(internalNodalPhi, domainNodeCoords, domainNodeTags,
                               domainNodeIndices, dom.elementVector, integrationType,
                               boundaryNodeToElements, solutionType);

            std::vector<std::vector<std::vector<std::vector<double>>>> elementNodalPhi;
            elementNodeData(elementNodalPhi, domainElementTags, domainElementTypes,
                            domainElementNodeTags, dom.elementVector, internalNodalPhi,
                            dom.physicalGroupTag, boundaryNodeToElements);

            std::vector<std::vector<std::vector<std::vector<double>>>> elementalPhi;
            std::vector<std::vector<std::vector<std::vector<double>>>> elementalElectricField;
            computeElementData(elementalPhi, elementalElectricField, domainElementTags,
                               domainElementTypes, domainEntityTags, dom.elementVector,
                               solutionType, integrationType, false, computeElectricField);

            if(computePhi)
            {
                displayData(phiViewTag, names[0], "results.msh", "ElementNodeData",
                            domainElementTags, elementNodalPhi, viewIndex);
            }
            if(computeElectricField)
                displayData(electricFieldViewTag, names[0], "results.msh", "ElementData",
                            domainElementTags, elementalElectricField,
                            viewIndex + (computePhi ? 1 : 0));
        }
    }

    double endTotal = omp_get_wtime();
    displayTime(startTotal, endTotal, "total multi-domain BEM time", showTime);
}
