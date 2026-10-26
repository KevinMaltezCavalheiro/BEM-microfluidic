#include "functionsBEM.hpp"

// ============================================================
//  solverBEM  – top-level dispatcher
// ============================================================
void solverBEM(std::map<int, double> &electrostaticPressure,
               std::map<int, std::pair<double,double>> &nodalDisplacements,
               const bool postProcessing,
               const int iteration,
               const bool untangleMesh,
               const BEMParams &params)
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

    extraViewsStruct extra;

    extra.fields = params.fields;

    if(postProcessing)
        extra.createViews();

    std::vector<domainData> domains(domainNames.size());

    for(std::size_t i = 0; i < domainNames.size(); ++i)
    {
        domains[i].name = domainNames[i];
        domains[i].physicalGroupDim = domainDimTags[i].first;
        domains[i].physicalGroupTag = domainDimTags[i].second;
        domains[i].epsilon = getEpsilon("BEM_domain_" + domainNames[i]);
    }

    std::vector<interfaceStruct> interfaces;
    detectInterfaces(interfaces, domainNames);

    multiDomainSolverBEM(
        electrostaticPressure, nodalDisplacements,
        postProcessing, iteration, domains, interfaces, extra, params);

    if(iteration == 1)
        std::cout << "--------------------[BEM SOLVER]--------------------\n\n";
}

// ============================================================
//  computeInclusionLoads
//  Force de Maxwell sur les inclusions + saut de contrainte normale (couplage FEM)
// ============================================================
static void computeInclusionLoads(std::map<int, double> &electrostaticPressure,
                                  const std::vector<domainData> &domains,
                                  const std::vector<interfaceStruct> &interfaces,
                                  const BEMParams &params)
{
    const bool wantForces = params.inclusionForces;
    const bool wantJump   = params.fields.count(Field::StressJump) > 0;

    if (!wantForces && !wantJump)
        return;

    for (const auto &iface : interfaces)
    {
        // On ne traite ici que les inclusions.
        if (iface.name.find("_inc") == std::string::npos)
            continue;

        const auto itA = std::find_if(domains.begin(), domains.end(),
            [&iface](const domainData &d) { return d.name == iface.domainA; });
        const auto itB = std::find_if(domains.begin(), domains.end(),
            [&iface](const domainData &d) { return d.name == iface.domainB; });

        if (itA == domains.end() || itB == domains.end())
            throw std::runtime_error("Domain not found for interface: " + iface.name);

        // Force électrique résultante (liste vide = toutes les inclusions)
        const bool selected = params.selectsInterface(iface.name);

        if (wantForces && selected)
        {
            const maxwellLoadStruct load =
                computeInclusionMaxwellForce(*itA, *itB, iface, params.channelThickness);

            std::cout << "[Maxwell] " << iface.name
                      << " : Fx = " << load.Fx << " N"
                      << ", Fy = " << load.Fy << " N\n";
        }
        // Saut de contrainte normale (toutes les inclusions, comme la vue)
        if (wantJump)
            computeElectrostaticPressureJump(electrostaticPressure, domains, iface);
    }
}

// ============================================================
//  multiDomainSolverBEM
// ============================================================
void multiDomainSolverBEM(std::map<int,double> &electrostaticPressure,
                          std::map<int,std::pair<double,double>> &nodalDisplacements,
                          const bool postProcessing,
                          const int iteration,
                          std::vector<domainData> &domains,
                          std::vector<interfaceStruct> &interfaces,
                          const extraViewsStruct &extra,
                          const BEMParams &params)
{
    const bool showTime                 = params.showTime;
    const bool solutionType             = params.analyticalKernels;
    const std::string &integrationType  = params.integration;

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

    distributeGlobalSolution(
        domains, interfaces, X_global,
        ifPotentialOffset, ifFluxOffset);

    // Pression sur tous les éléments, comme dans l'ancien solveur mono-domaine
    computeElectrostaticPressureMultiDomain(
        electrostaticPressure, domains);

    end = omp_get_wtime();
    displayTime(start, end, "distributeGlobalSolution", showTime);

    // ---- Étape 6 : pression électrostatique (Parallèle) ----
    for (auto& domain : domains)
    {
        computeBoundaryElectricField(domain.elementVector);
        computeCartesianElectricField(domain.elementVector);
    }

    // ---- Forces et saut de contrainte sur les inclusions ----
    computeInclusionLoads(electrostaticPressure, domains, interfaces, params);

    const Quadrature quad = loadQuadrature(integrationType);
    // ---- Étape 7 : post-processing (Parallèle par domaine) ----
    if(postProcessing)
    {
        start = omp_get_wtime();
        std::vector<std::string> modelNames;
        gmsh::model::list(modelNames);
        std::string modelName = modelNames[0];

        std::vector<DomainMesh> meshes(domains.size());
        for(std::size_t d = 0; d < domains.size(); d++)
            meshes[d] = loadDomainMesh(domains[d]);

        for(std::size_t d = 0; d < domains.size(); d++)
            postProcessDomain(domains[d], meshes[d], quad, extra, modelName, solutionType);

        if (extra.wants(Field::StressJump) && extra.stressJumpViewTag >= 0)
        {
            #pragma omp critical(gmsh_api)
            {
                std::vector<std::string> modelNames;
                gmsh::model::list(modelNames);

                if (modelNames.empty())
                {
                    std::cerr << "[ERROR] No Gmsh model available for stress-jump view.\n";
                }
                else
                {
                    const std::string modelName = modelNames[0];

                    for (const auto& iface : interfaces)
                    {
                        if (iface.name.find("_inc") == std::string::npos)
                            continue;

                        // Seulement les interfaces choisies dans parameters.hpp
                        if (!params.selectsInterface(iface.name))
                            continue;

                        const auto itA = std::find_if(
                            domains.begin(), domains.end(),
                            [&iface](const domainData& d)
                            {
                                return d.name == iface.domainA;
                            });

                        if (itA == domains.end())
                            continue;

                        displayInterfaceStressJump(
                            extra.stressJumpViewTag,
                            modelName,
                            iface,
                            *itA,
                            electrostaticPressure);
                    }
                }
            }
        }

        end = omp_get_wtime();
        displayTime(start, end, "Post-processing", showTime);
    }

    double endTotal = omp_get_wtime();
    displayTime(startTotal, endTotal, "total multi-domain BEM time", showTime);
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

    // 2. ASSEMBLAGE (boucle série sur les domaines ; parallélisme sur les lignes i de H et G)
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
        #pragma omp parallel for schedule(static)
        for(std::size_t i = 0; i < N; i++) {
            for(std::size_t j = 0; j < N; j++) {
                if(i == j) {
                    H(i,j) = -0.5;
                    G(i,j) = (dom.elementVector[i].length / (2.0*M_PI)) *
                             (std::log(dom.elementVector[i].length / 2.0) - 1.0);
                } else {
                    const elementStruct &ei = dom.elementVector[i];
                    const elementStruct &ej = dom.elementVector[j];

                    const Kernels k = computeAnalyticalKernels(
                        ei.midX,
                        ei.midY,
                        ej,
                        ej.normalOutward
                    );

                    H(i,j) = k.H;

                    if(solutionType) {
                        G(i,j) = k.G;
                    } else {
                        G(i,j) = computeNumericalG(
                            ei.midX,
                            ei.midY,
                            ej,
                            integrationType,
                            localCoord,
                            weights
                        );
                    }
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
    } // Fin de la boucle sur les domaines
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
            //domains[idxB].elementVector[jB].bcValue[1] = -qn_if; // Flux opposé sur l'autre face
            domains[idxB].elementVector[jB].bcValue[1] =
                -iface.kRatio * qn_if;
        }
    }
}

