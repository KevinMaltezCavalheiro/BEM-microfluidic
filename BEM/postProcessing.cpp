#include "functionsBEM.hpp"

static std::vector<elementStruct>
buildBIEVector(const std::vector<elementStruct> &elementVector)
{
    // Plus besoin de boucle ou d'OpenMP ici :
    // Le constructeur par copie est extrêmement rapide.
    return elementVector;
}
// Build a node→element map restricted to the BIE element list.
static std::map<std::size_t, std::vector<std::size_t>>
buildBIENodeToElements(const std::vector<elementStruct> &bieVec)
{
    std::map<std::size_t, std::vector<std::size_t>> m;

    for(std::size_t i = 0; i < bieVec.size(); i++)
    {
        m[bieVec[i].nodeTags[0]].push_back(i);
        m[bieVec[i].nodeTags[1]].push_back(i);

        if(bieVec[i].midNodeTag != 0)
            m[bieVec[i].midNodeTag].push_back(i);
    }

    return m;
}

// ============================================================
//  getModel  (unchanged from original)
// ============================================================
void getModel(std::vector<std::size_t> &domainNodeTags,
              std::map<std::size_t, std::size_t> &domainNodeIndices,
              std::vector<double> &domainNodeCoords,
              std::vector<int> &domainEntityTags,
              std::vector<std::vector<int>> &domainElementTypes,
              std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
              std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
              const int domainPhysicalGroupDim,
              const int domainPhysicalGroupTag)
{
    gmsh::model::mesh::getNodesForPhysicalGroup(domainPhysicalGroupDim, domainPhysicalGroupTag,
                                                domainNodeTags, domainNodeCoords);
    for(std::size_t index = 0; index < domainNodeTags.size(); index++)
        domainNodeIndices[domainNodeTags[index]] = index;

    gmsh::model::getEntitiesForPhysicalGroup(domainPhysicalGroupDim, domainPhysicalGroupTag,
                                             domainEntityTags);
    domainElementTypes.resize(domainEntityTags.size());
    domainElementTags.resize(domainEntityTags.size());

    for(std::size_t j = 0; j < domainEntityTags.size(); j++)
    {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags,
                                       domainPhysicalGroupDim, domainEntityTags[j]);

        domainElementTypes[j].resize(elementTypes.size());
        domainElementTags[j].resize(elementTypes.size());

        for(std::size_t k = 0; k < elementTypes.size(); k++)
        {
            int nbNodesPerElement = (int)(nodeTags[k].size() / elementTags[k].size());
            std::vector<std::size_t> elementNodeTags(nbNodesPerElement);

            domainElementTypes[j][k] = elementTypes[k];
            domainElementTags[j][k].resize(elementTags[k].size());

            for(std::size_t l = 0; l < elementTags[k].size(); l++)
            {
                domainElementTags[j][k][l] = elementTags[k][l];
                int n = 0;
                for(std::size_t mm = l*nbNodesPerElement; mm < (l+1)*nbNodesPerElement; mm++)
                    elementNodeTags[n++] = nodeTags[k][mm];
                domainElementNodeTags[elementTags[k][l]] = elementNodeTags;
            }
        }
    }
}

// ============================================================
//  computeNodeToElements  (unchanged from original)
// ============================================================
void computeNodeToElements(std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
                           const std::vector<elementStruct> &elementVector)
{
    for(std::size_t i = 0; i < elementVector.size(); i++)
    {
        auto ins = [&](std::size_t tag){
            if(!boundaryNodeToElements.count(tag)) boundaryNodeToElements[tag] = {i};
            else                                   boundaryNodeToElements[tag].push_back(i);
        };
        ins(elementVector[i].nodeTags[0]);
        ins(elementVector[i].nodeTags[1]);
        if(elementVector[i].midNodeTag != 0)
            boundaryNodeToElements[elementVector[i].midNodeTag] = {i};
    }
}

// ============================================================
//  computeInternalPhi
// ============================================================

void computeInternalPhi(std::map<int, double> &internalNodalPhi,
                        const std::vector<double> &domainNodeCoords,
                        const std::vector<std::size_t> &domainNodeTags,
                        std::map<std::size_t, std::size_t> &domainNodeIndices,
                        const std::vector<elementStruct> &elementVector,
                        const std::string integrationType,
                        std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
                        const bool &solutionType)
{
    // On construit le vecteur BIE une seule fois avant la boucle parallèle
    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);

        std::cout << "\n--- DIAGNOSTIC DOMAINE ---" << std::endl;
        int countA = 0, countB = 0, countE = 0;
        for(const auto& el : elementVector) {
            if(el.interfaceSide == 'A') countA++;
            else if(el.interfaceSide == 'B') countB++;
            else countE++;
        }
        std::cout << "Total elements: " << elementVector.size()
                  << " | E: " << countE << " | A: " << countA << " | B: " << countB << std::endl;

    std::cout << "Total nodes domain: " << domainNodeTags.size() << std::endl;

    // Vecteur temporaire pour éviter les conflits d'écriture dans la map
    // On l'initialise à 0.
    std::vector<double> internalNodalPhiTmp(domainNodeTags.size(), 0.0);

    std::vector<double> localCoord, weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    // Parallélisation sur les nœuds du domaine
    #pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        std::size_t nodeTag = domainNodeTags[i];

        // On ignore les nœuds qui sont sur la frontière
        if(boundaryNodeToElements.count(nodeTag)) continue;

        // Récupération des coordonnées (Lecture seule : OK)
        double x = domainNodeCoords[3 * domainNodeIndices[nodeTag]];
        double y = domainNodeCoords[3 * domainNodeIndices[nodeTag] + 1];

        double termGq = 0.0;
        double termHu = 0.0;

        // Boucle interne sur les éléments (Lecture seule sur bieVec : OK)
        for(const auto& el : bieVec)
        {
            double g = solutionType ? computeAnalyticalG(x, y, el)
                                    : computeNumericalG(x, y, el, integrationType, localCoord, weights);

            // Attention : computeAnalyticalH prend maintenant la normale en argument
            double h = computeAnalyticalH(x, y, el, el.normalOutward);

            termHu += h * el.bcValue[0];
            termGq += g * el.bcValue[1];
        }

        // Stockage dans le vecteur temporaire à l'indice i (propre à chaque thread)
        internalNodalPhiTmp[i] = termHu - termGq;
    }

    // Phase finale séquentielle : remplissage de la map
    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        // On ne remplit la map que si le nœud n'était pas sur la frontière
        if(!boundaryNodeToElements.count(domainNodeTags[i]))
        {
            internalNodalPhi[domainNodeTags[i]] = internalNodalPhiTmp[i];
        }
    }
}
// ============================================================
//  phiInitialization  (unchanged from original)
// ============================================================
void phiInitialization(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                       std::map<std::size_t, std::vector<std::size_t>> &position,
                       const std::vector<std::vector<int>> &domainElementTypes,
                       const std::vector<std::vector<std::vector<size_t>>> &domainElementTags)
{
    // Première dimension : tags d'entités (souvent petit, on ne parallélise pas ici)
    phi.resize(domainElementTypes.size());
    for(std::size_t i = 0; i < domainElementTypes.size(); i++)
    {
        // Deuxième dimension : types d'éléments
        phi[i].resize(domainElementTypes[i].size());

        // Parallélisation sur les types d'éléments au sein d'une entité
        #pragma omp parallel for schedule(guided) shared(phi, position)
        for(std::size_t j = 0; j < domainElementTypes[i].size(); j++)
        {
            int nbNodes;
            #pragma omp critical(gmsh_properties)
            {
                std::string elementName;
                int dim, order, nbPrimaryNodes;
                std::vector<double> localNodeCoord;
                gmsh::model::mesh::getElementProperties(domainElementTypes[i][j], elementName,
                                                        dim, order, nbNodes, localNodeCoord,
                                                        nbPrimaryNodes);
            }

            // Allocation de la troisième dimension (tags d'éléments)
            phi[i][j].resize(domainElementTags[i][j].size());

            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                // Allocation de la quatrième dimension (potentiel aux nœuds)
                phi[i][j][k].resize(nbNodes, HUGE_VAL);

                std::vector<std::size_t> pos = {i, j, k};

                // Section critique car std::map::operator[] modifie la structure interne
                #pragma omp critical(position_map_update)
                {
                    position[domainElementTags[i][j][k]] = pos;
                }
            }
        }
    }
}
// ============================================================
//  boundaryElementTags
//  BIE fix: only walk over 'E' and 'A' elements when searching for
//  2D elements adjacent to the boundary.
// ============================================================
void boundaryElementTags(std::vector<std::size_t> &domainBoundaryElementTags,
                         const std::vector<std::vector<int>> &domainElementTypes,
                         const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                         const std::vector<elementStruct> &elementVector,
                         const int domainPhysicalGroupTag)
{
    // 1. Récupération des tags d'entités (Surfaces) liées au groupe physique
    std::vector<int> entityTags; // Gmsh veut un vector<int> ici, pas de paires
    gmsh::model::getEntitiesForPhysicalGroup(2, domainPhysicalGroupTag, entityTags);

    std::set<std::size_t> elementInDomainSet;

    // 2. On boucle sur chaque entité (chaque surface) du groupe physique
    for (int tag : entityTags) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags;

        // On récupère les éléments de la surface d'index 'tag'
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2, tag);

        // On remplit notre Set avec tous les tags d'éléments trouvés
        for (const auto& group : elementTags) {
            for (std::size_t eTag : group) {
                elementInDomainSet.insert(eTag);
            }
        }
    }

    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);
    std::vector<std::size_t> globalLocalTags;

    // 2. Parallélisation avec protection de l'API
    #pragma omp parallel
    {
        std::vector<std::size_t> threadPrivateTags;

        #pragma omp for schedule(dynamic)
        for(std::size_t i = 0; i < bieVec.size(); i++)
        {
            std::vector<size_t> foundTags;

            // On protège l'accès à Gmsh.
            // Si le nommage (GMSH_API) échoue, utilise juste : #pragma omp critical
            #pragma omp critical
            {
                gmsh::model::mesh::getElementsByCoordinates(bieVec[i].x1, bieVec[i].y1, 0, foundTags, 2);
            }

            for(size_t tag : foundTags)
            {
                if(elementInDomainSet.find(tag) != elementInDomainSet.end())
                {
                    threadPrivateTags.push_back(tag);
                }
            }
        }

        // Fusion des résultats
        #pragma omp critical
        {
            globalLocalTags.insert(globalLocalTags.end(), threadPrivateTags.begin(), threadPrivateTags.end());
        }
    }

    // Remplissage du vecteur de sortie
    domainBoundaryElementTags = globalLocalTags;

    // Tri et suppression des doublons
    std::sort(domainBoundaryElementTags.begin(), domainBoundaryElementTags.end());
    domainBoundaryElementTags.erase(std::unique(domainBoundaryElementTags.begin(), domainBoundaryElementTags.end()), domainBoundaryElementTags.end());
}

// ============================================================
//  elementNodeData
//  BIE fix: use bieNodeMap so interface nodes reference only the
//  correct ('A'-side) element, not the mirror 'B'-side element.
// ============================================================
void elementNodeData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                     const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                     const std::vector<std::vector<int>> &domainElementTypes,
                     std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                     const std::vector<elementStruct> &elementVector,
                     std::map<int, double> &internalNodalPhi,
                     const int domainPhysicalGroupTag,
                     std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements)
{
    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);
    const auto bieNodeMap = buildBIENodeToElements(bieVec);

    std::map<std::size_t, std::vector<std::size_t>> position;
    phiInitialization(phi, position, domainElementTypes, domainElementTags);

    // ============================================================
    // 1. Reconstruction du potentiel aux nœuds BIE
    // ============================================================
    //
    // Le potentiel BEM est constant sur chaque élément 1D.
    // Au niveau d'un nœud partagé par plusieurs éléments BIE,
    // on reconstruit le potentiel nodal par moyenne.
    //
    std::map<std::size_t, double> boundaryNodePhi;

    for(const auto &[nodeTag, bieIndices] : bieNodeMap)
    {
        if(bieIndices.empty())
            continue;

        double sumPhi = 0.0;

        for(std::size_t bieIndex : bieIndices)
        {
            const auto &e = bieVec[bieIndex];
            sumPhi += e.bcValue[0];
        }

        boundaryNodePhi[nodeTag] =
            sumPhi / static_cast<double>(bieIndices.size());
    }


    // ============================================================
    // 2. Propagation aux éléments 2D
    // ============================================================
    //
    // IMPORTANT :
    // on parcourt TOUS les éléments 2D du domaine.
    //
    // Un élément peut toucher la frontière par un seul sommet
    // sans avoir lui-même d'arête BIE.
    //
    // Exemple :
    //
    //     5 ---- 63 ---- 64
    //            |
    //            |
    //          triangle
    //
    // Le triangle contenant 63 doit néanmoins recevoir phi(63).
    //
    std::vector<std::size_t> allElementTags;
    allElementTags.reserve(position.size());

    for(const auto &[elementTag, pos] : position)
        allElementTags.push_back(elementTag);

    #pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < allElementTags.size(); i++)
    {
        const std::size_t elementTag = allElementTags[i];

        const std::vector<std::size_t> &pos = position.at(elementTag);

        const int entityTagIndex   = pos[0];
        const int elementTypeIndex = pos[1];
        const int elementTagIndex  = pos[2];

        const std::vector<std::size_t> &nodeTags =
            domainElementNodeTags.at(elementTag);

        const int nbNodes =
            phi[entityTagIndex][elementTypeIndex][elementTagIndex].size();

        for(int j = 0; j < nbNodes; j++)
        {
            
            const std::size_t nodeTag = nodeTags[j];

            auto it = boundaryNodePhi.find(nodeTag);

            if(it != boundaryNodePhi.end())
            {
                phi[entityTagIndex][elementTypeIndex][elementTagIndex][j] =
                    it->second;
            }
        }
    }

    // ============================================================
    // 3. Association des nœuds internes
    // ============================================================
    //
    // Les nœuds qui n'ont pas reçu de valeur frontière sont encore
    // à HUGE_VAL et sont traités ici à partir de internalNodalPhi.
    //
    elementNodalPhiAssociation(phi,
                               domainElementTags,
                               domainElementNodeTags,
                               internalNodalPhi);
}
// ============================================================
//  elementNodalPhiAssociation  (unchanged from original)
// ============================================================
void elementNodalPhiAssociation(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                                const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                                std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                                std::map<int, double> &internalNodalPhi)
{
    // On boucle sur les entités et les types d'éléments
    for(std::size_t i = 0; i < domainElementTags.size(); i++)
    {
        for(std::size_t j = 0; j < domainElementTags[i].size(); j++)
        {
            // Parallélisation au niveau des éléments (k)
            // schedule(guided) est idéal car le nombre de noeuds (l) peut varier selon l'élément
            #pragma omp parallel for schedule(guided)
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                // Utilisation de .at() pour garantir un accès thread-safe en lecture seule
                const std::vector<std::size_t> &nodeTags = domainElementNodeTags[domainElementTags[i][j][k]];

                for(std::size_t l = 0; l < nodeTags.size(); l++)
                {
                    std::size_t nodeTag = nodeTags[l];

                    // On vérifie si la valeur est encore à l'état initial (HUGE_VAL)
                    // et si le potentiel nodal interne est disponible.
                    // .count() et .at() sur internalNodalPhi sont thread-safe en lecture seule.
                    if(phi[i][j][k][l] == HUGE_VAL && internalNodalPhi.count(nodeTag))
                    {
                        phi[i][j][k][l] = internalNodalPhi[nodeTag];
                    }
                }
            }
        }
    }
}

// ============================================================
//  computeElementData
//  BIE fix: only use 'E' and 'A' elements in the integral sum.
// ============================================================
void computeElementData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                        std::vector<std::vector<std::vector<std::vector<double>>>> &electricField,
                        const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                        const std::vector<std::vector<int>> &domainElementTypes,
                        const std::vector<int> &domainEntityTags,
                        const std::vector<elementStruct> &elementVector,
                        const bool solutionType,
                        const std::string integrationType,
                        const bool computePhi,
                        const bool computeElectricField)
{
    // Construction du vecteur BIE une seule fois
    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);

    std::vector<double> localCoord, weights;
    gmsh::model::mesh::getIntegrationPoints(1, integrationType, localCoord, weights);

    phi.resize(domainEntityTags.size());
    electricField.resize(domainEntityTags.size());

    for(std::size_t i = 0; i < domainEntityTags.size(); i++)
    {
        phi[i].resize(domainElementTypes[i].size());
        electricField[i].resize(domainElementTypes[i].size());

        for(std::size_t j = 0; j < domainElementTypes[i].size(); j++)
        {
            phi[i][j].resize(domainElementTags[i][j].size());
            electricField[i][j].resize(domainElementTags[i][j].size());

            // Récupération des barycentres pour le groupe actuel
            std::vector<double> barycenters;
            gmsh::model::mesh::getBarycenters(domainElementTypes[i][j], domainEntityTags[i],
                                              false, false, barycenters);

            // Parallélisation sur les éléments du domaine (k)
            #pragma omp parallel for schedule(guided)
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                double x = barycenters[3*k];
                double y = barycenters[3*k + 1];

                double currentPhi = 0.0, gradPhiX = 0.0, gradPhiY = 0.0;

                // Boucle interne sur les éléments de frontière (bieVec)
                for(std::size_t el = 0; el < bieVec.size(); el++)
                {
                    if(computePhi)
                    {
                        double h = computeAnalyticalH(x, y, bieVec[el], bieVec[el].normalOutward);
                        double g = solutionType ?
                            computeAnalyticalG(x, y, bieVec[el]) :
                            computeNumericalG(x, y, bieVec[el], integrationType, localCoord, weights);

                        currentPhi += h * bieVec[el].bcValue[0] - g * bieVec[el].bcValue[1];
                    }

                    if(computeElectricField)
                    {
                        double hX, gX, hY, gY;
                        if(!solutionType)
                        {
                            hX = computeNumericalGradHx(x, y, bieVec[el], integrationType, localCoord, weights);
                            gX = computeNumericalGradGx(x, y, bieVec[el], integrationType, localCoord, weights);
                            hY = computeNumericalGradHy(x, y, bieVec[el], integrationType, localCoord, weights);
                            gY = computeNumericalGradGy(x, y, bieVec[el], integrationType, localCoord, weights);
                        }
                        else
                        {
                            hX = computeAnalyticalGradHx(x, y, bieVec[el]);
                            gX = computeAnalyticalGradGx(x, y, bieVec[el]);
                            hY = computeAnalyticalGradHy(x, y, bieVec[el]);
                            gY = computeAnalyticalGradGy(x, y, bieVec[el]);
                        }
                        gradPhiX += hX * bieVec[el].bcValue[0] - gX * bieVec[el].bcValue[1];
                        gradPhiY += hY * bieVec[el].bcValue[0] - gY * bieVec[el].bcValue[1];
                    }
                }

                // Écriture des résultats (chaque k est unique, donc pas de conflit)
                phi[i][j][k].resize(1, currentPhi);
                electricField[i][j][k].resize(3);
                electricField[i][j][k][0] = -gradPhiX;
                electricField[i][j][k][1] = -gradPhiY;
                electricField[i][j][k][2] = 0.0;
            }
        }
    }
}

// ============================================================
//  displayData  (unchanged from original)
// ============================================================
void displayData(const int viewTag,
                 const std::string modelName,
                 const std::string fileName,
                 const std::string dataType,
                 const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
                 const std::vector<std::vector<std::vector<std::vector<double>>>> &data,
                 const int viewIndex)
{
    for(std::size_t i = 0; i < data.size(); i++)
        for(std::size_t j = 0; j < data[i].size(); j++)
            gmsh::view::addModelData(viewTag, 0, modelName, dataType,
                                     domainElementTags[i][j], data[i][j]);

    gmsh::option::setNumber("View[" + std::to_string(viewIndex) + "].IntervalsType", 3);
    gmsh::option::setNumber("View[" + std::to_string(viewIndex) + "].MaxRecursionLevel", 3);
    gmsh::option::setNumber("View[" + std::to_string(viewIndex) + "].TargetError", -0.0001);
    gmsh::option::setNumber("View[" + std::to_string(viewIndex) + "].Visible", 0);
}

void boundaryElementTags(std::vector<size_t> &domainBoundaryElTags,
                         const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                         const std::vector<std::vector<int>> &domainElementTypes,
                         const std::map<size_t, std::vector<size_t>> &domainElementNodeTags,
                         const std::map<size_t, std::vector<size_t>> &bieNodeMap)
{
    domainBoundaryElTags.clear();

    // On parcourt tous les éléments du domaine (Séquentiel pour éviter les conflits Gmsh)
    for(std::size_t i = 0; i < domainElementTags.size(); i++) {
        for(std::size_t j = 0; j < domainElementTags[i].size(); j++) {
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++) {

                size_t elTag = domainElementTags[i][j][k];
                const std::vector<size_t> &nodes = domainElementNodeTags.at(elTag);

                // Si UN SEUL des nœuds de l'élément est sur la frontière (bieNodeMap),
                // alors cet élément doit être traité par la boucle de post-proc.
                bool isOnBoundary = false;
                for(size_t nodeTag : nodes) {
                    if(bieNodeMap.count(nodeTag)) {
                        isOnBoundary = true;
                        break;
                    }
                }

                if(isOnBoundary) {
                    domainBoundaryElTags.push_back(elTag);
                }
            }
        }
    }
}
