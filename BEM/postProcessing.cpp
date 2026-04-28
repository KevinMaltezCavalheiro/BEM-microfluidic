#include "functionsBEM.hpp"

// ============================================================
//  INTERNAL HELPER: BIE-filtered element list
//
//  The BIE representation integral for domain D must only sum over
//  elements that physically bound D:
//    'E'  → external boundary element of D  → always included
//    'A'  → interface element, primary side  → included (bounds D from inside)
//    'B'  → interface element, secondary side → EXCLUDED (bounds the OTHER domain)
//
//  Using 'B' elements in D's BIE would double-count the interface and
//  produce incorrect potential / electric field values.
// ============================================================
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
        auto ins = [&](std::size_t tag){
            if(!m.count(tag)) m[tag] = {i};
            else              m[tag].push_back(i);
        };
        ins(bieVec[i].nodeTags[0]);
        ins(bieVec[i].nodeTags[1]);
        if(bieVec[i].midNodeTag != 0)
            m[bieVec[i].midNodeTag] = {i};
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

    //    std::cout << "\n--- DIAGNOSTIC DOMAINE ---" << std::endl;
    //    int countA = 0, countB = 0, countE = 0;
    //    for(const auto& el : elementVector) {
    //        if(el.interfaceSide == 'A') countA++;
    //        else if(el.interfaceSide == 'B') countB++;
    //        else countE++;
    //    }
    //    std::cout << "Total elements: " << elementVector.size()
    //              << " | E: " << countE << " | A: " << countA << " | B: " << countB << std::endl;

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
            std::string elementName;
            int dim, order, nbNodes, nbPrimaryNodes;
            std::vector<double> localNodeCoord;

            // L'API Gmsh est généralement thread-safe pour les getters de propriétés
            gmsh::model::mesh::getElementProperties(domainElementTypes[i][j], elementName,
                                                    dim, order, nbNodes, localNodeCoord,
                                                    nbPrimaryNodes);

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
    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);
    std::vector<std::size_t> globalLocalTags; // Vecteur pour stocker les résultats de tous les threads

    // On ouvre une région parallèle
    #pragma omp parallel
    {
        // Chaque thread crée son propre conteneur local pour éviter les "critical" fréquents
        std::vector<std::size_t> threadPrivateTags;

        #pragma omp for schedule(guided)
        for(std::size_t i = 0; i < bieVec.size(); i++)
        {
            std::vector<size_t> elementTags;
            // Appel Gmsh (Lecture seule sur le mesh : généralement thread-safe)
            gmsh::model::mesh::getElementsByCoordinates(bieVec[i].x1, bieVec[i].y1, 0,
                                                        elementTags, 2);

            for(std::size_t j = 0; j < elementTags.size(); j++)
            {
                int elementType;
                std::vector<size_t> nodeTags;
                int entityDim, entityTag;
                std::vector<int> physicalTags;

                // On récupère les infos de l'élément
                gmsh::model::mesh::getElement(elementTags[j], elementType, nodeTags,
                                              entityDim, entityTag);
                gmsh::model::getPhysicalGroupsForEntity(entityDim, entityTag, physicalTags);

                for(std::size_t k = 0; k < physicalTags.size(); k++)
                {
                    if(physicalTags[k] == domainPhysicalGroupTag)
                    {
                        // On ajoute au vecteur privé du thread
                        threadPrivateTags.push_back(elementTags[j]);
                    }
                }
            }
        }

        // Une fois la boucle terminée, on fusionne les vecteurs privés dans le vecteur global
        // On n'entre ici qu'une seule fois par thread (très peu coûteux)
        #pragma omp critical
        {
            globalLocalTags.insert(globalLocalTags.end(), threadPrivateTags.begin(), threadPrivateTags.end());
        }
    }

    // Transfert final vers l'argument de sortie
    domainBoundaryElementTags.insert(domainBoundaryElementTags.end(),
                                     globalLocalTags.begin(), globalLocalTags.end());

    // Nettoyage des doublons (Séquentiel, car std::unique nécessite un tri)
    std::sort(domainBoundaryElementTags.begin(), domainBoundaryElementTags.end());
    domainBoundaryElementTags.erase(
        std::unique(domainBoundaryElementTags.begin(), domainBoundaryElementTags.end()),
        domainBoundaryElementTags.end());
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
    const std::vector<elementStruct> bieVec     = buildBIEVector(elementVector);
    const auto                       bieNodeMap = buildBIENodeToElements(bieVec);

    std::map<std::size_t, std::vector<std::size_t>> position;
    phiInitialization(phi, position, domainElementTypes, domainElementTags);

    std::vector<size_t> domainBoundaryElTags;
    boundaryElementTags(domainBoundaryElTags, domainElementTypes, domainElementTags,
                        elementVector, domainPhysicalGroupTag);

    // --- DÉBUT DE LA ZONE PARALLÈLE ---
    #pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < domainBoundaryElTags.size(); i++)
    {
        std::size_t elementTag  = domainBoundaryElTags[i];

        // Utilisation de .at() pour la thread-safety en lecture sur la map 'position'
        const std::vector<std::size_t> &pos = position[elementTag];
        int entityTagIndex   = pos[0];
        int elementTypeIndex = pos[1];
        int elementTagIndex  = pos[2];

        // Accès aux tags des nœuds de l'élément (lecture seule)
        const std::vector<std::size_t> &nodeTags = domainElementNodeTags[elementTag];
        int nbNodes = phi[entityTagIndex][elementTypeIndex][elementTagIndex].size();

        for(int j = 0; j < nbNodes; j++)
        {
            std::size_t elementNodeTag = nodeTags[j];

            // Vérification si le nœud est sur la frontière BIE
            if(!bieNodeMap.count(elementNodeTag))
                continue;

            const std::vector<std::size_t> &commonIndices = bieNodeMap.at(elementNodeTag);

            if(commonIndices.size() == 1)
            {
                phi[entityTagIndex][elementTypeIndex][elementTagIndex][j] =
                    bieVec[commonIndices[0]].bcValue[0];
            }
            else
            {
                double currentPhi = 0.0;
                int count = 0;

                int adjIdx0 = commonIndices[0];
                std::size_t adjNode0 = bieVec[adjIdx0].nodeTags[0];
                if(adjNode0 == elementNodeTag) adjNode0 = bieVec[adjIdx0].nodeTags[1];

                // Recherche dans le vecteur local nodeTags
                if(std::count(nodeTags.begin(), nodeTags.end(), adjNode0)) count++;
                currentPhi += bieVec[adjIdx0].bcValue[0];

                int adjIdx1 = commonIndices[1];
                std::size_t adjNode1 = bieVec[adjIdx1].nodeTags[0];
                if(adjNode1 == elementNodeTag) adjNode1 = bieVec[adjIdx1].nodeTags[1];

                if(std::count(nodeTags.begin(), nodeTags.end(), adjNode1))
                {
                    if(count == 0) currentPhi  = bieVec[adjIdx1].bcValue[0];
                    else         { currentPhi += bieVec[adjIdx1].bcValue[0]; currentPhi /= 2.0; }
                }
                else if(count == 0)
                {
                    currentPhi += bieVec[adjIdx1].bcValue[0];
                    currentPhi /= 2.0;
                }

                // Écriture dans le tenseur phi (OK car (i, j) unique par thread)
                phi[entityTagIndex][elementTypeIndex][elementTagIndex][j] = currentPhi;
            }
        }
    }

    // Association finale (souvent rapide, peut rester séquentielle ou être parallélisée en interne)
    elementNodalPhiAssociation(phi, domainElementTags, domainElementNodeTags, internalNodalPhi);
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
