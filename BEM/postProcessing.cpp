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
    std::vector<elementStruct> out;
    out.reserve(elementVector.size());
    for(const auto &el : elementVector)
    {
        // On accepte TOUT ce qui appartient au domaine local.
        // 'E' (Extérieur), 'A' (Interface côté 1), 'B' (Interface côté 2)
        // La distinction A/B ne sert plus à "exclure",
        // mais juste à savoir d'où vient l'élément.
        out.push_back(el);
    }
    return out;
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
//  BIE fix: only use 'E' and 'A' elements in the integral sum.
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
    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);

        // --- PRINT DE CONTRÔLE AVANT LA BOUCLE ---
        // On veut savoir si le Domaine 2 a bien reçu ses éléments d'interface
        bool isDom2 = (bieVec.size() > 0 && bieVec[0].interfaceSide == 'B');
        if(isDom2) {
            std::cout << "\n[CHECK DOM2] Nombre d'éléments dans bieVec : " << bieVec.size() << std::endl;
            for(size_t k=0; k < bieVec.size(); k++) {
                std::cout << "  -> Elem " << k << " Tag: " << bieVec[k].tag
                          << " Phi: " << bieVec[k].bcValue[0]
                          << " Qn: " << bieVec[k].bcValue[1] << std::endl;
            }
        }

    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        if(boundaryNodeToElements.count(domainNodeTags[i])) continue;

        double x = domainNodeCoords[3 * domainNodeIndices[domainNodeTags[i]]];
        double y = domainNodeCoords[3 * domainNodeIndices[domainNodeTags[i]] + 1];

        double termGq = 0.0;
        double termHu = 0.0;
        double sH = 0.0;

        // --- LE PRINT QUE TU DEMANDES ---
        // On cible un point au centre de ton inclusion (x proche de 0, y proche de 0)
        bool targetNode = (std::abs(x) < 0.05 && std::abs(y) < 0.05);
        if(targetNode) {
            std::cout << "\n--- ANALYSE CRITIQUE NOEUD " << domainNodeTags[i] << " (x=" << x << ", y=" << y << ") ---" << std::endl;
        }

        for(const auto& el : bieVec)
        {
            double h = computeAnalyticalH(x, y, el);
            double g = solutionType ? computeAnalyticalG(x, y, el) : computeNumericalG(x, y, el, integrationType);

            if(targetNode) {
                std::cout << "  Elem " << el.tag
                          << " | h: " << h
                          << " | g: " << g
                          << " | phi: " << el.bcValue[0]
                          << " | qn: " << el.bcValue[1] << std::endl;
            }

            // Pour les éléments 'B', la normale stockée est inversée par rapport
            // à la normale sortante du domaine courant → compenser le signe de H.
            // G ne dépend pas de la normale → inchangé.
            double sign = (el.interfaceSide == 'B') ? -1.0 : 1.0;
            termHu += sign * h * el.bcValue[0];
            termGq += sign * g * el.bcValue[1];

            std::cout << "elTag: " << el.tag << std::endl;
            std::cout << "bcValue[1]: " << el.bcValue[1] << std::endl;
            std::cout << "g*bcValue[1]: " << g*el.bcValue[1] << std::endl;

        }

        if(targetNode) {
            std::cout << "  TOTAL -> sumH: " << sH << " | sumHu: " << termHu << " | sumGq: " << termGq << std::endl;
            std::cout << "  RESULTAT THEORIQUE (Hu-Gq): " << (termHu - termGq) << std::endl;
            std::cout << "---------------------------------------------------\n" << std::endl;
        }

        internalNodalPhi[domainNodeTags[i]] = termHu - termGq;
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
    phi.resize(domainElementTypes.size());
    for(std::size_t i = 0; i < domainElementTypes.size(); i++)
    {
        phi[i].resize(domainElementTypes[i].size());
        //#pragma omp parallel for schedule(guided) \
            shared(phi, domainElementTypes, domainElementTags, position)
        for(std::size_t j = 0; j < domainElementTypes[i].size(); j++)
        {
            std::string elementName;
            int dim, order, nbNodes, nbPrimaryNodes;
            std::vector<double> localNodeCoord;
            gmsh::model::mesh::getElementProperties(domainElementTypes[i][j], elementName,
                                                    dim, order, nbNodes, localNodeCoord,
                                                    nbPrimaryNodes);
            phi[i][j].resize(domainElementTags[i][j].size());
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                phi[i][j][k].resize(nbNodes, HUGE_VAL);
                std::vector<std::size_t> pos = {i, j, k};
                //#pragma omp critical
                position[domainElementTags[i][j][k]] = pos;
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

    std::vector<std::size_t> localElementTags;

    for(std::size_t i = 0; i < bieVec.size(); i++)
    {
        std::vector<size_t> elementTags;
        gmsh::model::mesh::getElementsByCoordinates(bieVec[i].x1, bieVec[i].y1, 0,
                                                    elementTags, 2);
        for(std::size_t j = 0; j < elementTags.size(); j++)
        {
            int elementType;
            std::vector<size_t> nodeTags;
            int entityDim, entityTag;
            std::vector<int> physicalTags;
            gmsh::model::mesh::getElement(elementTags[j], elementType, nodeTags,
                                          entityDim, entityTag);
            gmsh::model::getPhysicalGroupsForEntity(entityDim, entityTag, physicalTags);
            for(std::size_t k = 0; k < physicalTags.size(); k++)
                if(physicalTags[k] == domainPhysicalGroupTag)
                    localElementTags.push_back(elementTags[j]);
        }
    }

    domainBoundaryElementTags.insert(domainBoundaryElementTags.end(),
                                     localElementTags.begin(), localElementTags.end());
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

    //#pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < domainBoundaryElTags.size(); i++)
    {
        std::size_t elementTag  = domainBoundaryElTags[i];
        int entityTagIndex      = position[elementTag][0];
        int elementTypeIndex    = position[elementTag][1];
        int elementTagIndex     = position[elementTag][2];
        int nbNodes             = phi[entityTagIndex][elementTypeIndex][elementTagIndex].size();
        std::vector<std::size_t> nodeTags = domainElementNodeTags[elementTag];

        for(int j = 0; j < nbNodes; j++)
        {
            std::size_t elementNodeTag = nodeTags[j];
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

                phi[entityTagIndex][elementTypeIndex][elementTagIndex][j] = currentPhi;
            }
        }
    }

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
    for(std::size_t i = 0; i < domainElementTags.size(); i++)
        for(std::size_t j = 0; j < domainElementTags[i].size(); j++)
        {
            //#pragma omp parallel for schedule(guided)
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                std::vector<std::size_t> nodeTags =
                    domainElementNodeTags[domainElementTags[i][j][k]];
                for(std::size_t l = 0; l < nodeTags.size(); l++)
                {
                    std::size_t nodeTag = nodeTags[l];
                    if(phi[i][j][k][l] == HUGE_VAL && internalNodalPhi.count(nodeTag))
                        phi[i][j][k][l] = internalNodalPhi[nodeTag];
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
    const std::vector<elementStruct> bieVec = buildBIEVector(elementVector);

    std::vector<double> barycenters;
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

            gmsh::model::mesh::getBarycenters(domainElementTypes[i][j], domainEntityTags[i],
                                              false, false, barycenters);

            //#pragma omp parallel for schedule(static)
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                double x = barycenters[3*k];
                double y = barycenters[3*k + 1];

                double currentPhi = 0.0, gradPhiX = 0.0, gradPhiY = 0.0;
                double h, hX, hY, g, gX, gY;

                for(std::size_t el = 0; el < bieVec.size(); el++)
                {
                    if(computePhi)
                    {
                        h = computeAnalyticalH(x, y, bieVec[el]);
                        g = solutionType ?
                            computeAnalyticalG(x, y, bieVec[el]) :
                            computeNumericalG(x, y, bieVec[el], integrationType);
                        currentPhi += h * bieVec[el].bcValue[0] - g * bieVec[el].bcValue[1];
                    }
                    if(computeElectricField)
                    {
                        if(!solutionType)
                        {
                            hX = computeNumericalGradHx(x, y, bieVec[el], integrationType);
                            gX = computeNumericalGradGx(x, y, bieVec[el], integrationType);
                            hY = computeNumericalGradHy(x, y, bieVec[el], integrationType);
                            gY = computeNumericalGradGy(x, y, bieVec[el], integrationType);
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
