#include "functionsBEM.hpp"

// static std::vector<elementStruct>
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

DomainMesh loadDomainMesh(const domainData &dom)
{
    DomainMesh m;
    getModel(m.nodeTags, m.nodeIndices, m.nodeCoords, m.entityTags,
             m.elemTypes, m.elemTags, m.elemNodeTags,
             dom.physicalGroupDim, dom.physicalGroupTag);

    m.barycenters.resize(m.entityTags.size());
    for(std::size_t i = 0; i < m.entityTags.size(); i++)
    {
        m.barycenters[i].resize(m.elemTypes[i].size());
        for(std::size_t j = 0; j < m.elemTypes[i].size(); j++)
            gmsh::model::mesh::getBarycenters(m.elemTypes[i][j], m.entityTags[i],
                                              false, false, m.barycenters[i][j]);
    }

    m.nbNodes.resize(m.entityTags.size());
    for(std::size_t i = 0; i < m.entityTags.size(); i++)
    {
        m.nbNodes[i].resize(m.elemTypes[i].size());
        for(std::size_t j = 0; j < m.elemTypes[i].size(); j++)
        {
            std::string elementName;
            int dim, order, nbPrimaryNodes;
            std::vector<double> localNodeCoord;
            gmsh::model::mesh::getElementProperties(m.elemTypes[i][j], elementName, dim, order,
                                                    m.nbNodes[i][j], localNodeCoord, nbPrimaryNodes);
        }
    }
    return m;
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
            boundaryNodeToElements[elementVector[i].midNodeTag].push_back(i);
    }
}

// ============================================================
//  computeInternalPhi
// ============================================================
void computeInternalPhi(
    std::map<int, double> &internalNodalPhi,
    const std::vector<double> &domainNodeCoords,
    const std::vector<std::size_t> &domainNodeTags,
    const std::map<std::size_t, std::size_t> &domainNodeIndices,
    const std::vector<elementStruct> &elementVector,
    const Quadrature &quad,
    std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
    const bool &solutionType)
{
    const std::vector<elementStruct> &bieVec = elementVector;

    std::vector<double> internalNodalPhiTmp(domainNodeTags.size(), 0.0);

    #pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        const std::size_t nodeTag = domainNodeTags[i];

        // Ignorer les nœuds situés sur la frontière
        if(boundaryNodeToElements.count(nodeTag))
            continue;

        const std::size_t nodeIndex = domainNodeIndices.at(nodeTag);
        const double x = domainNodeCoords[3 * nodeIndex];
        const double y = domainNodeCoords[3 * nodeIndex + 1];

        double termGq = 0.0;
        double termHu = 0.0;

        for(const auto &el : bieVec)
        {
            double g, h;

            if(solutionType)
            {
                const Kernels kernels = computeAnalyticalKernels(
                    x, y, el, el.normalOutward);

                g = kernels.G;
                h = kernels.H;
            }
            else
            {
                g = computeNumericalG(
                    x, y, el, quad.type,
                    quad.localCoord, quad.weights);

                h = computeAnalyticalH(
                    x, y, el, el.normalOutward);
            }

            termHu += h * el.bcValue[0];
            termGq += g * el.bcValue[1];
        }

        internalNodalPhiTmp[i] = termHu - termGq;
    }

    // Remplissage séquentiel de la map
    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        if(!boundaryNodeToElements.count(domainNodeTags[i]))
        {
            internalNodalPhi[domainNodeTags[i]] =
                internalNodalPhiTmp[i];
        }
    }
}

// ============================================================
//  phiInitialization  (unchanged from original)
// ============================================================
void phiInitialization(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                       std::map<std::size_t, std::vector<std::size_t>> &position,
                       const std::vector<std::vector<int>> &domainElementTypes,
                       const std::vector<std::vector<int>> &nbNodesPerType,
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
            const int nNodes = nbNodesPerType[i][j];

            // Allocation de la troisième dimension (tags d'éléments)
            phi[i][j].resize(domainElementTags[i][j].size());

            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                // Allocation de la quatrième dimension (potentiel aux nœuds)
                phi[i][j][k].resize(nNodes, HUGE_VAL);

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
//  elementNodeData
//  BIE fix: use bieNodeMap so interface nodes reference only the
//  correct ('A'-side) element, not the mirror 'B'-side element.
// ============================================================
void elementNodeData(std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
                     const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                     const std::vector<std::vector<int>> &domainElementTypes,
                     const std::vector<std::vector<int>> &nbNodesPerType,
                     const std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
                     const std::vector<elementStruct> &elementVector,
                     std::map<int, double> &internalNodalPhi,
                     const int domainPhysicalGroupTag,
                     std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements)
{
    const std::vector<elementStruct> &bieVec = elementVector;
    const auto bieNodeMap = buildBIENodeToElements(bieVec);

    std::map<std::size_t, std::vector<std::size_t>> position;
    phiInitialization(phi, position, domainElementTypes, nbNodesPerType, domainElementTags);

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
                                const std::map<std::size_t, std::vector<std::size_t>> &domainElementNodeTags,
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
                const std::vector<std::size_t> &nodeTags = domainElementNodeTags.at(domainElementTags[i][j][k]);

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
void computeElementData(
    std::vector<std::vector<std::vector<std::vector<double>>>> &phi,
    std::vector<std::vector<std::vector<std::vector<double>>>> &electricField,
    const DomainMesh &mesh,
    const std::vector<elementStruct> &elementVector,
    const bool solutionType,
    const Quadrature &quad,
    const bool computePhi,
    const bool computeElectricField)
{
    const std::vector<elementStruct> &bieVec = elementVector;

    phi.resize(mesh.entityTags.size());
    electricField.resize(mesh.entityTags.size());

    for(std::size_t i = 0; i < mesh.entityTags.size(); i++)
    {
        phi[i].resize(mesh.elemTypes[i].size());
        electricField[i].resize(mesh.elemTypes[i].size());

        for(std::size_t j = 0; j < mesh.elemTypes[i].size(); j++)
        {
            phi[i][j].resize(mesh.elemTags[i][j].size());
            electricField[i][j].resize(mesh.elemTags[i][j].size());

            const std::vector<double> &barycenters = mesh.barycenters[i][j];

            #pragma omp parallel for schedule(guided)
            for(std::size_t k = 0; k < mesh.elemTags[i][j].size(); k++)
            {
                const double x = barycenters[3*k];
                const double y = barycenters[3*k + 1];

                double currentPhi = 0.0;
                double gradPhiX = 0.0;
                double gradPhiY = 0.0;

                for(std::size_t el = 0; el < bieVec.size(); el++)
                {
                    const elementStruct &e = bieVec[el];

                    if(computePhi)
                    {
                        double h, g;

                        if(solutionType)
                        {
                            const Kernels kernels = computeAnalyticalKernels(
                                x, y, e, e.normalOutward);

                            h = kernels.H;
                            g = kernels.G;
                        }
                        else
                        {
                            h = computeAnalyticalH(
                                x, y, e, e.normalOutward);

                            g = computeNumericalG(
                                x, y, e, quad.type,
                                quad.localCoord, quad.weights);
                        }

                        currentPhi += h * e.bcValue[0]
                                    - g * e.bcValue[1];
                    }


                    if(computeElectricField)
                    {
                        double hX, gX, hY, gY;

                        if(solutionType)
                        {
                            const FirstKernels kernels =
                                computeAnalyticalFirstKernels(x, y, e);

                            hX = kernels.Hx;
                            gX = kernels.Gx;
                            hY = kernels.Hy;
                            gY = kernels.Gy;
                        }
                        else
                        {
                            hX = computeNumericalGradHx(
                                x, y, e, quad.type, quad.localCoord, quad.weights);

                            gX = computeNumericalGradGx(
                                x, y, e, quad.type, quad.localCoord, quad.weights);

                            hY = computeNumericalGradHy(
                                x, y, e, quad.type, quad.localCoord, quad.weights);

                            gY = computeNumericalGradGy(
                                x, y, e, quad.type, quad.localCoord, quad.weights);
                        }

                        gradPhiX += hX * e.bcValue[0] - gX * e.bcValue[1];
                        gradPhiY += hY * e.bcValue[0] - gY * e.bcValue[1];
                    }

                }

                if(computePhi)
                    phi[i][j][k].assign(1, currentPhi);

                if(computeElectricField)
                {
                    electricField[i][j][k].resize(3);
                    electricField[i][j][k][0] = -gradPhiX;
                    electricField[i][j][k][1] = -gradPhiY;
                    electricField[i][j][k][2] = 0.0;
                }
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
                 const std::vector<std::vector<std::vector<std::vector<double>>>> &data)
{
    for(std::size_t i = 0; i < data.size(); i++)
        for(std::size_t j = 0; j < data[i].size(); j++)
            gmsh::view::addModelData(viewTag, 0, modelName, dataType,
                                     domainElementTags[i][j], data[i][j]);

    setViewDefaults(viewTag);
}

void computeInternalElectricField(
    std::map<std::size_t, std::vector<double>> &internalNodalE,
    const std::vector<double> &domainNodeCoords,
    const std::vector<std::size_t> &domainNodeTags,
    const std::map<std::size_t, std::size_t> &domainNodeIndices,
    const std::vector<elementStruct> &elementVector,
    const Quadrature &quad,
    std::map<std::size_t, std::vector<std::size_t>> &boundaryNodeToElements,
    const bool &solutionType)
{
    const std::vector<elementStruct> &bieVec = elementVector;

    std::vector<std::vector<double>> internalNodalETmp(
        domainNodeTags.size(), std::vector<double>(3, 0.0));

    #pragma omp parallel for schedule(guided)
    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        const std::size_t nodeTag = domainNodeTags[i];

        // Les nœuds de frontière sont traités dans elementNodeElectricField.
        if(boundaryNodeToElements.count(nodeTag))
            continue;

        const std::size_t nodeIndex = domainNodeIndices.at(nodeTag);
        const double x = domainNodeCoords[3 * nodeIndex];
        const double y = domainNodeCoords[3 * nodeIndex + 1];

        double gradPhiX = 0.0;
        double gradPhiY = 0.0;

        for(const auto &el : bieVec)
        {
            double hX, gX, hY, gY;

            if(solutionType)
            {
                // Noyaux d'ordre 1 depuis une seule géométrie
                const FirstKernels fk = computeAnalyticalFirstKernels(x, y, el);
                hX = fk.Hx;  gX = fk.Gx;
                hY = fk.Hy;  gY = fk.Gy;
            }
            else
            {
                hX = computeNumericalGradHx(
                    x, y, el, quad.type, quad.localCoord, quad.weights);
                gX = computeNumericalGradGx(
                    x, y, el, quad.type, quad.localCoord, quad.weights);
                hY = computeNumericalGradHy(
                    x, y, el, quad.type, quad.localCoord, quad.weights);
                gY = computeNumericalGradGy(
                    x, y, el, quad.type, quad.localCoord, quad.weights);
            }

            gradPhiX += hX * el.bcValue[0] - gX * el.bcValue[1];
            gradPhiY += hY * el.bcValue[0] - gY * el.bcValue[1];
        }

        internalNodalETmp[i][0] = -gradPhiX;
        internalNodalETmp[i][1] = -gradPhiY;
        internalNodalETmp[i][2] = 0.0;
    }

    for(std::size_t i = 0; i < domainNodeTags.size(); i++)
    {
        const std::size_t nodeTag = domainNodeTags[i];

        if(!boundaryNodeToElements.count(nodeTag))
            internalNodalE[nodeTag] = internalNodalETmp[i];
    }
}

void elementNodeElectricField(
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
        &electricFieldNode,
    const std::vector<std::vector<std::vector<std::size_t>>>
        &domainElementTags,
    const std::map<std::size_t, std::vector<std::size_t>>
        &domainElementNodeTags,
    const std::vector<elementStruct> &elementVector,
    const std::map<std::size_t, std::vector<double>> &internalNodalE)
{
    // Accumulateur nodal : [somme Ex, somme Ey, nombre de contributions].
    std::map<std::size_t, std::vector<double>> boundaryNodeE;

    for(const auto &el : elementVector)
    {
        if(el.Ex == HUGE_VAL || el.Ey == HUGE_VAL)
            continue;

        auto addContribution =
            [&](const std::size_t nodeTag)
        {
            auto &v = boundaryNodeE[nodeTag];

            if(v.empty())
                v = std::vector<double>(3, 0.0);

            v[0] += el.Ex;
            v[1] += el.Ey;
            v[2] += 1.0;
        };

        addContribution(el.nodeTags[0]);
        addContribution(el.nodeTags[1]);

        if(el.midNodeTag != 0)
            addContribution(el.midNodeTag);
    }

    // Moyenne des contributions aux nœuds de frontière.
    for(auto &[nodeTag, v] : boundaryNodeE)
    {
        if(v[2] > 0.0)
        {
            v[0] /= v[2];
            v[1] /= v[2];
        }
    }

    // Allocation et remplissage des données ElementNodeData.
    electricFieldNode.resize(domainElementTags.size());

    for(std::size_t i = 0; i < domainElementTags.size(); i++)
    {
        electricFieldNode[i].resize(domainElementTags[i].size());

        for(std::size_t j = 0; j < domainElementTags[i].size(); j++)
        {
            electricFieldNode[i][j].resize(
                domainElementTags[i][j].size());

            for(std::size_t k = 0;
                k < domainElementTags[i][j].size(); k++)
            {
                const std::size_t elementTag =
                    domainElementTags[i][j][k];

                const auto &nodeTags =
                    domainElementNodeTags.at(elementTag);

                electricFieldNode[i][j][k].resize(
                    nodeTags.size(), std::vector<double>(3, 0.0));

                for(std::size_t l = 0; l < nodeTags.size(); l++)
                {
                    const std::size_t nodeTag = nodeTags[l];

                    const auto itB = boundaryNodeE.find(nodeTag);

                    if(itB != boundaryNodeE.end())
                    {
                        electricFieldNode[i][j][k][l][0] =
                            itB->second[0];
                        electricFieldNode[i][j][k][l][1] =
                            itB->second[1];
                    }
                    else
                    {
                        const auto itI = internalNodalE.find(nodeTag);

                        if(itI != internalNodalE.end())
                        {
                            electricFieldNode[i][j][k][l][0] =
                                itI->second[0];
                            electricFieldNode[i][j][k][l][1] =
                                itI->second[1];
                            electricFieldNode[i][j][k][l][2] =
                                itI->second[2];
                        }
                    }
                }
            }
        }
    }
}

void displayNodeElectricField(
    const int viewTag,
    const std::string modelName,
    const std::vector<std::vector<std::vector<std::size_t>>> &domainElementTags,
    const std::vector<
        std::vector<
            std::vector<
                std::vector<
                    std::vector<double>>>>>
        &data)
{
    for(std::size_t i = 0; i < data.size(); i++)
    {
        for(std::size_t j = 0; j < data[i].size(); j++)
        {
            std::vector<std::vector<double>> elementData;
            elementData.reserve(data[i][j].size());

            for(std::size_t k = 0; k < data[i][j].size(); k++)
            {
                std::vector<double> values;
                values.reserve(data[i][j][k].size() * 3);

                for(const auto &nodeField : data[i][j][k])
                {
                    values.push_back(nodeField[0]); // Ex
                    values.push_back(nodeField[1]); // Ey
                    values.push_back(nodeField[2]); // Ez
                }

                elementData.push_back(std::move(values));
            }

            gmsh::view::addModelData(
                viewTag,
                0,
                modelName,
                "ElementNodeData",
                domainElementTags[i][j],
                elementData,
                0.0,
                3
            );
        }
    }
    setViewDefaults(viewTag);
}

void computeBoundaryElectricField(std::vector<elementStruct> &elementVector)
{
    const std::vector<elementStruct> &bieVec = elementVector;

    const auto bieNodeMap = buildBIENodeToElements(bieVec);

    // 1. Reconstruction du potentiel aux nœuds de frontière.
    std::map<std::size_t, double> boundaryNodePhi;

    for(const auto &[nodeTag, bieIndices] : bieNodeMap)
    {
        if(bieIndices.empty())
            continue;

        double sumPhi = 0.0;

        for(const std::size_t index : bieIndices)
        {
            sumPhi += bieVec[index].bcValue[0];
        }

        boundaryNodePhi[nodeTag] =
            sumPhi / static_cast<double>(bieIndices.size());
    }

    // 2. Calcul de En et Et sur chaque élément BIE.
    for(auto &el : elementVector)
    {
        if(el.length <= 1e-15)
            continue;

        const auto it0 = boundaryNodePhi.find(el.nodeTags[0]);
        const auto it1 = boundaryNodePhi.find(el.nodeTags[1]);

        if(it0 == boundaryNodePhi.end() ||
           it1 == boundaryNodePhi.end())
            continue;

        // Convention provisoire : En = -d(phi)/dn.
        el.En = -el.bcValue[1];
        el.Et = -(it1->second - it0->second) / el.length;
    }
}

void computeCartesianElectricField(
    std::vector<elementStruct> &elementVector)
{
    for(auto &el : elementVector)
    {
        if(el.length <= 1e-15 ||
           el.En == HUGE_VAL ||
           el.Et == HUGE_VAL)
            continue;

        // Normale orientée vers l'extérieur du domaine.
        const double nxOut =
            el.normalOutward ? el.nx : -el.nx;
        const double nyOut =
            el.normalOutward ? el.ny : -el.ny;

        // Tangente orientée du nœud 0 vers le nœud 1.
        const double tx = el.deltaX / el.length;
        const double ty = el.deltaY / el.length;

        el.Ex = el.En * nxOut + el.Et * tx;
        el.Ey = el.En * nyOut + el.Et * ty;
    }
}

//éviter calcul en double
//eviter le calcul des deux représentation de E (elementdata et elementnodedata)
//paralleliser E elementnodedata et les autres fonctions

// ---------- Maxwell nodal, dérivé de E déjà calculé ----------
void computeMaxwellStressNodal(
    NodalScalarData& Txx,
    NodalScalarData& Tyy,
    NodalScalarData& Txy,
    NodalScalarData& Tzz,
    const NodalVectorData& nodalE,
    const double epsilon)
{
    Txx.resize(nodalE.size());
    Tyy.resize(nodalE.size());
    Txy.resize(nodalE.size());
    Tzz.resize(nodalE.size());

    for (std::size_t i = 0; i < nodalE.size(); ++i)
    {
        Txx[i].resize(nodalE[i].size());
        Tyy[i].resize(nodalE[i].size());
        Txy[i].resize(nodalE[i].size());
        Tzz[i].resize(nodalE[i].size());

        for (std::size_t j = 0; j < nodalE[i].size(); ++j)
        {
            Txx[i][j].resize(nodalE[i][j].size());
            Tyy[i][j].resize(nodalE[i][j].size());
            Txy[i][j].resize(nodalE[i][j].size());
            Tzz[i][j].resize(nodalE[i][j].size());

            for (std::size_t k = 0; k < nodalE[i][j].size(); ++k)
            {
                const std::size_t nn = nodalE[i][j][k].size();

                Txx[i][j][k].resize(nn);
                Tyy[i][j][k].resize(nn);
                Txy[i][j][k].resize(nn);
                Tzz[i][j][k].resize(nn);

                for (std::size_t l = 0; l < nn; ++l)
                {
                    const double Ex = nodalE[i][j][k][l][0];
                    const double Ey = nodalE[i][j][k][l][1];

                    const double E2 = Ex * Ex + Ey * Ey;
                    const double halfEps = 0.5 * epsilon;

                    Txx[i][j][k][l] = halfEps * (Ex * Ex - Ey * Ey);
                    Tyy[i][j][k][l] = halfEps * (Ey * Ey - Ex * Ex);
                    Txy[i][j][k][l] = epsilon * Ex * Ey;
                    Tzz[i][j][k][l] = -halfEps * E2;
                }
            }
        }
    }
}


void computeElementGradE2(
    std::vector<std::vector<std::vector<std::vector<double>>>> &gradE2,
    const DomainMesh &mesh,
    const std::vector<elementStruct> &elementVector,
    const Quadrature &quad,
    const bool &solutionType,
    const double prefactor)
{
    const std::vector<elementStruct> &bieVec = elementVector;

    gradE2.resize(mesh.entityTags.size());

    for(std::size_t i = 0; i < mesh.entityTags.size(); i++)
    {
        gradE2[i].resize(mesh.elemTypes[i].size());

        for(std::size_t j = 0; j < mesh.elemTypes[i].size(); j++)
        {
            gradE2[i][j].resize(mesh.elemTags[i][j].size());

            const std::vector<double> &barycenters =
                mesh.barycenters[i][j];

            #pragma omp parallel for schedule(guided)
            for(std::size_t k = 0; k < mesh.elemTags[i][j].size(); k++)
            {
                const double x = barycenters[3 * k];
                const double y = barycenters[3 * k + 1];

                double px = 0.0, py = 0.0;
                double pxx = 0.0, pxy = 0.0, pyy = 0.0;

                for(const auto &el : bieVec)
                {
                    const double phi = el.bcValue[0];
                    const double q   = el.bcValue[1];

                    if(solutionType)
                    {
                        // Noyaux d'ordre 1 et 2 depuis une seule géométrie
                        const GradKernels gk = computeAnalyticalGradKernels(x, y, el);

                        px += gk.Hx * phi - gk.Gx * q;
                        py += gk.Hy * phi - gk.Gy * q;

                        pxx += gk.sk.Hxx * phi - gk.sk.Gxx * q;
                        pxy += gk.sk.Hxy * phi - gk.sk.Gxy * q;
                        pyy += gk.sk.Hyy * phi - gk.sk.Gyy * q;
                    }
                    else
                    {
                        px += computeNumericalGradHx(
                                  x, y, el, quad.type,
                                  quad.localCoord, quad.weights) * phi
                            - computeNumericalGradGx(
                                  x, y, el, quad.type,
                                  quad.localCoord, quad.weights) * q;

                        py += computeNumericalGradHy(
                                  x, y, el, quad.type,
                                  quad.localCoord, quad.weights) * phi
                            - computeNumericalGradGy(
                                  x, y, el, quad.type,
                                  quad.localCoord, quad.weights) * q;

                        const SecondKernels sk = computeNumericalSecondKernels(
                            x, y, el, quad.type,
                            quad.localCoord, quad.weights);

                        pxx += sk.Hxx * phi - sk.Gxx * q;
                        pxy += sk.Hxy * phi - sk.Gxy * q;
                        pyy += sk.Hyy * phi - sk.Gyy * q;
                    }
                }
                gradE2[i][j][k] = {
                    2.0 * prefactor * (px * pxx + py * pxy),
                    2.0 * prefactor * (px * pxy + py * pyy),
                    0.0
                };
            }
        }
    }
}


// ---------- affichage (tag-based, sans index de vue) ----------
void setViewDefaults(const int tag)
{
    gmsh::view::option::setNumber(tag, "IntervalsType", 3);
    gmsh::view::option::setNumber(tag, "MaxRecursionLevel", 3);
    gmsh::view::option::setNumber(tag, "TargetError", -0.0001);
    gmsh::view::option::setNumber(tag, "Visible", 0);
}

static void addNodalScalarView(
    const int tag,
    const std::string& modelName,
    const std::vector<std::vector<std::vector<std::size_t>>>& elTags,
    const NodalScalarData& data)
{
    if (tag < 0 || data.empty())
        return;

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        for (std::size_t j = 0; j < data[i].size(); ++j)
        {
            if (data[i][j].empty())
                continue;

            gmsh::view::addModelData(
                tag, 0, modelName, "ElementNodeData",
                elTags.at(i).at(j), data[i][j], 0.0, 1);
        }
    }

    setViewDefaults(tag);
}

void computeAndDisplayExtraViews(
    const extraViewsStruct& extra,
    const std::string& modelName,
    const DomainMesh &mesh,
    const NodalVectorData& nodalE,
    const std::vector<elementStruct>& elementVector,
    const Quadrature &quad,
    const double epsilon,
    const bool &solutionType)
{
    if (extra.wants(Field::Maxwell) && !nodalE.empty())
    {
        NodalScalarData Txx, Tyy, Txy, Tzz;

        computeMaxwellStressNodal(
            Txx, Tyy, Txy, Tzz, nodalE, epsilon);

        #pragma omp critical(gmsh_api)
        {
            addNodalScalarView(
                extra.txxViewTag, modelName, mesh.elemTags, Txx);

            addNodalScalarView(
                extra.tyyViewTag, modelName, mesh.elemTags, Tyy);

            addNodalScalarView(
                extra.txyViewTag, modelName, mesh.elemTags, Txy);

            addNodalScalarView(
                extra.tzzViewTag, modelName, mesh.elemTags, Tzz);
        }
    }

    if (extra.wants(Field::GradE2) && extra.gradE2ViewTag >= 0)
    {
        std::vector<std::vector<std::vector<std::vector<double>>>> gradE2;

        computeElementGradE2(
            gradE2, mesh, elementVector, quad, solutionType, 1.0);

        #pragma omp critical(gmsh_api)
        {
            for (std::size_t i = 0; i < gradE2.size(); ++i)
            {
                for (std::size_t j = 0; j < gradE2[i].size(); ++j)
                {
                    gmsh::view::addModelData(
                        extra.gradE2ViewTag, 0, modelName, "ElementData",
                        mesh.elemTags.at(i).at(j),
                        gradE2[i][j], 0.0, 3);
                }
            }

            setViewDefaults(extra.gradE2ViewTag);
        }
    }
}

static std::pair<double, double> maxwellTraction(
    const elementStruct& el,
    const double epsilon,
    const double nx,
    const double ny)
{
    const double Ex = el.Ex;
    const double Ey = el.Ey;

    const double En = Ex * nx + Ey * ny;
    const double E2 = Ex * Ex + Ey * Ey;

    const double Tnx = epsilon * En * Ex
                     - 0.5 * epsilon * E2 * nx;

    const double Tny = epsilon * En * Ey
                     - 0.5 * epsilon * E2 * ny;

    return {Tnx, Tny};
}

maxwellLoadStruct computeInclusionMaxwellForce(
    const domainData& domainA,
    const domainData& domainB,
    const interfaceStruct& iface,
    const double thickness)
{
    maxwellLoadStruct result;

    if (iface.name.find("_inc") == std::string::npos)
        return result;

    if (iface.elemIndicesA.size() != iface.elemIndicesB.size())
        throw std::runtime_error(
            "Mismatched interface element counts: " + iface.name);

    for (std::size_t m = 0; m < iface.elemIndicesA.size(); ++m)
    {
        const auto& elA =
            domainA.elementVector.at(iface.elemIndicesA[m]);

        const auto& elB =
            domainB.elementVector.at(iface.elemIndicesB[m]);

        // Normale commune dirigée de B vers A.
        const double nx = elB.normalOutward ? elB.nx : -elB.nx;
        const double ny = elB.normalOutward ? elB.ny : -elB.ny;

        const auto tA = maxwellTraction(
            elA, domainA.epsilon, nx, ny);

        const auto tB = maxwellTraction(
            elB, domainB.epsilon, nx, ny);

        result.Fx += (tA.first - tB.first)
                   * elB.length * thickness;

        result.Fy += (tA.second - tB.second)
                   * elB.length * thickness;
    }

    return result;
}

double maxwellNormalStress(
    const elementStruct& el,
    const double epsilon,
    const double nx,
    const double ny)
{
    const double En = el.Ex * nx + el.Ey * ny;
    const double E2 = el.Ex * el.Ex + el.Ey * el.Ey;

    return 0.5 * epsilon * (2.0 * En * En - E2);
}

void displayInterfaceStressJump(
    const int viewTag,
    const std::string& modelName,
    const interfaceStruct& iface,
    const domainData& domainA,
    const std::map<int, double>& electrostaticPressure)
{
    std::vector<std::size_t> elementTags;
    std::vector<std::vector<double>> values;

    elementTags.reserve(iface.nbElements());
    values.reserve(iface.nbElements());

    for (std::size_t i = 0; i < iface.nbElements(); ++i)
    {
        const auto& elA =
            domainA.elementVector.at(iface.elemIndicesA.at(i));

        const auto it = electrostaticPressure.find(
            static_cast<int>(elA.tag));

        if (it == electrostaticPressure.end())
            continue;

        elementTags.push_back(elA.tag);
        values.push_back({it->second});
    }

    if (elementTags.empty())
    {
        std::cerr << "[WARNING] No interface stress data for "
                  << iface.name << '\n';
        return;
    }

    gmsh::view::addModelData(
        viewTag, 0, modelName, "ElementData",
        elementTags, values, 0.0, 1);

    setViewDefaults(viewTag);

    gmsh::view::option::setNumber(
        viewTag, "LineWidth", 15);
}

// ============================================================
//  Un champ = une fonction : calcul + publication dans sa vue
// ============================================================
static void postProcessPhi(const domainData &dom,
                           const DomainMesh &mesh,
                           const Quadrature &quad,
                           const extraViewsStruct &extra,
                           const std::string &modelName,
                           const bool solutionType,
                           std::map<std::size_t, std::vector<std::size_t>> &bNodeToElems)
{
    const double t0 = omp_get_wtime();

    std::map<int, double> internalPhi;
    std::vector<std::vector<std::vector<std::vector<double>>>> nPhi;

    computeInternalPhi(internalPhi,
                       mesh.nodeCoords, mesh.nodeTags, mesh.nodeIndices,
                       dom.elementVector, quad, bNodeToElems, solutionType);

    elementNodeData(nPhi,
                    mesh.elemTags, mesh.elemTypes, mesh.nbNodes, mesh.elemNodeTags,
                    dom.elementVector, internalPhi, dom.physicalGroupTag, bNodeToElems);

    #pragma omp critical(gmsh_api)
    {
        displayData(extra.phiViewTag, modelName, "results.msh", "ElementNodeData",
                    mesh.elemTags, nPhi);
    }

    displayTime(t0, omp_get_wtime(), "post-processing [" + dom.name + "] Phi", true);
}

static void postProcessENodal(const domainData &dom,
                              const DomainMesh &mesh,
                              const Quadrature &quad,
                              const extraViewsStruct &extra,
                              const std::string &modelName,
                              const bool solutionType,
                              std::map<std::size_t, std::vector<std::size_t>> &bNodeToElems,
                              NodalVectorData &nEField)   // sortie : relu par le tenseur de Maxwell
{
    const double t0 = omp_get_wtime();

    std::map<std::size_t, std::vector<double>> internalNodalE;

    computeInternalElectricField(internalNodalE,
                                 mesh.nodeCoords, mesh.nodeTags, mesh.nodeIndices,
                                 dom.elementVector, quad, bNodeToElems, solutionType);

    elementNodeElectricField(nEField,
                             mesh.elemTags, mesh.elemNodeTags,
                             dom.elementVector, internalNodalE);

    #pragma omp critical(gmsh_api)
    {
        displayNodeElectricField(extra.electricFieldNodeViewTag, modelName,
                                 mesh.elemTags, nEField);
    }

    displayTime(t0, omp_get_wtime(), "post-processing [" + dom.name + "] E nodal", true);
}

static void postProcessEBary(const domainData &dom,
                             const DomainMesh &mesh,
                             const Quadrature &quad,
                             const extraViewsStruct &extra,
                             const std::string &modelName,
                             const bool solutionType)
{
    const double t0 = omp_get_wtime();

    std::vector<std::vector<std::vector<std::vector<double>>>> ePhi;   // non publié
    std::vector<std::vector<std::vector<std::vector<double>>>> eField;

    computeElementData(ePhi, eField, mesh, dom.elementVector,
                       solutionType, quad, false, true);

    #pragma omp critical(gmsh_api)
    {
        displayData(extra.electricFieldViewTag, modelName, "results.msh", "ElementData",
                    mesh.elemTags, eField);
    }

    displayTime(t0, omp_get_wtime(), "post-processing [" + dom.name + "] E barycentres", true);
}

// ============================================================
//  postProcessDomain : liste des champs demandés, rien d'autre
// ============================================================
void postProcessDomain(const domainData &dom,
                       const DomainMesh &mesh,
                       const Quadrature &quad,
                       const extraViewsStruct &extra,
                       const std::string &modelName,
                       const bool solutionType)
{
    // Lien nœud -> éléments de frontière (nécessaire pour Phi et pour E nodal)
    std::map<std::size_t, std::vector<std::size_t>> bNodeToElems;
    if (extra.wants(Field::Phi) || extra.wants(Field::ENodal))
        computeNodeToElements(bNodeToElems, dom.elementVector);

    if (extra.wants(Field::Phi))
        postProcessPhi(dom, mesh, quad, extra, modelName, solutionType, bNodeToElems);

    // nEField est aussi lu par le tenseur de Maxwell
    NodalVectorData nEField;
    if (extra.wants(Field::ENodal))
        postProcessENodal(dom, mesh, quad, extra, modelName, solutionType,
                          bNodeToElems, nEField);

    if (extra.wants(Field::EBary))
        postProcessEBary(dom, mesh, quad, extra, modelName, solutionType);

    // Tenseur de Maxwell et grad|E|²
    if (extra.wants(Field::Maxwell) || extra.wants(Field::GradE2))
    {
        const double t0 = omp_get_wtime();
        computeAndDisplayExtraViews(extra, modelName, mesh, nEField,
                                    dom.elementVector, quad, dom.epsilon, solutionType);
        displayTime(t0, omp_get_wtime(),
                    "post-processing [" + dom.name + "] Maxwell + grad|E|^2", true);
    }
}
