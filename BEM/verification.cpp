#include "functionsBEM.hpp"

// Displays the time taken between the values of "start" and "end" (in seconds) if "showTime" is set to 'true'. "functionName" is
// the name of the function that has been timed.
void displayTime(const double start,
                 const double end,
                 std::string functionName, const bool showTime)
{
    if(showTime)
    {
        double duration = end - start;
        std::cout << functionName << "\n";
        std::cout << "Elapsed time: " << 1000000*duration << " [µs]\tor " << duration << " [s]\n\n";
    }
}

void convergenceFun(const std::vector<int> &domainEntityTags,
                    const std::vector<std::vector<int>> &domainElementTypes,
                    const std::vector<std::vector<std::vector<size_t>>> &domainElementTags,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> &elementalPhi,
                    const std::vector<std::vector<std::vector<std::vector<double>>>> &elementalElectricField)
{
    double l = 1.0;
    double phiLeft = 100.0;
    double phiRight = 200.0;
    double exactE = (phiLeft - phiRight) / l;

    double xLeft = 0.25;
    double xRight = 0.75;
    double yBottom = 0.25;
    double yTop = 0.75;

    double errorPhi = 0.0;
    double errorE = 0.0;
    double errorEWindow = 0.0;
    int count = 0;
    int countWindow = 0;

    std::vector<double> barycenters;

    for(std::size_t i = 0; i < domainEntityTags.size(); i++)
    {
        for(std::size_t j = 0; j < domainElementTypes[i].size(); j++)
        {
            gmsh::model::mesh::getBarycenters(domainElementTypes[i][j], domainEntityTags[i], false, false, barycenters);
            for(std::size_t k = 0; k < domainElementTags[i][j].size(); k++)
            {
                double numericalPhi = elementalPhi[i][j][k][0];
                double numericalE = elementalElectricField[i][j][k][0];

                double x = barycenters[3*k];
                double y = barycenters[3*k + 1];

                double exactPhi = phiLeft + (phiRight - phiLeft) * x / l;

                double relErrorPhi = abs( (exactPhi - numericalPhi) / exactPhi );
                double relErrorE = abs( (exactE - numericalE) / exactE);

                if(x > xLeft && x < xRight && y > yBottom && y < yTop)
                {
                    errorEWindow += relErrorE*relErrorE;
                    countWindow++;
                }

                errorPhi += relErrorPhi*relErrorPhi;
                errorE += relErrorE*relErrorE;
                count++;
            }
        }
    }
    errorPhi = sqrt(errorPhi / count);
    errorE = sqrt(errorE / count);
    errorEWindow = sqrt(errorEWindow / countWindow);

    std::cout << "\tRelative error on the electric potential: " << errorPhi << ".\n";
    std::cout << "\tRelative error on the electric field: " << errorE << ".\n";
    std::cout << "\tRelative error on the electric field (without the edges of the domain): " << errorEWindow << ".\n\n";
}

void printInterfaceFieldDiagnostics(
    const interfaceStruct& iface,
    const domainData& domainA,
    const domainData& domainB,
    const double epsA,
    const double epsB)
{
    const auto& indicesA = iface.elemIndicesA;
    const auto& indicesB = iface.elemIndicesB;
    const auto& elementsA = domainA.elementVector;
    const auto& elementsB = domainB.elementVector;

    

    if(indicesA.size() != indicesB.size())
    {
        std::cerr << "[FIELD CHECK] " << iface.name
                  << ": mismatched interface element counts: "
                  << indicesA.size() << " vs " << indicesB.size()
                  << '\n';
        return;
    }

    double maxRelEt = 0.0;
    double maxRelDn = 0.0;
    std::size_t validPairs = 0;

    std::cout << "\n[FIELD CHECK] " << iface.name << '\n';

    for(std::size_t i = 0; i < indicesA.size(); ++i)
    {
        if(indicesA[i] >= elementsA.size() ||
           indicesB[i] >= elementsB.size())
        {
            std::cerr << "  Invalid element index at pair " << i << '\n';
            continue;
        }

        const auto& a = elementsA[indicesA[i]];
        const auto& b = elementsB[indicesB[i]];
        const double xm = 0.5 * (a.x1 + a.x2);
        const double ym = 0.5 * (a.y1 + a.y2);

        if(a.Ex == HUGE_VAL || a.Ey == HUGE_VAL ||
           b.Ex == HUGE_VAL || b.Ey == HUGE_VAL ||
           a.length <= 1e-15)
            continue;

        // Repère commun défini par le côté A.
        const double nx = a.normalOutward ? a.nx : -a.nx;
        const double ny = a.normalOutward ? a.ny : -a.ny;
        const double tx = a.deltaX / a.length;
        const double ty = a.deltaY / a.length;

        const double EtA = a.Ex * tx + a.Ey * ty;
        const double EtB = b.Ex * tx + b.Ey * ty;

        const double DnA = epsA * (a.Ex * nx + a.Ey * ny);
        const double DnB = epsB * (b.Ex * nx + b.Ey * ny);

        const double relEt =
            std::abs(EtA - EtB) /
            std::max({std::abs(EtA), std::abs(EtB), 1e-12});

        const double relDn =
            std::abs(DnA - DnB) /
            std::max({std::abs(DnA), std::abs(DnB), 1e-30});

        maxRelEt = std::max(maxRelEt, relEt);
        maxRelDn = std::max(maxRelDn, relDn);
        ++validPairs;

        // Afficher seulement les cinq premières paires.
        if(i < 5)
        {
            std::cout
                << "  pair " << i
                << " | tags A/B: " << a.tag << "/" << b.tag
                << " | midpoint = (" << xm << ", " << ym << ")"
                << "\n    Et A/B = " << EtA << " / " << EtB
                << "\n    Dn A/B = " << DnA << " / " << DnB
                << "\n    relative errors: Et = " << relEt
                << ", Dn = " << relDn << '\n';
        }

        if(i == 0 && iface.name == "BEM_interface_1_2_jux")
        {
            std::cout << "[INTERFACE NODES]\n";

            for(int n = 0; n < 2; ++n)
            {
                const auto nodeA = a.nodeTags[n];
                const auto nodeB = b.nodeTags[n];

                std::cout
                    << "  n=" << n
                    << " | A tag=" << nodeA
                    << " phi=" << a.bcValue[0]
                    << " | B tag=" << nodeB
                    << " phi=" << b.bcValue[0]
                    << '\n';
            }
        }
    }

    std::cout
        << "  Valid pairs: " << validPairs
        << "\n  Maximum relative Et mismatch: " << maxRelEt
        << "\n  Maximum relative Dn mismatch: " << maxRelDn
        << "\n";
}

//Je suis en train de finir la lecture de ta partie multi-domaine et de postProcessing_1.cpp pour que le plan de refonte colle à ton code actuel.
//
//Ton intuition est la bonne : le problème n’est pas la physique mais le fait que les options décident comment on calcule, alors qu’elles devraient seulement décider ce qu’on affiche. Je n’ai rien compilé : tout ce qui suit est un plan et des squelettes, à adapter à tes signatures.
//
//1. Ce que je vois dans ton code actuel
//singleDomainSolverBEM et multiDomainSolverBEM font la même chose : assemblage, résolution, puis post-traitement dupliqué, avec des variantes de viewIndex à chaque endroit.
//Les décisions sont dispersées : computePhi, computeElectricField, electricFieldAtNodes, visualisation, postProcessing, extra.*, plus quatre tags de vue, nbViews et viewIndex. Chaque ajout (Maxwell, ∇|E|²) a ajouté une couche.
//Appels Gmsh dans la boucle parallèle, hors section critique : getModel, et getIntegrationPoints appelé à l’intérieur de computeInternalPhi et computeInternalElectricField. L’API Gmsh n’est pas thread-safe. Ça marche aujourd’hui, mais c’est une source de plantages aléatoires.
//Doublons de calcul : computeBoundaryElectricField tourne deux fois en multi-domaine (étape 6, puis dans la branche nodale), et buildBIEVector recopie le vecteur d’éléments à chaque appel.
//Mono-domaine : le potentiel est calculé (elementNodalPhi) mais jamais publié dans la branche visualisation.
//Parallélisme mal placé : la boucle parallèle sur les domaines enferme les boucles internes sur les nœuds (le parallélisme imbriqué est désactivé par défaut). Avec 2 domaines, seuls 2 threads travaillent sur la partie la plus coûteuse.
//Un DIAGNOSTIC DOMAINE s’affiche à chaque appel de computeInternalPhi.
//2. L’idée
//
//Un pipeline unique, toujours le même pour n ≥ 1 domaine (le mono-domaine est le cas n = 1 sans interface) :
//
//buildProblem → solveBEM → computeBoundaryFields → computeLoads → [post-traitement]
//                                                                   1. loadDomainMesh (série, seul à LIRE Gmsh)
//                                                                   2. computeFields  (zéro appel Gmsh)
//                                                                   3. publish        (série, seul à ÉCRIRE dans Gmsh)
//Étapes toujours exécutées et peu coûteuses : computeBoundaryFields (En, Et, Ex, Ey) et computeLoads (pression, forces). Elles ne dépendent d’aucune option.
//Seule option utilisateur : un std::set<Field> qui dit quels champs publier. Les dépendances (Txx a besoin de E) se résolvent automatiquement.
//Aucun bool ne traverse les fonctions de calcul.
//3. Structures
//cpp
//// ---- config.hpp ----
//struct SolverConfig {
//    std::string integration = "Gauss30";
//    bool analytic = true;          // false : noyaux numériques (tests)
//    bool timing   = true;
//};
//
//enum class Field   { Phi, E, Txx, Tyy, Txy, Tzz, GradE2 };
//enum class Support { Node, Element };      // ElementNodeData ou ElementData
//
//struct FieldSpec { Field id; const char *name; Support support; int nComp; };
//const std::vector<FieldSpec>& fieldSpecs();   // la table : 1 ligne par champ
//
//struct PostConfig { std::set<Field> publish = {Field::Phi, Field::E}; };
//cpp
//// ---- mesh/données : aucune de ces structures n'appelle Gmsh après sa construction ----
//using ElemData = std::vector<std::vector<std::vector<std::vector<double>>>>;  // [entité][type][élément][valeurs]
//
//struct Quadrature { std::vector<double> localCoord, weights; };  // 1 seul getIntegrationPoints
//
//struct DomainMesh {                       // enveloppe l'actuel getModel(...) + getBarycenters
//    std::vector<std::size_t> nodeTags;
//    std::map<std::size_t, std::size_t> nodeIndex;
//    std::vector<double> nodeCoords;
//    std::vector<int> entityTags;
//    std::vector<std::vector<int>> elemTypes;
//    std::vector<std::vector<std::vector<std::size_t>>> elemTags;
//    std::map<std::size_t, std::vector<std::size_t>> elemNodeTags;
//    std::vector<std::vector<std::vector<double>>> barycenters;   // [entité][type] -> xyz plat
//};
//DomainMesh loadDomainMesh(const domainData &d);
//
//struct DomainCtx {                        // tout ce qu'un calcul de champ a le droit de lire
//    const DomainMesh &mesh;
//    const domainData &dom;
//    const std::vector<elementStruct> &bie;                              // construit UNE fois
//    const std::map<std::size_t,std::vector<std::size_t>> &nodeToElems;  // idem
//    const Quadrature &quad;
//    const SolverConfig &cfg;
//};
//cpp
//// ---- fields.cpp : un champ = une entrée de table ----
//using FieldStore = std::map<Field, ElemData>;
//
//struct FieldComputer {
//    Field id;
//    std::vector<Field> deps;
//    std::function<ElemData(const DomainCtx&, const FieldStore&)> run;
//};
//
//FieldStore computeFields(const DomainCtx &ctx, const std::set<Field> &requested)
//{
//    FieldStore store;
//    std::function<void(Field)> need = [&](Field f) {
//        if(store.count(f)) return;
//        const FieldComputer &c = findComputer(f);
//        for(Field d : c.deps) need(d);
//        store[f] = c.run(ctx, store);
//    };
//    for(Field f : requested) need(f);
//    return store;
//}
//
//Exemple de table, en réutilisant tes fonctions existantes sans réécrire les maths :
//
//cpp
//{Field::E,   {},         [](auto &c, auto &) { /* computeInternalElectricField + elementNodeElectricField, puis aplatir */ }},
//{Field::Txx, {Field::E}, [](auto &c, auto &s) { return maxwellTxx(s.at(Field::E), c.dom.epsilon); }},
//{Field::GradE2, {},      [](auto &c, auto &) { /* ton computeElementGradE2, avec c.mesh.barycenters */ }},
//
//Aplatis les données nodales. NodalVectorData à 5 niveaux disparaît. Les valeurs par élément sont un vector<double> de taille nNoeuds × nComp, dans l’ordre nœud par nœud, ce qu’addModelData attend. displayNodeElectricField n’a alors plus de raison d’exister.
//
//cpp
//// ---- views.cpp : le SEUL code qui écrit dans les vues ----
//class ViewRegistry {
//public:
//    void publish(Field f, const std::string &model, int step, double time,
//                 const DomainMesh &mesh, const ElemData &data)
//    {
//        const FieldSpec &s = spec(f);
//        const int tag = tagOf(f);                         // crée la vue au 1er appel
//        const char *type = (s.support == Support::Node) ? "ElementNodeData" : "ElementData";
//        for(std::size_t i = 0; i < data.size(); i++)
//            for(std::size_t j = 0; j < data[i].size(); j++)
//                gmsh::view::addModelData(tag, step, model, type,
//                                         mesh.elemTags[i][j], data[i][j], time, s.nComp);
//    }
//private:
//    int tagOf(Field f);                                   // gmsh::view::add + options par tag
//    std::map<Field,int> tags_;
//};
//nbViews, viewIndex et viewIndexBase disparaissent : les options sont posées par tag via gmsh::view::option::setNumber.
//Dans la boucle de couplage, step = iteration ajoute des pas de temps à la même vue au lieu de créer de nouvelles vues.
//Si le maillage change entre itérations (untangleMesh, remaillage), les tags d’éléments changent : il faudra alors vider le registre (views.reset()).
//4. Le chef d’orchestre
//
//Il remplace solverBEM, singleDomainSolverBEM et multiDomainSolverBEM :
//
//cpp
//void runBEM(Problem &P, const SolverConfig &cfg, const PostConfig &post,
//            ViewRegistry &views, const int iteration, Loads &loads)
//{
//    solveBEM(P, cfg);                    // assemblage + résolution + distribution (n >= 1)
//    computeBoundaryFields(P);            // En, Et, Ex, Ey : UNE fois
//    loads = computeLoads(P, thickness);  // pression + forces (ta computeInclusionMaxwellForce)
//    if(post.publish.empty()) return;
//
//    const Quadrature quad = loadQuadrature(cfg);                 // phase 1 : série
//    std::vector<DomainMesh> meshes;
//    for(auto &d : P.domains) meshes.push_back(loadDomainMesh(d));
//
//    std::vector<FieldStore> stores(P.domains.size());            // phase 2 : zéro Gmsh
//    for(std::size_t d = 0; d < P.domains.size(); d++)            // série ici, boucles internes //
//        stores[d] = computeFields(makeCtx(P, meshes, d, quad, cfg), post.publish);
//
//    const std::string model = currentModelName();                // phase 3 : série
//    for(std::size_t d = 0; d < P.domains.size(); d++)
//        for(Field f : post.publish)
//            views.publish(f, model, iteration, iteration, meshes[d], stores[d].at(f));
//}
//Phase 2 en série sur les domaines. Maintenant que Gmsh n’y est plus, tu peux choisir. Je garderais les domaines en série et le omp parallel for sur les nœuds à l’intérieur : c’est là qu’est le coût, et il y a bien plus de nœuds que de domaines.
//main.cpp devient :
//cpp
//SolverConfig cfg;  PostConfig post;  ViewRegistry views;  Loads loads;
//post.publish = {Field::Phi, Field::E, Field::Txx, Field::Tyy, Field::Txy, Field::GradE2};
//Problem P = buildProblem(nodalDisplacements);
//runBEM(P, cfg, post, views, 1, loads);
//
//Cette sélection peut aussi être lue depuis des paramètres ONELAB (cases à cocher dans le GUI), sans recompiler.
//
//5. Correspondance ancien → nouveau
//Ancien    Nouveau
//solverBEM (détection domaines/interfaces, création des vues)    buildProblem + ViewRegistry (vues créées à la demande)
//singleDomainSolverBEM    supprimé : cas n = 1 de runBEM
//multiDomainSolverBEM étapes 1-5    solveBEM
//étape 6 + doublons de computeBoundaryElectricField    computeBoundaryFields (une fois)
//computeElectrostaticPressure / …MultiDomain / …Jump, ta force MST    computeLoads
//bloc post-traitement (étape 7)    phases 2 et 3
//displayData, displayNodeElectricField, computeAndDisplayExtraViews, addNodalScalarView, setViewDefaults    ViewRegistry::publish
//computeElementData (mode barycentres)    un Field de support Element, ou supprimé
//nbViews, viewIndex, extraViewsStruct, computePhi, computeElectricField, electricFieldAtNodes, visualisation, postProcessing    PostConfig::publish
//showTime, solutionType, integrationType codés en dur par fonction    SolverConfig
//convergenceFun (branche legacy)    fonction de test séparée, appliquée à un FieldStore
//
//Je suppose que tu abandonnes le mode barycentres pour E (electricFieldAtNodes = false). Si tu veux le garder, c’est un Field de plus avec Support::Element, sans aucun if.
//
//6. Ordre de migration
//
//À chaque étape : un commit, puis un test de non-régression.
//
//Référence. Pour 3 cas (carré Dirichlet seul, 2 carrés accolés en juxtaposition, goutte _inc), sauvegarde dans un fichier les bcValue triés par domaine/tag, les min/max/norme L² de φ et E nodaux, et les forces. Un petit script compare avec 1e-10 de tolérance relative. Rien de mathématique ne change, donc l’écart doit être ~1e-12.
//Nettoyage sans risque : supprimer le debug codé en dur (DIAGNOSTIC DOMAINE, DEBUG DISTRIBUTE), le code commenté, buildBIEVector en const &, un seul appel de computeBoundaryFields.
//Quadrature et DomainMesh. Les fonctions computeInternal* et computeElementGradE2 reçoivent ces données au lieu d’appeler Gmsh. Les omp critical(gmsh_api) disparaissent des calculs.
//ViewRegistry et données aplaties. Remplace les quatre fonctions d’affichage. Vérifie les vues visuellement, puis par gmsh::view::getModelData.
//Table FieldComputer et PostConfig.
//solveBEM pour n = 1. Avant de supprimer singleDomainSolverBEM, fais passer le cas mono-domaine par le chemin multi et compare les bcValue. Points d’attention :
//fillSystem et assembleGlobalSystem doivent avoir les mêmes conventions de signe.
//computeElectrostaticPressure écrivait x dans bcValue : c’est maintenant distributeGlobalSolution qui s’en charge.
//Avec 0 interface, nbExternalElements = elementVector.size(), donc la matrice reste carrée.
//computeLoads et main.
//Découpage en fichiers : kernels.cpp, assembly.cpp, fields.cpp, views.cpp, loads.cpp.
//
//Pour ajouter un champ plus tard (par exemple une densité de force) : une ligne dans fieldSpecs(), une entrée dans la table de calcul, et l’ajouter au std::set dans main. Rien d’autre à toucher.
//
//Pour que je te produise le code de l’étape 2 (Quadrature/DomainMesh et le changement de signatures de computeInternalPhi, computeInternalElectricField et computeElementGradE2), il me suffit de ton feu vert. C’est l’étape qui élimine le risque thread-safety et qui conditionne les suivantes. Je la ferais sur tes fichiers _1.
