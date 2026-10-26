#pragma once
// ============================================================
//  PARAMÈTRES DU CODE BEM : le seul fichier à modifier
// ============================================================
#include <algorithm>
#include <set>
#include <string>
#include <vector>

// Champs publiables dans Gmsh (un champ = une vue)
enum class Field
{
    Phi,        // potentiel aux nœuds
    ENodal,     // champ E aux nœuds
    EBary,      // champ E aux barycentres des éléments
    ENorm,      // norme |E| aux nœuds (vue scalaire)
    Maxwell,    // tenseur de Maxwell Txx, Tyy, Txy, Tzz aux nœuds
    GradE2,     // grad |E|² aux barycentres (vecteur ; le plus coûteux)
    GradE2Norm, // norme |grad |E|²| aux barycentres (vue scalaire)
    StressJump  // saut de contrainte normale sur les interfaces "_inc"
                // (active aussi son calcul, utile pour le couplage FEM)
};

struct BEMParams
{
    // ---------------- Calcul ----------------
    std::string integration       = "Gauss30";  // quadrature Gmsh
    bool        analyticalKernels = true;       // true : noyaux analytiques, false : numériques
    bool        showTime          = true;       // affiche les chronos

    // ---------------- Visualisation ----------------
//    std::set<Field> fields = {Field::Phi, Field::ENodal, Field::ENorm, Field::Maxwell, Field::GradE2, Field::GradE2Norm, Field::StressJump};

    std::set<Field> fields = {Field::Phi, Field::ENodal, Field::ENorm, Field::Maxwell, Field::GradE2, Field::GradE2Norm};

    // ---------------- Forces sur les inclusions ----------------
    bool   inclusionForces  = false;    // affiche [Maxwell] ... Fx, Fy
    double channelThickness = 30e-6;   // épaisseur hors plan [m]
    // Interfaces concernées par la force ET par la vue du saut de contrainte (StressJump).
    // Vide = toutes les "_inc".
    // Exemple pour ta puce : { "BEM_interface_2_7_inc" } (la goutte d'eau)
    std::vector<std::string> inclusionInterfaces = {"BEM_interface_2_7_inc"};

    // Cette interface est-elle sélectionnée ci-dessus ?
    bool selectsInterface(const std::string &name) const
    {
        return inclusionInterfaces.empty()
            || std::find(inclusionInterfaces.begin(), inclusionInterfaces.end(), name)
               != inclusionInterfaces.end();
    }
};
