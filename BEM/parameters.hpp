#pragma once
// ============================================================
//  PARAMÈTRES DU CODE BEM : le seul fichier à modifier
// ============================================================
#include <set>
#include <string>
#include <vector>

// Champs publiables dans Gmsh (un champ = une vue)
enum class Field
{
    Phi,        // potentiel aux nœuds
    ENodal,     // champ E aux nœuds
    EBary,      // champ E aux barycentres des éléments
    Maxwell,    // tenseur de Maxwell Txx, Tyy, Txy, Tzz (demande ENodal : ajouté automatiquement)
    GradE2,     // grad |E|² aux barycentres (le plus coûteux)
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
    std::set<Field> fields = {Field::Phi, Field::ENodal, Field::Maxwell, Field::GradE2, Field::StressJump};

    // ---------------- Forces sur les inclusions ----------------
    bool   inclusionForces  = true;    // affiche [Maxwell] ... Fx, Fy
    double channelThickness = 30e-6;   // épaisseur hors plan [m]
    // Interfaces concernées par la force. Vide = toutes les "_inc".
    // Exemple pour ta puce : { "BEM_interface_2_7_inc" } (la goutte d'eau)
    std::vector<std::string> inclusionInterfaces = {"BEM_interface_2_7_inc"};
};
