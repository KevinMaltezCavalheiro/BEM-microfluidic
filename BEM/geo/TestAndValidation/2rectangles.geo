// ============================================================
//  Validation_BEM_Full_Logic_V2.geo
// ============================================================
// NOMENCLATURE DES GROUPES PHYSIQUES :
//
// _E    : Frontière extérieure
//         → condition imposée sur un seul domaine
//
// _jux  : Interface juxtaposée
//         → domaines A et B côte à côte
//
// _inc  : Interface incluse
//         → domaine A entoure le domaine B
//
// _H    : Trou / Hole
//         → frontière d'un trou dans un domaine
//
// FORMAT :
// BEM_interface_{ID_A}_{ID_B}_{TYPE}
// ============================================================


// ============================================================
// PARAMÈTRES GÉOMÉTRIQUES
// ============================================================

scale = 10;

// Toutes les dimensions géométriques sont définies
// dans une géométrie de référence puis multipliées par scale.

eps0 = 8.8541878128e-12;
lc = 0.04 * scale;

// ============================================================
// POINTS POUR LES BLOCS
// ============================================================

// Géométrie de référence :
// Domaine 1 : x = [0, 1]
// Domaine 2 : x = [1, 2.5]
// hauteur = 1

Point(1) = {0      * scale, 0 * scale, 0, lc};
Point(2) = {1      * scale, 0 * scale, 0, lc};
Point(3) = {2.5    * scale, 0 * scale, 0, lc};
Point(4) = {2.5    * scale, 1 * scale, 0, lc};
Point(5) = {1      * scale, 1 * scale, 0, lc};
Point(6) = {0      * scale, 1 * scale, 0, lc};


// ============================================================
// CERCLE — INCLUSION DANS LE DOMAINE 2
// ============================================================

cx = 1.75 * scale;
cy = 0.30 * scale;
r  = 0.15 * scale;

Point(10) = {cx,       cy,       0, lc};
Point(11) = {cx + r,   cy,       0, lc};
Point(12) = {cx,       cy + r,   0, lc};
Point(13) = {cx - r,   cy,       0, lc};
Point(14) = {cx,       cy - r,   0, lc};

Circle(20) = {11, 10, 12};
Circle(21) = {12, 10, 13};
Circle(22) = {13, 10, 14};
Circle(23) = {14, 10, 11};


// ============================================================
// CARRÉ — TROU DANS LE DOMAINE 2
// ============================================================

tx = 1.6  * scale;
ty = 0.7  * scale;
ts = 0.2  * scale;

Point(30) = {tx,      ty,      0, lc};
Point(31) = {tx + ts, ty,      0, lc};
Point(32) = {tx + ts, ty + ts, 0, lc};
Point(33) = {tx,      ty + ts, 0, lc};

Line(30) = {30, 31};
Line(31) = {31, 32};
Line(32) = {32, 33};
Line(33) = {33, 30};


// ============================================================
// LIGNES DES DOMAINES
// ============================================================

// Domaine 1
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};

// Domaine 2
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 5};


// ============================================================
// SURFACES
// ============================================================

// ------------------------------------------------------------
// Domaine 1 — rectangle gauche
// ------------------------------------------------------------

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};


// ------------------------------------------------------------
// Domaine 2 — rectangle droit avec inclusion + trou
// ------------------------------------------------------------

Curve Loop(2) = {5, 6, 7, -2};

Curve Loop(3) = {20, 21, 22, 23};
Curve Loop(4) = {30, 31, 32, 33};

Plane Surface(2) = {2, 3, 4};


// ------------------------------------------------------------
// Domaine 3 — goutte
// ------------------------------------------------------------

Curve Loop(5) = {20, 21, 22, 23};
Plane Surface(3) = {5};


// ============================================================
// GROUPES PHYSIQUES
// ============================================================

// ------------------------------------------------------------
// Domaines
// ------------------------------------------------------------

Physical Surface("BEM_domain_1", 1) = {1};
Physical Surface("BEM_domain_2", 2) = {2};
Physical Surface("BEM_domain_3", 3) = {3};


// ------------------------------------------------------------
// Frontières extérieures
// ------------------------------------------------------------

Physical Curve("BEM_left_E", 100) =
{
    4
};

Physical Curve("BEM_right_E", 101) =
{
    6
};

Physical Curve("BEM_top_bottom_E", 102) =
{
    1, 3, 5, 7
};


// ------------------------------------------------------------
// Interfaces
// ------------------------------------------------------------

// Domaine 1 ↔ Domaine 2
Physical Curve("BEM_interface_1_2_jux", 103) =
{
    2
};

// Domaine 2 ↔ Domaine 3
Physical Curve("BEM_interface_2_3_inc", 104) =
{
    20:23
};


// ------------------------------------------------------------
// Trou / électrode
// ------------------------------------------------------------

Physical Curve("BEM_hole_electrode_H", 105) =
{
    30:33
};


// ============================================================
// MATÉRIAUX
// ============================================================

SetNumber("Materials/BEM_domain_1/Epsilon",
          1.0 * eps0);

SetNumber("Materials/BEM_domain_2/Epsilon",
          2.0 * eps0);

SetNumber("Materials/BEM_domain_3/Epsilon",
          80.0 * eps0);


// ============================================================
// SOLVEUR
// ============================================================

SetNumber("Non_linear_solver", 0);


// ============================================================
// CONDITIONS AUX LIMITES
// ============================================================

// Électrode gauche
SetNumber(
    "Boundary Conditions/BEM_left_E/BEM_domain_1/dirichlet",
    200.0
);

// Électrode droite
SetNumber(
    "Boundary Conditions/BEM_right_E/BEM_domain_2/dirichlet",
    100.0
);


// Parois horizontales — Neumann nul
SetNumber(
    "Boundary Conditions/BEM_top_bottom_E/BEM_domain_1/neumann",
    0.0
);

SetNumber(
    "Boundary Conditions/BEM_top_bottom_E/BEM_domain_2/neumann",
    0.0
);


// Électrode dans le trou
SetNumber(
    "Boundary Conditions/BEM_hole_electrode_H/BEM_domain_2/dirichlet",
    150.0
);

//Mesh.RecombineAll = 1;
//Mesh.ElementOrder = 2;
