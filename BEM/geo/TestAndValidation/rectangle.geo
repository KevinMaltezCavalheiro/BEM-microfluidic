// ============================================================
// Simple_Square_BEM.geo
// ============================================================

scale = 10;

eps0 = 8.8541878128e-12;
lc = 0.004 * scale;

// ============================================================
// POINTS POUR LE CARRÉ
// ============================================================

Point(1) = {0 * scale, 0 * scale, 0, lc};
Point(2) = {1 * scale, 0 * scale, 0, lc};
Point(3) = {1 * scale, 1 * scale, 0, lc};
Point(4) = {0 * scale, 1 * scale, 0, lc};

// ============================================================
// LIGNES DU DOMAINE
// ============================================================

Line(1) = {1, 2}; // Bas
Line(2) = {2, 3}; // Droite
Line(3) = {3, 4}; // Haut
Line(4) = {4, 1}; // Gauche

// ============================================================
// SURFACE
// ============================================================

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ============================================================
// GROUPES PHYSIQUES
// ============================================================

// Domaines
Physical Surface("BEM_domain_1", 1) = {1};

// Frontières extérieures
Physical Curve("BEM_left_E", 100) =
{
    4
};

Physical Curve("BEM_right_E", 101) =
{
    2
};

Physical Curve("BEM_top_bottom_E", 102) =
{
    1, 3
};

// ============================================================
// MATÉRIAUX
// ============================================================

SetNumber("Materials/BEM_domain_1/Epsilon", 1.0 * eps0);

// ============================================================
// SOLVEUR
// ============================================================

SetNumber("Non_linear_solver", 0);

// ============================================================
// CONDITIONS AUX LIMITES
// ============================================================

SetNumber("Boundary Conditions/BEM_left_E/BEM_domain_1/dirichlet", 200.0);
SetNumber("Boundary Conditions/BEM_right_E/BEM_domain_1/dirichlet", 100.0);
SetNumber("Boundary Conditions/BEM_top_bottom_E/BEM_domain_1/neumann", 0.0);


