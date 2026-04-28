// ============================================================
//  Minimal_BEM_1Domain_Hole.geo
//  Carré (PDMS) avec trou circulaire interne
//  Orientation naturelle — le flip est géré par le code via getBoundary
// ============================================================

scale = 1;
epsilon_0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;

// ---- Paramètres ----
L = 1.0 * scale;
R = 0.2 * scale;
lc = 0.1 * scale;

// ============================================================
// Géométrie
// ============================================================

// ---- Points carré extérieur ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// ---- Points cercle (CCW naturel) ----
Point(5) = { 0,  0, 0, lc}; // Centre
Point(6) = { R,  0, 0, lc}; // Droite
Point(7) = { 0,  R, 0, lc}; // Haut
Point(8) = {-R,  0, 0, lc}; // Gauche
Point(9) = { 0, -R, 0, lc}; // Bas

// ============================================================
// Lignes
// ============================================================

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cercle CCW : 6 → 7 → 8 → 9 → 6
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// ============================================================
// Surfaces
// ============================================================

Curve Loop(1) = {1, 2, 3, 4};

// Trou : on inverse le signe pour que GMSH sache que c'est un trou
// getBoundary retournera des tags négatifs → ccw=false → flip dans le code
Curve Loop(2) = {-5, -6, -7, -8};

Plane Surface(1) = {1, 2};

// ============================================================
// Groupes physiques
// ============================================================

Physical Curve("Left", 101)      = {4};
Physical Curve("Right", 102)     = {2};
Physical Curve("TopBottom", 103) = {1, 3};
Physical Curve("Hole", 104)      = {5, 6, 7, 8};

Physical Surface("BEM_domain_1", 9) = {1};

// ============================================================
// Propriétés matériau & Conditions aux limites
// ============================================================

SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS);

SetNumber("Boundary Conditions/Left/BEM_domain_1/dirichlet",     1.0);
SetNumber("Boundary Conditions/Right/BEM_domain_1/dirichlet",    0.0);
SetNumber("Boundary Conditions/TopBottom/BEM_domain_1/neumann",  0);
SetNumber("Boundary Conditions/Hole/BEM_domain_1/dirichlet",     0.5);

Mesh.Algorithm = 6;
