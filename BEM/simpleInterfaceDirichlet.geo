// ============================================================
//  Minimal_BEM_1Domain_Hole_V2_FIXED.geo
// ============================================================

scale = 1;
epsilon_0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;

L = 1.0 * scale;
R = 0.2 * scale;
lc = 0.05 * scale;

// ---- 1. CADRE EXTÉRIEUR ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};

// ---- 2. LE TROU (Électrode) ----
Point(5) = { 0,  0, 0, lc};
Point(6) = { R,  0, 0, lc};
Point(7) = { 0,  R, 0, lc};
Point(8) = {-R,  0, 0, lc};
Point(9) = { 0, -R, 0, lc};

Circle(5) = {6, 5, 7}; Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9}; Circle(8) = {9, 5, 6};

// ---- SURFACES ----
Curve Loop(1) = {1, 2, 3, 4}; // Boucle extérieure
Curve Loop(2) = {5, 6, 7, 8}; // Boucle du trou

// Domaine 1 : Le PDMS est la surface 1 perçée par la boucle 2
Plane Surface(1) = {1, 2};

// ---- GROUPES PHYSIQUES (Convention Nomenclature) ----
Physical Curve("BEM_outside_E", 100) = {1, 2, 3, 4};
Physical Curve("BEM_Hole_H", 106)      = {5, 6, 7, 8};

Physical Surface("BEM_domain_1", 1) = {1};

// ---- PARAMÈTRES ET BCs ----
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * epsilon_0);

SetNumber("Boundary Conditions/BEM_outside_E/BEM_domain_1/dirichlet", 1.0);
SetNumber("Boundary Conditions/BEM_Hole_H/BEM_domain_1/dirichlet", 0.3);

Mesh.Algorithm = 6;
