// ============================================================
//  Minimal_BEM_1Domain_Hole.geo - CORRIGÉ
// ============================================================

scale = 1;
epsilon_r_PDMS = 2.8;

L = 1.0 * scale;
R = 0.2 * scale;
lc = 0.05 * scale;

// ---- Points carré ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// ---- Points cercle (SENS NATUREL ANTI-HORAIRE) ----
Point(5) = { 0,  0, 0, lc};
Point(6) = { R,  0, 0, lc};
Point(7) = { 0,  R, 0, lc}; // Retour au sens naturel
Point(8) = {-R,  0, 0, lc};
Point(9) = { 0, -R, 0, lc};

// ---- Lignes ----
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// ---- Loops & Surfaces ----
Curve Loop(1) = {1, 2, 3, 4};
// IMPORTANT : Pas de signe moins ici. Laisse Gmsh gérer la boucle.
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};

// ============================================================
// Groupes physiques
// ============================================================
Physical Curve("Left", 101)      = {4};
Physical Curve("Right", 102)     = {2};
Physical Curve("TopBottom", 103) = {1, 3};

// C'est ce NOM "BEM_Hole" qui déclenche l'inversion dans ton C++
Physical Curve("BEM_Hole", 104)  = {5, 6, 7, 8};

Physical Surface("BEM_domain_1", 9) = {1};

// ============================================================
// BCs
// ============================================================
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS);
SetNumber("Boundary Conditions/Left/BEM_domain_1/dirichlet",    1.0);
SetNumber("Boundary Conditions/Right/BEM_domain_1/dirichlet",   0.0);
SetNumber("Boundary Conditions/TopBottom/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_Hole/BEM_domain_1/dirichlet", 0.3);

Mesh.Algorithm = 6;
