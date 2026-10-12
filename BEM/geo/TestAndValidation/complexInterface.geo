// ============================================================
//  BEM_MultiDomain_With_Hole_V2.geo
// ============================================================
// NOMENCLATURE DES GROUPES PHYSIQUES :
// ------------------------------------
// _E    : Frontière Extérieure (Side A unique, normale directe)
// _inc  : Interface Incluse (Side A entoure Side B, normale inversée pour A)
// _H    : Trou / Hole (Side A entoure un vide, normale inversée pour A)
//
// FORMAT : BEM_interface_{ID_A}_{ID_B}_{TYPE} ou BEM_{NOM}_{TYPE}
// ============================================================

scale = 1;
epsilon_0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL  = 80;

L = 1.0 * scale;
R_oil = 0.2 * scale;
R_hole = 0.1 * scale;
lc = 0.01 * scale; // Légèrement augmenté pour la clarté du debug

// ---- Points Carré Extérieur (Domaine 1) ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// ---- Points Cercle Huile (Interface 1-2) ----
Point(5) = {-0.2, 0.2, 0, lc};
Point(6) = {-0.2+R_oil, 0.2, 0, lc};
Point(7) = {-0.2, 0.2+R_oil, 0, lc};
Point(8) = {-0.2-R_oil, 0.2, 0, lc};
Point(9) = {-0.2, 0.2-R_oil, 0, lc};

// ---- Points Cercle Trou (Hole dans Domaine 1) ----
Point(10) = {0.2, -0.2, 0, lc};
Point(11) = {0.2+R_hole, -0.2, 0, lc};
Point(12) = {0.2, -0.2+R_hole, 0, lc};
Point(13) = {0.2-R_hole, -0.2, 0, lc};
Point(14) = {0.2, -0.2-R_hole, 0, lc};

// ---- Lignes ----
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};

Circle(5) = {6, 5, 7}; Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9}; Circle(8) = {9, 5, 6};

Circle(9) = {11, 10, 12}; Circle(10) = {12, 10, 13};
Circle(11) = {13, 10, 14}; Circle(12) = {14, 10, 11};

// ---- Loops & Surfaces ----
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {9, 10, 11, 12};

Plane Surface(1) = {1, 2, 3}; // Domaine 1 (PDMS)
Plane Surface(2) = {2};       // Domaine 2 (Huile)

// ============================================================
// GROUPES PHYSIQUES (Structure V2)
// ============================================================

// Frontières Externes (_E)
Physical Curve("BEM_Left_E", 101)      = {4};
Physical Curve("BEM_Right_E", 102)     = {2};
Physical Curve("BEM_TopBottom_E", 103) = {1, 3};

// Le Trou (_H) - Sera traité comme une électrode interne
Physical Curve("BEM_Hole_Electrode_H", 104) = {9, 10, 11, 12};

// L'interface Incluse (_inc) entre PDMS (1) et Huile (2)
Physical Curve("BEM_interface_1_2_inc", 105) = {5, 6, 7, 8};

// Domaines
Physical Surface("BEM_domain_1", 1) = {1};
Physical Surface("BEM_domain_2", 2) = {2};

// ============================================================
// MATÉRIAUX & BCs
// ============================================================
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * epsilon_0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL * epsilon_0);

SetNumber("Non_linear_solver", 0);

// Conditions aux limites
SetNumber("Boundary Conditions/BEM_Left_E/BEM_domain_1/dirichlet", 1.0);
SetNumber("Boundary Conditions/BEM_Right_E/BEM_domain_1/dirichlet", 0.0);
SetNumber("Boundary Conditions/BEM_TopBottom_E/BEM_domain_1/neumann", 0.0);

// Potentiel imposé sur l'électrode circulaire interne
SetNumber("Boundary Conditions/BEM_Hole_Electrode_H/BEM_domain_1/dirichlet", 0.5);

Mesh.Algorithm = 6; // Frontal-Delaunay
