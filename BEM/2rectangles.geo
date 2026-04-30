// ============================================================
//  Validation_BEM_Full_Logic_V2.geo
// ============================================================
// NOMENCLATURE DES GROUPES PHYSIQUES :
// ------------------------------------
// _E    : Frontière Extérieure (Side A unique, normale directe)
// _jux  : Interface Juxtaposée (Side A & B, domaines côte à côte)
// _inc  : Interface Incluse (Side A entoure Side B, normale inversée pour A)
// _H    : Trou / Hole (Side A entoure un vide, normale inversée pour A)
//
// FORMAT : BEM_interface_{ID_A}_{ID_B}_{TYPE}
// ============================================================

scale = 1.0;
eps0 = 8.8541878128e-12;
lc = 0.04;

// --- Points pour les blocs ---
Point(1) = {0, 0, 0, lc}; Point(2) = {1, 0, 0, lc}; Point(3) = {2.5, 0, 0, lc};
Point(4) = {2.5, 1, 0, lc}; Point(5) = {1, 1, 0, lc}; Point(6) = {0, 1, 0, lc};

// --- Cercle (Inclusion dans le Domaine 2) ---
cx = 1.75; cy = 0.3; r = 0.15;
Point(10) = {cx, cy, 0, lc};
Point(11) = {cx + r, cy, 0, lc};
Point(12) = {cx, cy + r, 0, lc};
Point(13) = {cx - r, cy, 0, lc};
Point(14) = {cx, cy - r, 0, lc};

Circle(20) = {11, 10, 12}; Circle(21) = {12, 10, 13};
Circle(22) = {13, 10, 14}; Circle(23) = {14, 10, 11};

// --- Carré (Trou dans le Domaine 2) ---
tx = 1.6; ty = 0.7; ts = 0.2;
Point(30) = {tx, ty, 0, lc}; Point(31) = {tx+ts, ty, 0, lc};
Point(32) = {tx+ts, ty+ts, 0, lc}; Point(33) = {tx, ty+ts, 0, lc};

Line(30) = {30, 31}; Line(31) = {31, 32};
Line(32) = {32, 33}; Line(33) = {33, 30};

// --- Lignes des domaines ---
Line(1) = {1, 2}; Line(2) = {2, 5}; Line(3) = {5, 6}; Line(4) = {6, 1}; // Bloc 1
Line(5) = {2, 3}; Line(6) = {3, 4}; Line(7) = {4, 5};                  // Bloc 2

// --- Surfaces ---
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; // Domaine 1 (Gauche)

Curve Loop(2) = {5, 6, 7, -2};    // Contour extérieur du Domaine 2
Curve Loop(3) = {20, 21, 22, 23}; // Contour de la goutte
Curve Loop(4) = {30, 31, 32, 33}; // Contour du trou
Plane Surface(2) = {2, 3, 4};     // Domaine 2 (percé par goutte et trou)

Curve Loop(5) = {20, 21, 22, 23};
Plane Surface(3) = {5};           // Domaine 3 (La Goutte)

// --- GROUPES PHYSIQUES ---
Physical Surface("BEM_domain_1", 1) = {1};
Physical Surface("BEM_domain_2", 2) = {2};
Physical Surface("BEM_domain_3", 3) = {3};

// Frontières Extérieures
Physical Curve("BEM_left_E", 100) = {4};
Physical Curve("BEM_right_E", 101) = {6};
Physical Curve("BEM_top_bottom_E", 102) = {1, 3, 5, 7};

// Interfaces
Physical Curve("BEM_interface_1_2_jux", 103) = {2};      // Juxtaposition
Physical Curve("BEM_interface_2_3_inc", 104) = {20:23}; // Inclusion goutte

// Trous (Holes)
Physical Curve("BEM_hole_electrode_H", 105) = {30:33};

// --- Paramètres ---
SetNumber("Materials/BEM_domain_1/Epsilon", 1.0 * eps0);
SetNumber("Materials/BEM_domain_2/Epsilon", 2.0 * eps0);
SetNumber("Materials/BEM_domain_3/Epsilon", 80.0 * eps0);

SetNumber("Non_linear_solver", 0);

// --- BCs ---
SetNumber("Boundary Conditions/BEM_left_E/BEM_domain_1/dirichlet", 200.0);
SetNumber("Boundary Conditions/BEM_right_E/BEM_domain_2/dirichlet", 100.0);
SetNumber("Boundary Conditions/BEM_top_bottom_E/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_top_bottom_E/BEM_domain_2/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_hole_electrode_H/BEM_domain_2/dirichlet", 150.0);
