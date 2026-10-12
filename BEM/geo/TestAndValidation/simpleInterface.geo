// ============================================================
//  Minimal_BEM_Debug_V2.geo
//  But : Tester la continuité et le saut de flux (Inclusion)
// ============================================================
// NOMENCLATURE DES GROUPES PHYSIQUES :
// ------------------------------------
// _E    : Frontière Extérieure (Side A unique, normale directe)
// _inc  : Interface Incluse (Side A entoure Side B, normale inversée pour A)
//
// FORMAT : BEM_{NOM}_{TYPE} ou BEM_interface_{ID_A}_{ID_B}_{TYPE}
// ============================================================

scale = 1;
epsilon_0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL  = 80; // On met un fort contraste pour bien voir la déviation

// ---- Paramètres simples ----
L = 1.0 * scale;
R = 0.2 * scale;
lc = 0.02 * scale;

// ---- Points Carré Extérieur (Domaine 1) ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// ---- Points Cercle Intérieur (Interface) ----
Point(5) = {0,  0, 0, lc}; // Centre
Point(6) = {R,  0, 0, lc};
Point(7) = {0,  R, 0, lc};
Point(8) = {-R, 0, 0, lc};
Point(9) = {0, -R, 0, lc};

// ---- Lignes ----
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cercle
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// ---- Loops & Surfaces ----
Curve Loop(1) = {1, 2, 3, 4}; // Contour extérieur
Curve Loop(2) = {5, 6, 7, 8}; // Contour inclusion

Plane Surface(1) = {1, 2};    // Domaine 1 (PDMS, percé)
Plane Surface(2) = {2};       // Domaine 2 (Huile/Goutte)

// ---- GROUPES PHYSIQUES V2 ----

// Frontières extérieures (_E)
Physical Curve("BEM_Left_E", 101)  = {4};
Physical Curve("BEM_Right_E", 102) = {2};
Physical Curve("BEM_TopBottom_E", 103) = {1, 3};

// Interface Incluse (1 = PDMS entoure 2 = OIL)
Physical Curve("BEM_interface_1_2_inc", 8) = {5, 6, 7, 8};

// Domaines
Physical Surface("BEM_domain_1", 1)  = {1};
Physical Surface("BEM_domain_2", 2) = {2};

// ---- Propriétés et BCs ----
// Utilisation de la permittivité absolue (eps_r * eps_0)
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * epsilon_0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL  * epsilon_0);

SetNumber("Non_linear_solver", 0);

// BCs indexées sur les nouveaux noms
SetNumber("Boundary Conditions/BEM_Left_E/BEM_domain_1/dirichlet", 1.0);
SetNumber("Boundary Conditions/BEM_Right_E/BEM_domain_1/dirichlet", 0.0);
SetNumber("Boundary Conditions/BEM_TopBottom_E/BEM_domain_1/neumann", 0.0);
