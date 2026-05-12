// ============================================================
//  Rectangle_BEM_Simple_Validation_V2.geo
// ============================================================

scale = 1.0;
eps0 = 8.8541878128e-12;

Lx = 1.0;
Ly = 1.0;
lc = 0.05;

Point(1) = {0,  0,  0, lc};
Point(2) = {Lx, 0,  0, lc};
Point(3) = {Lx, Ly, 0, lc};
Point(4) = {0,  Ly, 0, lc};

Line(1) = {1, 2}; // Bas
Line(2) = {2, 3}; // Droite
Line(3) = {3, 4}; // Haut
Line(4) = {4, 1}; // Gauche

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Maillage Transfini
nx = 20;
ny = 20;
Transfinite Curve {1, 3} = nx + 1;
Transfinite Curve {2, 4} = ny + 1;
Transfinite Surface {1};

// --- NOUVEAUX GROUPES PHYSIQUES (Structure V2) ---
// On ajoute le suffixe _E pour indiquer qu'il s'agit de frontières extérieures
Physical Curve("BEM_left_E", 1) = {4};
Physical Curve("BEM_right_E", 2) = {2};
Physical Curve("BEM_top_bottom_E", 3) = {1, 3};

Physical Surface("BEM_domain_1", 10) = {1};

// --- Paramètres ---
SetNumber("Materials/BEM_domain_1/Epsilon", 1.0 * eps0);
SetNumber("Non_linear_solver", 0);

// --- BCs mises à jour ---
SetNumber("Boundary Conditions/BEM_left_E/BEM_domain_1/dirichlet", 200.0);
SetNumber("Boundary Conditions/BEM_right_E/BEM_domain_1/dirichlet", 100.0);
SetNumber("Boundary Conditions/BEM_top_bottom_E/BEM_domain_1/neumann", 0.0);
