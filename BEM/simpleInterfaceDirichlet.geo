// ============================================================
//  Debug_BEM_NoInterface.geo – Carré avec Trou Circulaire (Dirichlet)
//  But : Tester le Domaine 1 seul, sans couplage d'interface.
// ============================================================

scale = 1e-6;
epsilon_0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;

// ---- Paramètres ----
L = 1000.0 * scale;
R = 200.0 * scale;
lc = 25.0 * scale;

// ---- Points Carré Extérieur ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// ---- Points Cercle Intérieur (Le trou) ----
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

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// ---- Loops & Surfaces ----
Curve Loop(1) = {1, 2, 3, 4}; // Extérieur
Curve Loop(2) = {5, 6, 7, 8}; // Le trou

// Surface 1 : On ne définit QUE le domaine 1 (le PDMS)
// Le cercle est soustrait (trou) car il est listé dans les boucles de la surface
Plane Surface(1) = {1, 2};

// ---- Groupes Physiques ----

// Un seul domaine BEM ici
Physical Surface("BEM_domain_1", 9) = {1};

// Bord extérieur
Physical Curve("outside", 5) = {1, 2, 3, 4};

// Le cercle n'est plus une "BEM_interface", c'est maintenant un simple "hole"
Physical Curve("hole1", 8) = {5, 6, 7, 8};

// ---- Propriétés et BCs ----
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * epsilon_0);

// On impose une différence de potentiel pour voir si le champ se calcule bien
SetNumber("Boundary Conditions/outside/BEM_domain_1/dirichlet", 1000.0);
SetNumber("Boundary Conditions/hole1/BEM_domain_1/dirichlet", 0.0);
