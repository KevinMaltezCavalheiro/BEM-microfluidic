// ============================================================
// Thick_Electrodes_Holes.geo
// Électrodes considérées comme des trous dans le domaine BEM
// ============================================================
// NOMENCLATURE DES GROUPES PHYSIQUES :
// _E : Frontière extérieure
// _H : Trou / Hole (Frontière d'une électrode/trou dans un domaine)
// FORMAT : BEM_interface_{ID_A}_{ID_B}_{TYPE}
// ============================================================

// ============================================================
// PARAMÈTRES GÉOMÉTRIQUES & PHYSIQUES
// ============================================================

scale = 1e-6;
eps0  = 8.8541878128e-12;
epsilon_r = 1.0;

// Dimensions du problème
L_gap    = 400 * scale;  // Gap central entre électrodes (400 um)
W_elec   = 400 * scale;  // Largeur de chaque électrode (400 um)
H_elec   = 200 * scale;  // Épaisseur / Hauteur des électrodes (200 um)
L_pad    = 1000 * scale; // Isolation latérale de chaque côté (1000 um)

// Hauteurs du domaine fluide étendu
H_top    = 3000 * scale; // Extension du domaine vers le haut (3000 um)
H_bottom = 1000 * scale; // Extension du domaine vers le bas sous les électrodes (1000 um)

// Calcul des dimensions globales
L_total = 2 * L_pad + 2 * W_elec + L_gap;
x_min   = -L_total / 2;
x_max   =  L_total / 2;

// Tailles de maille
lc      = 20 * scale;    // Maille large en périphérie
lc_elec = 4 * scale;     // Maille fine aux abords des électrodes et dans le gap

// ============================================================
// 1. CONTOUR EXTÉRIEUR DU DOMAINE FLUIDE (Grande boîte)
// ============================================================

Point(1) = {x_min, -H_bottom, 0, lc};
Point(2) = {x_max, -H_bottom, 0, lc};
Point(3) = {x_max,  H_top,    0, lc};
Point(4) = {x_min,  H_top,    0, lc};

Line(1) = {1, 2}; // Bottom
Line(2) = {2, 3}; // Right
Line(3) = {3, 4}; // Top
Line(4) = {4, 1}; // Left

// Boucle du contour extérieur (sens anti-horaire)
Curve Loop(1) = {1, 2, 3, 4};

// ============================================================
// 2. ÉLECTRODE GAUCHE (TROU / HOLE 1) - Dir. = 1000 V
// ============================================================
// Positionnée entre y = 0 et y = H_elec
x1_left  = -L_gap/2 - W_elec;
x1_right = -L_gap/2;

Point(10) = {x1_left,  0,      0, lc_elec};
Point(11) = {x1_right, 0,      0, lc_elec};
Point(12) = {x1_right, H_elec, 0, lc_elec};
Point(13) = {x1_left,  H_elec, 0, lc_elec};

Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 10};

// Boucle du trou 1
Curve Loop(2) = {10, 11, 12, 13};

// ============================================================
// 3. ÉLECTRODE DROITE (TROU / HOLE 2) - Dir. = 0 V
// ============================================================
// Positionnée entre y = 0 et y = H_elec
x2_left  = L_gap/2;
x2_right = L_gap/2 + W_elec;

Point(20) = {x2_left,  0,      0, lc_elec};
Point(21) = {x2_right, 0,      0, lc_elec};
Point(22) = {x2_right, H_elec, 0, lc_elec};
Point(23) = {x2_left,  H_elec, 0, lc_elec};

Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 20};

// Boucle du trou 2
Curve Loop(3) = {20, 21, 22, 23};

// ============================================================
// 4. SURFACE DU DOMAINE FLUIDE
// ============================================================
// Domaine 1 = Contour Extérieur (1) MOINS Trou 1 (2) MOINS Trou 2 (3)
Plane Surface(1) = {1, 2, 3};

// ============================================================
// 5. GROUPES PHYSIQUES
// ============================================================

// Domaine principal
Physical Surface("BEM_domain_1", 1) = {1};

// Frontières extérieures (Boîte fermée avec Neumann = 0)
Physical Curve("BEM_outside_E", 100) = {1, 2, 3, 4};

// Trou Électrode Gauche (Dirichlet = 1000 V)
Physical Curve("BEM_hole_electrode_1_H", 101) = {10:13};

// Trou Électrode Droite (Dirichlet = 0 V)
Physical Curve("BEM_hole_electrode_2_H", 102) = {20:23};

// ============================================================
// 6. MATÉRIAUX ET CONDITIONS AUX LIMITES
// ============================================================

SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r * eps0);

// Conditions aux limites extérieures
SetNumber("Boundary Conditions/BEM_outside_E/BEM_domain_1/neumann", 0.0);

// Potentiels imposés sur les trous d'électrodes
SetNumber("Boundary Conditions/BEM_hole_electrode_1_H/BEM_domain_1/dirichlet", 1000.0);
SetNumber("Boundary Conditions/BEM_hole_electrode_2_H/BEM_domain_1/dirichlet", 0.0);
