// ============================================================
//  BEM_MultiDomain_With_Hole.geo
// ============================================================

scale = 1;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL  = 80;

L = 1.0 * scale;
R_oil = 0.2 * scale;  // L'interface (huile)
R_hole = 0.1 * scale; // Le trou (vide)
lc = 0.01 * scale;

// ---- Points Carré Extérieur (Domaine 1) ----
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// ---- Points Cercle Huile (Interface 1-2) ----
Point(5) = {-0.2, 0.2, 0, lc}; // Centre décalé pour l'huile
Point(6) = {-0.2+R_oil, 0.2, 0, lc};
Point(7) = {-0.2, 0.2+R_oil, 0, lc};
Point(8) = {-0.2-R_oil, 0.2, 0, lc};
Point(9) = {-0.2, 0.2-R_oil, 0, lc};

// ---- Points Cercle Trou (Hole dans Domaine 1) ----
Point(10) = {0.2, -0.2, 0, lc}; // Centre décalé pour le trou
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
Curve Loop(1) = {1, 2, 3, 4};  // Frontière externe
Curve Loop(2) = {5, 6, 7, 8};  // Interface Huile
Curve Loop(3) = {9, 10, 11, 12}; // Trou

// Surface 1 : PDMS (Cadre avec deux découpes)
Plane Surface(1) = {1, 2, 3};

// Surface 2 : Huile (Disque plein)
Plane Surface(2) = {2};

// ============================================================
// Groupes Physiques (La clé pour ton code C++)
// ============================================================

// Bords Externes (interfaceSide = 'E', normalOutward = true)
Physical Curve("Left", 101)      = {4};
Physical Curve("Right", 102)     = {2};
Physical Curve("TopBottom", 103) = {1, 3};

// Le Trou (Ton code va détecter "Hole" -> interfaceSide = 'H' -> flip normale)
Physical Curve("BEM_Hole_1", 104) = {9, 10, 11, 12};

// L'interface (interfaceSide = 'A'/'B' -> gestion automatique du signe)
Physical Curve("BEM_interface_1_2", 105) = {5, 6, 7, 8};

// Domaines
Physical Surface("BEM_domain_1", 9)  = {1};
Physical Surface("BEM_domain_2", 10) = {2};

// ============================================================
// Matériaux & BCs
// ============================================================
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL);

SetNumber("Boundary Conditions/Left/BEM_domain_1/dirichlet", 1.0);
SetNumber("Boundary Conditions/Right/BEM_domain_1/dirichlet", 0.0);
SetNumber("Boundary Conditions/TopBottom/BEM_domain_1/neumann", 0.0);

// On impose un potentiel sur le trou (ex: une électrode interne)
SetNumber("Boundary Conditions/BEM_Hole_1/BEM_domain_1/dirichlet", 0.5);

Mesh.Algorithm = 6; // Frontal-Delaunay
