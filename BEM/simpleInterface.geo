// ============================================================
//  Minimal_BEM_Debug.geo – Carré (PDMS) avec Cercle (Oil)
//  But : Tester la continuité à l'interface 1_2
// ============================================================

scale = 1;
epsilon_0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL  = 80;

// ---- Paramètres simples ----
L = 1.0 * scale; // Côté du carré
R = 0.2 * scale;  // Rayon du cercle central
lc = 0.01 * scale;  // Taille des mailles

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
// Carré extérieur
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Cercle (Interface)
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// ---- Loops & Surfaces ----
Curve Loop(1) = {1, 2, 3, 4}; // Extérieur
Curve Loop(2) = {5, 6, 7, 8}; // Interface

// Surface 1 : Le "cadre" (PDMS)
Plane Surface(1) = {1, 2};

// Surface 2 : Le disque (Oil)
Plane Surface(2) = {2};

// ---- Remplace tes Groupes Physiques actuels par ceux-ci ----

// Côtés individuels pour BCs différenciées
Physical Curve("Left", 101)  = {4}; // Bord gauche (phi = 1)
Physical Curve("Right", 102) = {2}; // Bord droit (phi = 0)
Physical Curve("TopBottom", 103) = {1, 3}; // Haut et Bas (Isolants ou qn=0)

// L'interface reste inchangée
Physical Curve("BEM_interface_1_2", 8) = {5, 6, 7, 8};

// Domaines restent inchangés
Physical Surface("BEM_domain_1", 9)  = {1};
Physical Surface("BEM_domain_2", 10) = {2};

// ---- Propriétés et BCs ----
//SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * epsilon_0);
//SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL  * epsilon_0);
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL);

// ---- Nouvelles BCs pour tester le gradient ----

// Bord Gauche : 10.0 Volt
SetNumber("Boundary Conditions/Left/BEM_domain_1/dirichlet", 1.0);

// Bord Droit : 0.0 Volt
SetNumber("Boundary Conditions/Right/BEM_domain_1/dirichlet", 0.0);

// Haut et Bas : Flux nul (Symétrie/Isolant)
// Si ton code ne gère pas encore 'neumann', utilise Dirichlet 0.5 (moyenne)
// ou laisse-les en Dirichlet 1.0 pour voir si le champ dévie.
SetNumber("Boundary Conditions/TopBottom/BEM_domain_1/neumann", 0);

//Mesh.Tangents = 50; // Set pixel size for 1D orientation vectors
