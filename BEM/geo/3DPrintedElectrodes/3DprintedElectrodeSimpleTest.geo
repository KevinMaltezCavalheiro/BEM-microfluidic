// ============================================================
//  Coplanar_Electrodes_Validation.geo
//  Électrodes planes sur la frontière inférieure
// ============================================================

scale = 1e-6;
//epsilon_r = 2.8; // Ex: PDMS ou milieu fluide
epsilon_r = 1; // Ex: PDMS ou milieu fluide
eps0 = 8.8541878128e-12;

// ---- DIMENSIONS ----
L_gap   = 400 * scale; // Largeur du gap central entre électrodes
W_elec  = 400 * scale; // Largeur de chaque électrode
L_pad   = 1000 * scale; // Marge d'isolation sur les côtés pour éviter les effets de bord
H_domain= 3000 * scale; // Hauteur du domaine de calcul

// Calcul de la largeur totale : L_pad + W_elec + L_gap + W_elec + L_pad
L_total = 2 * L_pad + 2 * W_elec + L_gap;
x_min   = -L_total / 2;

lc = 15 * scale;       // Taille de maille générale (loin des zones critiques)
lc_elec = 5 * scale;   // Taille de maille très fine au niveau des électrodes et du gap

// ---- 1. POINTS DE LA FRONTIÈRE INFÉRIEURE (y = 0) ----
// De la gauche vers la droite
Point(1) = {x_min,                             0, 0, lc};
Point(2) = {x_min + L_pad,                     0, 0, lc_elec};
Point(3) = {x_min + L_pad + W_elec,            0, 0, lc_elec};
Point(4) = {x_min + L_pad + W_elec + L_gap,    0, 0, lc_elec};
Point(5) = {x_min + L_pad + 2*W_elec + L_gap,  0, 0, lc_elec};
Point(6) = {x_min + L_total,                   0, 0, lc};

// ---- 2. POINTS SUPÉRIEURS ET LATÉRAUX ----
Point(7) = {x_min + L_total, H_domain, 0, lc};
Point(8) = {x_min,           H_domain, 0, lc};

// ---- 3. LIGNES ----
// Bords inférieurs découpés :
Line(1) = {1, 2}; // Isolation gauche
Line(2) = {2, 3}; // Électrode 1 (V1)
Line(3) = {3, 4}; // Gap central (Neumann / Isolat)
Line(4) = {4, 5}; // Électrode 2 (V2)
Line(5) = {5, 6}; // Isolation droite

// Bords extérieur (côtés et haut) :
Line(6) = {6, 7}; // Côté droit
Line(7) = {7, 8}; // Haut
Line(8) = {8, 1}; // Côté gauche

// ---- 4. SURFACES ----
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// ---- 5. GROUPES PHYSIQUES ----
// Électrodes planes (Dirichlet)
Physical Curve("BEM_Electrode_1_E", 100) = {2};
Physical Curve("BEM_Electrode_2_E", 101) = {4};

// Bords isolants extérieurs et gap central (Neumann, dn V = 0)
Physical Curve("BEM_outside_top_E", 102) = {1, 3, 5, 6, 7, 8};

// Domaine physique
Physical Surface("BEM_domain_1", 1) = {1};

// ---- PARAMÈTRES ET BCs ----
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r * eps0);

// Conditions aux limites
SetNumber("Boundary Conditions/BEM_Electrode_1_E/BEM_domain_1/dirichlet", 1000.0);
SetNumber("Boundary Conditions/BEM_Electrode_2_E/BEM_domain_1/dirichlet", 0.0);
SetNumber("Boundary Conditions/BEM_outside_top_E/BEM_domain_1/neumann", 0.0);
