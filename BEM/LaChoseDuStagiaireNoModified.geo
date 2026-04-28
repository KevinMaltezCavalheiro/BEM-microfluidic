// ============================================================
//  Electrodes_Channel_Coupled_FIXED.geo
// ============================================================

scale = 1e-6;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL = 9.62;
eps0 = 8.8541878128e-12;

Lout_x = 1000*scale; Lout_y = 200*scale;
lc = 5*scale;

// ---- Géométrie : Tout en CCW naturel ----

// Carré Extérieur
Point(1) = {-Lout_x/2, -Lout_y/2, 0, lc};
Point(2) = { Lout_x/2, -Lout_y/2, 0, lc};
Point(3) = { Lout_x/2,  Lout_y/2, 0, lc};
Point(4) = {-Lout_x/2,  Lout_y/2, 0, lc};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};

// Electrode 1 (Trou)
cx1 = 350*scale; w1 = 100*scale; h1 = 60*scale;
Point(5) = {cx1-w1/2, -h1/2, 0, lc};
Point(6) = {cx1+w1/2, -h1/2, 0, lc};
Point(7) = {cx1+w1/2,  h1/2, 0, lc};
Point(8) = {cx1-w1/2,  h1/2, 0, lc};
Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8}; Line(8) = {8, 5};
Curve Loop(2) = {5, 6, 7, 8};

// Electrode 2 (Trou)
cx2 = -350*scale;
Point(9)  = {cx2-w1/2, -h1/2, 0, lc};
Point(10) = {cx2+w1/2, -h1/2, 0, lc};
Point(11) = {cx2+w1/2,  h1/2, 0, lc};
Point(12) = {cx2-w1/2,  h1/2, 0, lc};
Line(9) = {9, 10}; Line(10) = {10, 11}; Line(11) = {11, 12}; Line(12) = {12, 9};
Curve Loop(3) = {9, 10, 11, 12};

// Canal / Interface
Point(13) = {-150*scale, -30*scale, 0, lc};
Point(14) = { 150*scale, -30*scale, 0, lc};
Point(15) = { 150*scale,  30*scale, 0, lc};
Point(16) = {-150*scale,  30*scale, 0, lc};
Line(13) = {13, 14}; Line(14) = {14, 15}; Line(15) = {15, 16}; Line(16) = {16, 13};
Curve Loop(4) = {13, 14, 15, 16};

// ---- Surfaces ----
// Domaine 1 : PDMS est le carré MOINS les électrodes et MOINS le canal
Plane Surface(1) = {1, 2, 3, 4};

// Domaine 2 : Huile est l'INTÉRIEUR de la boucle 4
Plane Surface(2) = {4};

// ============================================================
// Groupes Physiques
// ============================================================
Physical Curve("BEM_outside", 100) = {1, 2, 3, 4};
Physical Curve("BEM_Hole_E1", 101) = {5, 6, 7, 8};
Physical Curve("BEM_Hole_E2", 102) = {9, 10, 11, 12};
Physical Curve("BEM_interface_1_2", 103) = {13, 14, 15, 16};

Physical Surface("BEM_domain_1", 1) = {1};
Physical Surface("BEM_domain_2", 2) = {2};

// ============================================================
// BCs
// ============================================================
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL * eps0);

SetNumber("Boundary Conditions/BEM_outside/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_Hole_E1/BEM_domain_1/dirichlet", 2000.0);
SetNumber("Boundary Conditions/BEM_Hole_E2/BEM_domain_1/dirichlet", 0.0);
// Pour l'interface, on ne met PAS de condition Dirichlet/Neumann, le solver couple les deux domaines
