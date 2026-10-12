// ============================================================
//  Electrodes_Homogeneous_Validation.geo
// ============================================================

scale = 1e-6;
//epsilon_r_PDMS   = 2.8;
epsilon_r_PDMS   = 4.0;
epsilon_r_RESIN  = 4.0;
eps0 = 8.8541878128e-12;

Lout_x = 1000*scale; Lout_y = 6000*scale;
H_interface = -120*scale;
lc = 10*scale;

// ---- 1. CADRE EXTÉRIEUR ----
Point(1) = {-Lout_x/2, -Lout_y/2, 0, lc};
Point(2) = { Lout_x/2, -Lout_y/2, 0, lc};
Point(3) = { Lout_x/2,  Lout_y/3, 0, lc};
Point(4) = {-Lout_x/2,  Lout_y/3, 0, lc};
Point(100) = {-Lout_x/2, H_interface, 0, lc};
Point(101) = { Lout_x/2, H_interface, 0, lc};

Line(1) = {1, 2}; Line(2) = {2, 101}; Line(3) = {101, 3};
Line(4) = {3, 4}; Line(5) = {4, 100}; Line(6) = {100, 1};
Line(7) = {100, 101}; // INTERFACE PLATE

// ---- 2. ÉLECTRODES (Trous dans la Résine) ----
cx1 = 350*scale; w1 = 100*scale; h1 = 60*scale; y_elec = -151*scale;
Point(5) = {cx1-w1/2, y_elec-h1/2, 0, lc}; Point(6) = {cx1+w1/2, y_elec-h1/2, 0, lc};
Point(7) = {cx1+w1/2, y_elec+h1/2, 0, lc}; Point(8) = {cx1-w1/2, y_elec+h1/2, 0, lc};
Line(8) = {5, 6}; Line(9) = {6, 7}; Line(10) = {7, 8}; Line(11) = {8, 5};

cx2 = -350*scale;
Point(9)  = {cx2-w1/2, y_elec-h1/2, 0, lc}; Point(10) = {cx2+w1/2, y_elec-h1/2, 0, lc};
Point(11) = {cx2+w1/2, y_elec+h1/2, 0, lc}; Point(12) = {cx2-w1/2, y_elec+h1/2, 0, lc};
Line(12) = {9, 10}; Line(13) = {10, 11}; Line(14) = {11, 12}; Line(15) = {12, 9};

// ---- SURFACES ----
Curve Loop(1) = {7, 3, 4, 5};
Plane Surface(1) = {1};        // Domaine supérieur homogène (ex-PDMS, sans canal)

Curve Loop(2) = {1, 2, -7, 6};
Curve Loop(3) = {8:11};
Curve Loop(4) = {12:15};
Plane Surface(5) = {2, 3, 4};  // RÉSINE (percée par électrodes)

// ---- GROUPES PHYSIQUES AVEC CONVENTION DE NOMMAGE ----
Physical Curve("BEM_outside_top_E", 100) = {3, 4, 5};
Physical Curve("BEM_outside_bottom_E", 101) = {1, 2, 6};

// Interface Juxtaposée (Side A: 5, Side B: 1)
Physical Curve("BEM_interface_1_5_jux", 102) = {7};

// Trous Électrodes (H)
Physical Curve("BEM_Hole_E1_H", 106) = {8:11};
Physical Curve("BEM_Hole_E2_H", 107) = {12:15};

Physical Surface("BEM_domain_1", 1) = {1}; // Milieu supérieur homogène
Physical Surface("BEM_domain_5", 5) = {5}; // RÉSINE

// ---- PARAMÈTRES ET BCs ----
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_5/Epsilon", epsilon_r_RESIN * eps0);

SetNumber("Boundary Conditions/BEM_outside_top_E/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_outside_bottom_E/BEM_domain_5/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_Hole_E1_H/BEM_domain_5/dirichlet", 2000.0);
SetNumber("Boundary Conditions/BEM_Hole_E2_H/BEM_domain_5/dirichlet", 0.0);
