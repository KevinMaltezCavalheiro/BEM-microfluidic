// ============================================================
//  Electrodes_Channel_T_Shape_Two_Particles_V2.geo
// ============================================================

scale = 1e-6;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL = 9.62;
epsilon_r_WATER = 80.0;
eps0 = 8.8541878128e-12;

Lout_x = 1000*scale; Lout_y = 200*scale;
lc = 2*scale;
lc_part = 1*scale;

// ---- 1. CADRE EXTÉRIEUR ----
Point(1) = {-Lout_x/2, -Lout_y/2, 0, lc};
Point(2) = { Lout_x/2, -Lout_y/2, 0, lc};
Point(3) = { Lout_x/2,  Lout_y/2, 0, lc};
Point(4) = {-Lout_x/2,  Lout_y/2, 0, lc};

Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};

// ---- 2. ÉLECTRODES (Trous dans le PDMS) ----
cx1 = 350*scale; w1 = 100*scale; h1 = 60*scale;
Point(5) = {cx1-w1/2, -h1/2, 0, lc}; Point(6) = {cx1+w1/2, -h1/2, 0, lc};
Point(7) = {cx1+w1/2,  h1/2, 0, lc}; Point(8) = {cx1-w1/2,  h1/2, 0, lc};
Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8}; Line(8) = {8, 5};

cx2 = -350*scale;
Point(9)  = {cx2-w1/2, -h1/2, 0, lc}; Point(10) = {cx2+w1/2, -h1/2, 0, lc};
Point(11) = {cx2+w1/2,  h1/2, 0, lc}; Point(12) = {cx2-w1/2,  h1/2, 0, lc};
Line(9) = {9, 10}; Line(10) = {10, 11}; Line(11) = {11, 12}; Line(12) = {12, 9};

// ---- 3. CANAL EN T (Inclusion dans le PDMS) ----
cx3 = 0.0; Lx3 = 300.0*scale; Ly3 = 60.0*scale; w_cheminee = 80*scale;
Point(13) = {cx3 - Lx3/2,      -Ly3/2, 0, lc}; Point(14) = {cx3 + Lx3/2,      -Ly3/2, 0, lc};
Point(15) = {cx3 + Lx3/2,       0,      0, lc}; Point(16) = {cx3 + w_cheminee/2, 0,      0, lc};
Point(17) = {cx3 + w_cheminee/2, Ly3/2,  0, lc}; Point(18) = {cx3 - w_cheminee/2, Ly3/2,  0, lc};
Point(19) = {cx3 - w_cheminee/2, 0,      0, lc}; Point(20) = {cx3 - Lx3/2,       0,      0, lc};
Line(13) = {13, 14}; Line(14) = {14, 15}; Line(15) = {15, 16}; Line(16) = {16, 17};
Line(17) = {17, 18}; Line(18) = {18, 19}; Line(19) = {19, 20}; Line(20) = {20, 13};

// ---- 4. PARTICULES (Inclusions dans l'Huile) ----
p1_rad = 20*scale; p1_cx = -19*scale; p1_cy = 9*scale;
Point(21) = {p1_cx, p1_cy, 0, lc_part}; Point(22) = {p1_cx + p1_rad, p1_cy, 0, lc_part};
Point(23) = {p1_cx, p1_cy + p1_rad, 0, lc_part}; Point(24) = {p1_cx - p1_rad, p1_cy, 0, lc_part};
Point(25) = {p1_cx, p1_cy - p1_rad, 0, lc_part};
Circle(21) = {22, 21, 23}; Circle(22) = {23, 21, 24}; Circle(23) = {24, 21, 25}; Circle(24) = {25, 21, 22};

p2_rad = 20*scale; p2_cx = 19.5*scale; p2_cy = -8.5*scale;
Point(26) = {p2_cx, p2_cy, 0, lc_part}; Point(27) = {p2_cx + p2_rad, p2_cy, 0, lc_part};
Point(28) = {p2_cx, p2_cy + p2_rad, 0, lc_part}; Point(29) = {p2_cx - p2_rad, p2_cy, 0, lc_part};
Point(30) = {p2_cx, p2_cy - p2_rad, 0, lc_part};
Circle(25) = {27, 26, 28}; Circle(26) = {28, 26, 29}; Circle(27) = {29, 26, 30}; Circle(28) = {30, 26, 27};

// ---- SURFACES ----
Curve Loop(1) = {1, 2, 3, 4};        // Cadre extérieur
Curve Loop(2) = {5, 6, 7, 8};        // Électrode 1
Curve Loop(3) = {9, 10, 11, 12};     // Électrode 2
Curve Loop(4) = {13, 14, 15, 16, 17, 18, 19, 20}; // Canal

Plane Surface(1) = {1, 2, 3, 4};     // PDMS percé par E1, E2 et le Canal

Curve Loop(5) = {21, 22, 23, 24};    // Particule 1
Curve Loop(6) = {25, 26, 27, 28};    // Particule 2
Plane Surface(2) = {4, 5, 6};        // HUILE percée par P1 et P2
Plane Surface(3) = {5};              // EAU 1
Plane Surface(4) = {6};              // EAU 2

// ---- GROUPES PHYSIQUES (NOUVELLE STRUCTURE) ----
Physical Curve("BEM_outside_E", 100) = {1, 2, 3, 4};
Physical Curve("BEM_Hole_E1_H", 101) = {5, 6, 7, 8};
Physical Curve("BEM_Hole_E2_H", 102) = {9, 10, 11, 12};

// Interface Inclusion : Side A = PDMS (1), Side B = HUILE (2)
Physical Curve("BEM_interface_1_2_inc", 103) = {13:20};

// Interface Inclusion : Side A = HUILE (2), Side B = EAU (3 ou 4)
Physical Curve("BEM_interface_2_3_inc", 104) = {21:24};
Physical Curve("BEM_interface_2_4_inc", 105) = {25:28};

Physical Surface("BEM_domain_1", 1) = {1}; // PDMS
Physical Surface("BEM_domain_2", 2) = {2}; // OIL
Physical Surface("BEM_domain_3", 3) = {3}; // WATER 1
Physical Surface("BEM_domain_4", 4) = {4}; // WATER 2

// ---- PARAMÈTRES ET BCs ----
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL * eps0);
SetNumber("Materials/BEM_domain_3/Epsilon", epsilon_r_WATER * eps0);
SetNumber("Materials/BEM_domain_4/Epsilon", epsilon_r_WATER * eps0);

SetNumber("Boundary Conditions/BEM_outside_E/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_Hole_E1_H/BEM_domain_1/dirichlet", 2000.0);
SetNumber("Boundary Conditions/BEM_Hole_E2_H/BEM_domain_1/dirichlet", 0.0);
