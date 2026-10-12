// ============================================================
//  Electrodes_Channel_MultiWell_PDMS_Glass_FC40_Resin_v7.geo
// ============================================================

scale = 1e-6;

epsilon_r_PDMS  = 2.8;
epsilon_r_RESIN = 4.0;
epsilon_r_FC40  = 1.9;  // Huile FC-40 entre Résine et Verre
epsilon_r_GLASS = 5.0;
epsilon_r_OIL   = 9.62; // Huile du canal
epsilon_r_WATER = 80.0;

eps0 = 8.8541878128e-12;

Lout_x = 1000*scale;
Lout_y = 600*scale;

// --- Niveaux en Y des sous-couches ---
H_interface_resin_fc40 = -120*scale; // Interface RESIN / FC-40
H_fc40_thickness       = 3*scale;    // Épaisseur FC-40 (3 µm)
H_interface_fc40_glass = H_interface_resin_fc40 + H_fc40_thickness; // -117 µm

H_interface_glass_pdms = -35*scale;  // Interface GLASS / PDMS (-35 µm)

// Tailles de maillage relaxées pour accélérer la résolution
lc       = 12*scale;
lc_thin  = 5*scale;  // Couche FC-40 allégée
lc_part  = 3*scale;  // Gouttes
lc_canal = 4*scale;  // Canal & Puits


// ---- 1. CADRE EXTÉRIEUR ET STRATES HORIZONTAILES ----

Point(1) = {-Lout_x/2, -Lout_y/2, 0, lc};
Point(2) = { Lout_x/2, -Lout_y/2, 0, lc};
Point(3) = { Lout_x/2,  Lout_y/3, 0, lc};
Point(4) = {-Lout_x/2,  Lout_y/3, 0, lc};

// Points d'interfaces latérales
Point(100) = {-Lout_x/2, H_interface_resin_fc40, 0, lc_thin};
Point(101) = { Lout_x/2, H_interface_resin_fc40, 0, lc_thin};

Point(104) = {-Lout_x/2, H_interface_fc40_glass, 0, lc_thin};
Point(105) = { Lout_x/2, H_interface_fc40_glass, 0, lc_thin};

Point(102) = {-Lout_x/2, H_interface_glass_pdms, 0, lc};
Point(103) = { Lout_x/2, H_interface_glass_pdms, 0, lc};


// RESIN - contour
Line(1) = {1, 2};
Line(2) = {2, 101};
Line(6) = {100, 1};

// FC-40 - côté droit & gauche
Line(48) = {101, 105};
Line(49) = {104, 100};

// GLASS - côté droit & gauche
Line(3) = {105, 103};
Line(7) = {102, 104};

// PDMS - côté droit & gauche
Line(32) = {103, 3};
Line(5)  = {4, 102};

// Top
Line(4)  = {3, 4};


// Interfaces horizontales orientées vers la droite (Gauche -> Droite)
Line(35) = {100, 101}; // Gauche = Domaine 8 (FC40), Droite = Domaine 5 (RESIN)
Line(50) = {104, 105}; // Gauche = Domaine 6 (GLASS), Droite = Domaine 8 (FC40)
Line(34) = {102, 103}; // Gauche = Domaine 1 (PDMS), Droite = Domaine 6 (GLASS)


// ============================================================
// 2. ÉLECTRODES (Trous dans la Résine)
// ============================================================

cx1 = 150*scale;
w1 = 100*scale;
h1 = 60*scale;
y_elec = -151*scale;

Point(5) = {cx1-w1/2, y_elec-h1/2, 0, lc};
Point(6) = {cx1+w1/2, y_elec-h1/2, 0, lc};
Point(7) = {cx1+w1/2, y_elec+h1/2, 0, lc};
Point(8) = {cx1-w1/2, y_elec+h1/2, 0, lc};

Line(8)  = {5, 6};
Line(9)  = {6, 7};
Line(10) = {7, 8};
Line(11) = {8, 5};


cx2 = -150*scale;

Point(9)  = {cx2-w1/2, y_elec-h1/2, 0, lc};
Point(10) = {cx2+w1/2, y_elec-h1/2, 0, lc};
Point(11) = {cx2+w1/2, y_elec+h1/2, 0, lc};
Point(12) = {cx2-w1/2, y_elec+h1/2, 0, lc};

Line(12) = {9, 10};
Line(13) = {10, 11};
Line(14) = {11, 12};
Line(15) = {12, 9};


// ============================================================
// 3. CANAL (30 µm) ET PUITS (30 µm SUPPLÉMENTAIRES)
// ============================================================

Lx3 = 700.0*scale;
H_canal = 30.0*scale;
H_puit  = 30.0*scale;
w_well  = 60.0*scale;

y_canal_bottom = H_interface_glass_pdms + 1*scale;
y_canal_top    = y_canal_bottom + H_canal;
y_well_top     = y_canal_top + H_puit;

x_w1 = -260.0*scale;
x_w2 = 0.0*scale;
x_w3 = 260.0*scale;

Point(13) = {-Lx3/2, y_canal_bottom, 0, lc_canal};
Point(14) = { Lx3/2, y_canal_bottom, 0, lc_canal};
Point(15) = { Lx3/2, y_canal_top,    0, lc_canal};

// Puits 3 (x = 260 µm)
Point(16) = {x_w3+w_well/2, y_canal_top, 0, lc_canal};
Point(17) = {x_w3+w_well/2, y_well_top,  0, lc_canal};
Point(18) = {x_w3-w_well/2, y_well_top,  0, lc_canal};
Point(19) = {x_w3-w_well/2, y_canal_top, 0, lc_canal};

// Puits 2 (x = 0 µm)
Point(20) = {x_w2+w_well/2, y_canal_top, 0, lc_canal};
Point(21) = {x_w2+w_well/2, y_well_top,  0, lc_canal};
Point(22) = {x_w2-w_well/2, y_well_top,  0, lc_canal};
Point(23) = {x_w2-w_well/2, y_canal_top, 0, lc_canal};

// Puits 1 (x = -260 µm)
Point(24) = {x_w1+w_well/2, y_canal_top, 0, lc_canal};
Point(25) = {x_w1+w_well/2, y_well_top,  0, lc_canal};
Point(26) = {x_w1-w_well/2, y_well_top,  0, lc_canal};
Point(27) = {x_w1-w_well/2, y_canal_top, 0, lc_canal};

Point(28) = {-Lx3/2, y_canal_top, 0, lc_canal};

Line(16) = {13, 14};
Line(17) = {14, 15};
Line(18) = {15, 16};
Line(19) = {16, 17};
Line(20) = {17, 18};
Line(21) = {18, 19};
Line(22) = {19, 20};
Line(23) = {20, 21};
Line(36) = {21, 22};
Line(37) = {22, 23};
Line(38) = {23, 24};
Line(39) = {24, 25};
Line(40) = {25, 26};
Line(41) = {26, 27};
Line(42) = {27, 28};
Line(43) = {28, 13};


// ============================================================
// 4. GOUTTES D'EAU
// ============================================================

p_rad = 20*scale;
p_cy  = 5.0*scale;

// Goutte 1 (x = -260 µm)
p1_cx = -260.0*scale;
Point(201) = {p1_cx, p_cy, 0, lc_part};
Point(202) = {p1_cx + p_rad, p_cy, 0, lc_part};
Point(203) = {p1_cx, p_cy + p_rad, 0, lc_part};
Point(204) = {p1_cx - p_rad, p_cy, 0, lc_part};
Point(205) = {p1_cx, p_cy - p_rad, 0, lc_part};

Circle(24) = {202, 201, 203};
Circle(25) = {203, 201, 204};
Circle(26) = {204, 201, 205};
Circle(27) = {205, 201, 202};

// Goutte 2 (x = 0 µm)
p2_cx = 0.0*scale;
Point(206) = {p2_cx, p_cy, 0, lc_part};
Point(207) = {p2_cx + p_rad, p_cy, 0, lc_part};
Point(208) = {p2_cx, p_cy + p_rad, 0, lc_part};
Point(209) = {p2_cx - p_rad, p_cy, 0, lc_part};
Point(210) = {p2_cx, p_cy - p_rad, 0, lc_part};

Circle(28) = {207, 206, 208};
Circle(29) = {208, 206, 209};
Circle(30) = {209, 206, 210};
Circle(31) = {210, 206, 207};

// Goutte 3 (x = +260 µm)
p3_cx = 260.0*scale;
Point(211) = {p3_cx, p_cy, 0, lc_part};
Point(212) = {p3_cx + p_rad, p_cy, 0, lc_part};
Point(213) = {p3_cx, p_cy + p_rad, 0, lc_part};
Point(214) = {p3_cx - p_rad, p_cy, 0, lc_part};
Point(215) = {p3_cx, p_cy - p_rad, 0, lc_part};

Circle(44) = {212, 211, 213};
Circle(45) = {213, 211, 214};
Circle(46) = {214, 211, 215};
Circle(47) = {215, 211, 212};


// ============================================================
// 5. SURFACES (DOMAINES)
// ============================================================

// PDMS (Domaine 1)
Curve Loop(1) = {34, 32, 4, 5};
Curve Loop(5) = {16:23, 36:43};
Plane Surface(1) = {1, 5};

// HUILE DU CANAL (Domaine 2)
Curve Loop(6) = {24:27};
Curve Loop(7) = {28:31};
Curve Loop(9) = {44:47};
Plane Surface(2) = {5, 6, 7, 9};

// EAU (Domaines 3, 4, 7)
Plane Surface(3) = {6};
Plane Surface(4) = {7};
Plane Surface(7) = {9};

// RESIN (Domaine 5)
Curve Loop(2) = {1, 2, -35, 6};
Curve Loop(3) = {8:11};
Curve Loop(4) = {12:15};
Plane Surface(5) = {2, 3, 4};

// GLASS (Domaine 6)
Curve Loop(8) = {50, 3, -34, 7};
Plane Surface(6) = {8};

// COUCHE MINCE HUILE FC-40 (Domaine 8)
Curve Loop(10) = {35, 48, -50, 49};
Plane Surface(8) = {10};


// ============================================================
// 6. GROUPES PHYSIQUES (MAPPING COURBES)
// ============================================================

Physical Curve("BEM_outside_top_E", 100) = {32, 4, 5};
Physical Curve("BEM_outside_bottom_E", 101) = {1, 2, 6};
Physical Curve("BEM_outside_glass_E", 108) = {3, 7};
Physical Curve("BEM_outside_fc40_E", 111)  = {48, 49};

// Interfaces jux corrigées selon l'orientation des lignes de 100->101, 104->105 et 102->103 :
// - Line(35) [100->101]: Domaine à gauche = 8 (FC40), à droite = 5 (RESIN) -> BEM_interface_8_5_jux
// - Line(50) [104->105]: Domaine à gauche = 6 (GLASS), à droite = 8 (FC40) -> BEM_interface_6_8_jux
// - Line(34) [102->103]: Domaine à gauche = 1 (PDMS), à droite = 6 (GLASS) -> BEM_interface_1_6_jux

Physical Curve("BEM_interface_8_5_jux", 102) = {35}; // FC40 / RESIN
Physical Curve("BEM_interface_6_8_jux", 112) = {50}; // GLASS / FC40
Physical Curve("BEM_interface_1_6_jux", 109) = {34}; // PDMS / GLASS

Physical Curve("BEM_interface_1_2_inc", 103) = {16:23, 36:43};
Physical Curve("BEM_interface_2_3_inc", 104) = {24:27};
Physical Curve("BEM_interface_2_4_inc", 105) = {28:31};
Physical Curve("BEM_interface_2_7_inc", 110) = {44:47};

Physical Curve("BEM_Hole_E1_H", 106) = {8:11};
Physical Curve("BEM_Hole_E2_H", 107) = {12:15};


// ============================================================
// 7. DOMAINES BEM
// ============================================================

Physical Surface("BEM_domain_1", 1) = {1}; // PDMS
Physical Surface("BEM_domain_2", 2) = {2}; // OIL CANAL
Physical Surface("BEM_domain_3", 3) = {3}; // WATER 1
Physical Surface("BEM_domain_4", 4) = {4}; // WATER 2
Physical Surface("BEM_domain_7", 7) = {7}; // WATER 3
Physical Surface("BEM_domain_5", 5) = {5}; // RESIN
Physical Surface("BEM_domain_6", 6) = {6}; // GLASS
Physical Surface("BEM_domain_8", 8) = {8}; // FC-40 LAYER


// ============================================================
// 8. PARAMÈTRES MATÉRIAUX
// ============================================================

SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL * eps0);
SetNumber("Materials/BEM_domain_3/Epsilon", epsilon_r_WATER * eps0);
SetNumber("Materials/BEM_domain_4/Epsilon", epsilon_r_WATER * eps0);
SetNumber("Materials/BEM_domain_7/Epsilon", epsilon_r_WATER * eps0);
SetNumber("Materials/BEM_domain_5/Epsilon", epsilon_r_RESIN * eps0);
SetNumber("Materials/BEM_domain_6/Epsilon", epsilon_r_GLASS * eps0);
SetNumber("Materials/BEM_domain_8/Epsilon", epsilon_r_FC40 * eps0);


// ============================================================
// 9. CONDITIONS AUX LIMITES
// ============================================================

SetNumber("Boundary Conditions/BEM_outside_top_E/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_outside_bottom_E/BEM_domain_5/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_outside_glass_E/BEM_domain_6/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_outside_fc40_E/BEM_domain_8/neumann", 0.0);

SetNumber("Boundary Conditions/BEM_Hole_E1_H/BEM_domain_5/dirichlet", 2000.0);
SetNumber("Boundary Conditions/BEM_Hole_E2_H/BEM_domain_5/dirichlet", 0.0);
