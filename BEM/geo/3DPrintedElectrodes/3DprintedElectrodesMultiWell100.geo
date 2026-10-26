// ============================================================
//  Electrodes_Channel_MultiWell_PDMS_Glass_Resin_v5.geo
//  Version optimisée pour BEM (Raffinement Glass & Inter-électrodes)
// ============================================================

scale = 1e-6;

// --- PROPRIÉTÉS DIÉLECTRIQUES ---
epsilon_r_PDMS  = 2.8;
epsilon_r_RESIN = 4.0;
epsilon_r_GLASS = 6.7;
epsilon_r_OIL   = 5.9; // 9.62 pour FC40 et 5.9 pour HFE7500
epsilon_r_WATER = 80.0;

eps0 = 8.8541878128e-12;


// ============================================================
// 1. GEOMÉTRIE DU DOMAINE DE SIMULATION
// ============================================================

Lout_x        = 3600*scale;  // Largeur totale : 3.6 mm (-1800 µm à +1800 µm)
Lout_y_top    = 1000*scale;  // Sommet PDMS : +1000 µm (~1 mm au-dessus des puits)
Lout_y_bottom = -1200*scale; // Profondeur Résine : -1200 µm (~900 µm sous les électrodes)

H_interface           = -120*scale; // Interface RESIN / GLASS
H_interface_glass_pdms = -35*scale;  // Interface GLASS / PDMS (verre = 85 µm)

// --- TAILLES DE MAILLAGE (BEM) ---
lc_ext   = 25*scale; // Bords extérieurs très grossiers
lc       = 10*scale; // Électrodes et maillage général
lc_glass = 5*scale;  // Maillage fin dédié à la couche de verre (Lines 3, 7, 34, 35)
lc_canal = 5*scale;  // Canal d'huile et puits


// ============================================================
// 2. CADRE EXTÉRIEUR ET COUCHE DE VERRE
// ============================================================

Point(1) = {-Lout_x/2, Lout_y_bottom, 0, lc_ext};
Point(2) = { Lout_x/2, Lout_y_bottom, 0, lc_ext};
Point(3) = { Lout_x/2, Lout_y_top,    0, lc_ext};
Point(4) = {-Lout_x/2, Lout_y_top,    0, lc_ext};

// Sommets de la couche de verre affectés de lc_glass pour raffiner les lignes 3, 7, 34, 35
Point(100) = {-Lout_x/2, H_interface, 0, lc_glass};
Point(101) = { Lout_x/2, H_interface, 0, lc_glass};

Point(102) = {-Lout_x/2, H_interface_glass_pdms, 0, lc_glass};
Point(103) = { Lout_x/2, H_interface_glass_pdms, 0, lc_glass};

// RESIN - Bas & côtés
Line(1) = {1, 2};
Line(2) = {2, 101};

// GLASS - Côté droit (Line 3)
Line(3) = {101, 103};

// PDMS - Côté droit
Line(32) = {103, 3};

// PDMS - Sommet
Line(4) = {3, 4};

// PDMS - Côté gauche
Line(5) = {4, 102};

// GLASS - Côté gauche (Line 7)
Line(7) = {102, 100};

// RESIN - Côté gauche
Line(6) = {100, 1};

// Interfaces horizontales
Line(35) = {100, 101}; // RESIN / GLASS
Line(34) = {102, 103}; // GLASS / PDMS


// ============================================================
// 3. ÉLECTRODES (Trous dans la Résine)
// ============================================================

w1 = 400*scale;            // Largeur d'électrode = 400 µm
h1 = 200*scale;            // Profondeur d'électrode = 200 µm

y_elec_top    = -121*scale;
y_elec_bottom = y_elec_top - h1;

// Électrode Droite (de +200 µm à +600 µm)
cx1 = 400*scale;

Point(5) = {cx1-w1/2, y_elec_bottom, 0, lc};
Point(6) = {cx1+w1/2, y_elec_bottom, 0, lc};
Point(7) = {cx1+w1/2, y_elec_top,    0, lc};
Point(8) = {cx1-w1/2, y_elec_top,    0, lc};

Line(8)  = {5, 6};
Line(9)  = {6, 7};
Line(10) = {7, 8};
Line(11) = {8, 5};

// Électrode Gauche (de -600 µm à -200 µm)
cx2 = -400*scale;

Point(9)  = {cx2-w1/2, y_elec_bottom, 0, lc};
Point(10) = {cx2+w1/2, y_elec_bottom, 0, lc};
Point(11) = {cx2+w1/2, y_elec_top,    0, lc};
Point(12) = {cx2-w1/2, y_elec_top,    0, lc};

Line(12) = {9, 10};
Line(13) = {10, 11};
Line(14) = {11, 12};
Line(15) = {12, 9};


// ============================================================
// 4. CANAL D'HUILE ET PUITS
// ============================================================

Lx3     = 1400.0*scale; // Longueur du canal = 1400 µm (-700 µm à +700 µm)
H_canal = 30.0*scale;   // Hauteur du canal
H_puit  = 30.0*scale;   // Hauteur supplémentaire du puit
w_well  = 60.0*scale;   // Largeur du puit

y_canal_bottom = H_interface_glass_pdms + 1*scale;
y_canal_top    = y_canal_bottom + H_canal;
y_well_top     = y_canal_top + H_puit;

x_w1 = -260.0*scale;
x_w2 = 0.0*scale;
x_w3 = 260.0*scale;

Point(13) = {-Lx3/2, y_canal_bottom, 0, lc_canal};
Point(14) = { Lx3/2, y_canal_bottom, 0, lc_canal};
Point(15) = { Lx3/2, y_canal_top,    0, lc_canal};

// Puits 3 (Droit, x = +260 µm)
Point(16) = {x_w3+w_well/2, y_canal_top, 0, lc_canal};
Point(17) = {x_w3+w_well/2, y_well_top,  0, lc_canal};
Point(18) = {x_w3-w_well/2, y_well_top,  0, lc_canal};
Point(19) = {x_w3-w_well/2, y_canal_top, 0, lc_canal};

// Puits 2 (Centre, x = 0 µm)
Point(20) = {x_w2+w_well/2, y_canal_top, 0, lc_canal};
Point(21) = {x_w2+w_well/2, y_well_top,  0, lc_canal};
Point(22) = {x_w2-w_well/2, y_well_top,  0, lc_canal};
Point(23) = {x_w2-w_well/2, y_canal_top, 0, lc_canal};

// Puits 1 (Gauche, x = -260 µm)
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
// 5. CHAMP DE RAFFINEMENT DANS LA ZONE INTER-ÉLECTRODES
// ============================================================

// Zone de fort gradient entre x = -200 µm et x = +200 µm
Field[1] = Box;
Field[1].VIn = 2*scale;            // Taille de maillage fine entre électrodes (2 µm)
Field[1].VOut = lc_ext;            // Taille de maillage par défaut hors zone
Field[1].XMin = -800*scale;        // Englobe la région inter-électrodes
Field[1].XMax = 800*scale;
Field[1].YMin = H_interface - 300*scale;
Field[1].YMax = H_interface_glass_pdms + 100*scale;
Field[1].Thickness = 100*scale;    // Transition douce du maillage

Background Field = 1;


// ============================================================
// 6. SURFACES
// ============================================================

// PDMS
Curve Loop(1) = {34, 32, 4, 5};
Curve Loop(5) = {16:23, 36:43};
Plane Surface(1) = {1, 5};

// HUILE (CANAL)
Plane Surface(2) = {5};

// RÉSINE
Curve Loop(2) = {1, 2, -35, 6};
Curve Loop(3) = {8:11};
Curve Loop(4) = {12:15};
Plane Surface(5) = {2, 3, 4};

// VERRE (GLASS)
Curve Loop(8) = {35, 3, -34, 7};
Plane Surface(6) = {8};


// ============================================================
// 7. GROUPES PHYSIQUES
// ============================================================

Physical Curve("BEM_outside_top_E", 100)    = {32, 4, 5};
Physical Curve("BEM_outside_bottom_E", 101) = {1, 2, 6};
Physical Curve("BEM_outside_glass_E", 108)  = {3, 7};

Physical Curve("BEM_interface_6_5_jux", 102) = {35};
Physical Curve("BEM_interface_1_6_jux", 109) = {34};

Physical Curve("BEM_interface_1_2_inc", 103) = {16:23, 36:43};

Physical Curve("BEM_Hole_E1_H", 106) = {8:11};
Physical Curve("BEM_Hole_E2_H", 107) = {12:15};


// ============================================================
// 8. DOMAINES BEM
// ============================================================

Physical Surface("BEM_domain_1", 1) = {1}; // PDMS
Physical Surface("BEM_domain_2", 2) = {2}; // OIL
Physical Surface("BEM_domain_5", 5) = {5}; // RESIN
Physical Surface("BEM_domain_6", 6) = {6}; // GLASS


// ============================================================
// 9. PARAMÈTRES MATÉRIAUX
// ============================================================

SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL * eps0);
SetNumber("Materials/BEM_domain_5/Epsilon", epsilon_r_RESIN * eps0);
SetNumber("Materials/BEM_domain_6/Epsilon", epsilon_r_GLASS * eps0);


// ============================================================
// 10. CONDITIONS AUX LIMITES
// ============================================================

SetNumber("Boundary Conditions/BEM_outside_top_E/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_outside_bottom_E/BEM_domain_5/neumann", 0.0);
SetNumber("Boundary Conditions/BEM_outside_glass_E/BEM_domain_6/neumann", 0.0);

SetNumber("Boundary Conditions/BEM_Hole_E1_H/BEM_domain_5/dirichlet", 1000.0);
SetNumber("Boundary Conditions/BEM_Hole_E2_H/BEM_domain_5/dirichlet", 0.0);
