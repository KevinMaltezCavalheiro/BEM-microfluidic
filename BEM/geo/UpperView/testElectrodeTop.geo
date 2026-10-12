// Généré par dxf_to_geo_v4.py

lc = 1.0;

// ── Points ──────────────────────────────────────────────
Point(1) = {1288.785545, 830.7032745, 0, lc};
Point(2) = {1582.579369, 830.7032745, 0, lc};
Point(3) = {1582.269923, 849.6397623, 0, lc};
Point(4) = {1288.404173, 850.5653622, 0, lc};
Point(5) = {1289.867559, 774.3512613, 0, lc};
Point(6) = {1583.734087, 773.3851092, 0, lc};
Point(7) = {1583.425646, 789.4489101, 0, lc};
Point(8) = {1289.559237, 790.4088584, 0, lc};
Point(9) = {1358.917403, 812.615323, 0, lc};
Point(10) = {1347.214161, 812.615323, 0, lc};
Point(11) = {1353.065782, 812.615323, 0, lc};
Point(12) = {1322.041165, 812.5302062, 0, lc};
Point(13) = {1310.206524, 812.5302062, 0, lc};
Point(14) = {1316.123845, 812.5302062, 0, lc};
Point(15) = {1395.532423, 812.5550224, 0, lc};
Point(16) = {1383.769278, 812.5550224, 0, lc};
Point(17) = {1389.65085, 812.5550224, 0, lc};
Point(18) = {1431.635513, 812.4953896, 0, lc};
Point(19) = {1420.040193, 812.4953896, 0, lc};
Point(20) = {1425.837853, 812.4953896, 0, lc};
Point(21) = {1468.013474, 812.4447538, 0, lc};
Point(22) = {1456.153354, 812.4447538, 0, lc};
Point(23) = {1462.083414, 812.4447538, 0, lc};
Point(24) = {1504.47721, 812.1349371, 0, lc};
Point(25) = {1492.459452, 812.1349371, 0, lc};
Point(26) = {1498.468331, 812.1349371, 0, lc};
Point(27) = {1541.187942, 811.8658715, 0, lc};
Point(28) = {1529.047161, 811.8658715, 0, lc};
Point(29) = {1535.117551, 811.8658715, 0, lc};
Point(30) = {1578.627276, 811.6379789, 0, lc};
Point(31) = {1566.241139, 811.6379789, 0, lc};
Point(32) = {1572.434207, 811.6379789, 0, lc};

// ── Courbes ─────────────────────────────────────────────
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Circle(9) = {9, 11, 10};
Circle(10) = {10, 11, 9};
Circle(11) = {12, 14, 13};
Circle(12) = {13, 14, 12};
Circle(13) = {15, 17, 16};
Circle(14) = {16, 17, 15};
Circle(15) = {18, 20, 19};
Circle(16) = {19, 20, 18};
Circle(17) = {21, 23, 22};
Circle(18) = {22, 23, 21};
Circle(19) = {24, 26, 25};
Circle(20) = {25, 26, 24};
Circle(21) = {27, 29, 28};
Circle(22) = {28, 29, 27};
Circle(23) = {30, 32, 31};
Circle(24) = {31, 32, 30};
Line(25) = {1, 8};
Line(26) = {7, 2};

// ── Physical Curve (une par paire/domaine réelle) ──────
Physical Curve("BEM_interface_PDMS1_Oil_jux", 100) = {1};
Physical Curve("PDMS1_E", 101) = {2, 4};
Physical Curve("PDMS1_1000_E", 102) = {3};
Physical Curve("PDMS2_0_E", 103) = {5};
Physical Curve("PDMS2_E", 104) = {6, 8};
Physical Curve("BEM_interface_PDMS2_Oil_jux", 105) = {7};
Physical Curve("BEM_interface_Oil_Well2_inc", 106) = {9:10};
Physical Curve("BEM_interface_Oil_Well1_inc", 107) = {11:12};
Physical Curve("BEM_interface_Oil_Well3_inc", 108) = {13:14};
Physical Curve("BEM_interface_Oil_Well4_inc", 109) = {15:16};
Physical Curve("BEM_interface_Oil_Well5_inc", 110) = {17:18};
Physical Curve("BEM_interface_Oil_Well7_inc", 111) = {19:20};
Physical Curve("BEM_interface_Oil_Well8_inc", 112) = {21:22};
Physical Curve("BEM_interface_Oil_Well9_inc", 113) = {23:24};
Physical Curve("Oil_E", 114) = {25:26};

// ── Curve Loop ──────────────────────────────────────────
Curve Loop(1) = {25, -7, 26, -1};
Curve Loop(2) = {1:4};
Curve Loop(3) = {5:8};
Curve Loop(4) = {11:12};
Curve Loop(5) = {9:10};
Curve Loop(6) = {13:14};
Curve Loop(7) = {15:16};
Curve Loop(8) = {17:18};
Curve Loop(9) = {19:20};
Curve Loop(10) = {21:22};
Curve Loop(11) = {23:24};

// ── Plane Surface / Physical Surface ───────────────────
Plane Surface(1) = {1, 5, 4, 6, 7, 8, 9, 10, 11};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Physical Surface("BEM_domain_Oil", 1) = {1};
Physical Surface("BEM_domain_PDMS1", 2) = {2};
Physical Surface("BEM_domain_PDMS2", 3) = {3};
Physical Surface("BEM_domain_Well1", 4) = {4};
Physical Surface("BEM_domain_Well2", 5) = {5};
Physical Surface("BEM_domain_Well3", 6) = {6};
Physical Surface("BEM_domain_Well4", 7) = {7};
Physical Surface("BEM_domain_Well5", 8) = {8};
Physical Surface("BEM_domain_Well7", 9) = {9};
Physical Surface("BEM_domain_Well8", 10) = {10};
Physical Surface("BEM_domain_Well9", 11) = {11};

// ── Matériaux (HYPOTHÈSES DE REGROUPEMENT À VÉRIFIER) ──
eps0 = 8.8541878128e-12;
epsilon_r_OIL = 9.62;  // À COMPLÉTER — domaines : Oil
epsilon_r_PDMS = 2.8;  // À COMPLÉTER — domaines : PDMS1, PDMS2
epsilon_r_WELL = 80.0;  // À COMPLÉTER — domaines : Well1, Well2, Well3, Well4, Well5, Well7, Well8, Well9

SetNumber("Materials/BEM_domain_Oil/Epsilon", epsilon_r_OIL * eps0);
SetNumber("Materials/BEM_domain_PDMS1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_PDMS2/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_Well1/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well2/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well3/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well4/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well5/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well7/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well8/Epsilon", epsilon_r_WELL * eps0);
SetNumber("Materials/BEM_domain_Well9/Epsilon", epsilon_r_WELL * eps0);

// ── Conditions aux limites ──────────────────────────────
SetNumber("Boundary Conditions/PDMS1_E/BEM_domain_PDMS1/neumann", 0.0);
SetNumber("Boundary Conditions/PDMS1_1000_E/BEM_domain_PDMS1/dirichlet", 1000.0);
SetNumber("Boundary Conditions/PDMS2_0_E/BEM_domain_PDMS2/dirichlet", 0.0);
SetNumber("Boundary Conditions/PDMS2_E/BEM_domain_PDMS2/neumann", 0.0);
SetNumber("Boundary Conditions/Oil_E/BEM_domain_Oil/neumann", 0.0);

Mesh.Algorithm = 6;
Mesh 2;
