/* Electrode.geo - Correction Définitive du Chainage */

u = 1e-6;
lc = 1 * u;

eps0 = 8.8541878128e-12;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL = 9.6;

// ---- 1. POINTS ----
Point (6) = {95.3964 *u, 85.9858 *u, 0 *u, lc} ;
Point (5) = {95.3964 *u, -88.0821 *u, 0 *u, lc} ;
Point (7) = {7.77115 *u, 85.9858 *u, 0 *u, lc} ;
Point (8) = {5.68121 *u, 82.1795 *u, 0 *u, lc} ;
Point (9) = {3.76411 *u, 79.1136 *u, 0 *u, lc} ;
Point (10) = {1.50621 *u, 75.1961 *u, 0 *u, lc} ;
Point (11) = {-1.0958 *u, 69.9508 *u, 0 *u, lc} ;
Point (12) = {-3.95399 *u, 64.6636 *u, 0 *u, lc} ;
Point (13) = {-6.72755 *u, 60.0111 *u, 0 *u, lc} ;
Point (14) = {-8.72983 *u, 55.3271 *u, 0 *u, lc} ;
Point (15) = {-10.7353 *u, 51.0186 *u, 0 *u, lc} ;
Point (16) = {-13.9363 *u, 45.7295 *u, 0 *u, lc} ;
Point (17) = {-16.7481 *u, 42.323 *u, 0 *u, lc} ;
Point (18) = {-20.162 *u, 39.0828 *u, 0 *u, lc} ;
Point (19) = {-22.6756 *u, 36.9537 *u, 0 *u, lc} ;
Point (20) = {-24.9335 *u, 35.5911 *u, 0 *u, lc} ;
Point (21) = {-27.021 *u, 34.4414 *u, 0 *u, lc} ;
Point (22) = {-30.0331 *u, 33.102 *u, 0 *u, lc} ;
Point (23) = {-31.9561 *u, 32.4895 *u, 0 *u, lc} ;
Point (24) = {-33.2451 *u, 32.236 *u, 0 *u, lc} ;
Point (25) = {-34.7877 *u, 32.2571 *u, 0 *u, lc} ;
Point (26) = {-35.9499 *u, 32.4683 *u, 0 *u, lc} ;
Point (27) = {-37.0065 *u, 33.0809 *u, 0 *u, lc} ;
Point (28) = {-37.603 *u, 33.73 *u, 0 *u, lc} ;
Point (29) = {-38.1281 *u, 35.3573 *u, 0 *u, lc} ;
Point (30) = {-38.3574 *u, 37.4512 *u, 0 *u, lc} ;
Point (31) = {-38.0515 *u, 39.5287 *u, 0 *u, lc} ;
Point (32) = {-37.1974 *u, 44.6993 *u, 0 *u, lc} ;
Point (33) = {-36.3846 *u, 50.1832 *u, 0 *u, lc} ;
Point (34) = {-35.231 *u, 57.2363 *u, 0 *u, lc} ;
Point (35) = {-34.416 *u, 63.0757 *u, 0 *u, lc} ;
Point (36) = {-33.3999 *u, 68.7628 *u, 0 *u, lc} ;
Point (37) = {-32.4492 *u, 75.0725 *u, 0 *u, lc} ;
Point (38) = {-31.5009 *u, 80.3534 *u, 0 *u, lc} ;
Point (39) = {-30.5331 *u, 85.9858 *u, 0 *u, lc} ;
Point (40) = {-61.2474 *u, 85.9858 *u, 0 *u, lc} ;
Point (41) = {-60.6378 *u, 82.3057 *u, 0 *u, lc} ;
Point (42) = {-60.048 *u, 78.2774 *u, 0 *u, lc} ;
Point (43) = {-59.065 *u, 73.6596 *u, 0 *u, lc} ;
Point (44) = {-57.8855 *u, 67.961 *u, 0 *u, lc} ;
Point (45) = {-56.6071 *u, 59.9698 *u, 0 *u, lc} ;
Point (46) = {-55.6242 *u, 54.0747 *u, 0 *u, lc} ;
Point (47) = {-54.7395 *u, 47.6884 *u, 0 *u, lc} ;
Point (48) = {-53.6049 *u, 41.8405 *u, 0 *u, lc} ;
Point (49) = {-52.7516 *u, 37.2251 *u, 0 *u, lc} ;
Point (50) = {-52.2497 *u, 33.9642 *u, 0 *u, lc} ;
Point (51) = {-52.1995 *u, 31.7066 *u, 0 *u, lc} ;
Point (52) = {-52.7516 *u, 29.8504 *u, 0 *u, lc} ;
Point (53) = {-53.8056 *u, 28.2953 *u, 0 *u, lc} ;
Point (54) = {-55.1608 *u, 26.8906 *u, 0 *u, lc} ;
Point (55) = {-56.265 *u, 26.1381 *u, 0 *u, lc} ;
Point (56) = {-58.3228 *u, 25.9875 *u, 0 *u, lc} ;
Point (57) = {-60.983 *u, 25.9374 *u, 0 *u, lc} ;
Point (58) = {-63.7936 *u, 25.9374 *u, 0 *u, lc} ;
Point (59) = {-67.2569 *u, 26.0879 *u, 0 *u, lc} ;
Point (60) = {-70.4189 *u, 26.5394 *u, 0 *u, lc} ;
Point (61) = {-73.0288 *u, 27.1414 *u, 0 *u, lc} ;
Point (62) = {-76.492 *u, 28.7969 *u, 0 *u, lc} ;
Point (63) = {-79.1599 *u, 30.562 *u, 0 *u, lc} ;
Point (64) = {-81.4687 *u, 33.2209 *u, 0 *u, lc} ;
Point (65) = {-82.8741 *u, 35.4283 *u, 0 *u, lc} ;
Point (1) = {-84.7813 *u, 39.1407 *u, 0 *u, lc} ;
Point (2) = {-86.8894 *u, 43.5554 *u, 0 *u, lc} ;
Point (3) = {-89.0759 *u, 48.1207 *u, 0 *u, lc} ;
Point (4) = {-89.0759 *u, -88.0821 *u, 0 *u, lc} ;

// Interface y = 0
Point(100) = {-89.0759 * u, 0, 0, lc};
Point(101) = { 95.3964 * u, 0, 0, lc};

// ---- 2. LIGNES ----
// Electrode Droite (V+)
Line(1) = {7,8}; Line(2) = {8,9}; Line(3) = {9,10}; Line(4) = {10,11}; Line(5) = {11,12};
Line(6) = {12,13}; Line(7) = {13,14}; Line(8) = {14,15}; Line(9) = {15,16}; Line(10) = {16,17};
Line(11) = {17,18}; Line(12) = {18,19}; Line(13) = {19,20}; Line(14) = {20,21}; Line(15) = {21,22};
Line(16) = {22,23}; Line(17) = {23,24}; Line(18) = {24,25}; Line(19) = {25,26}; Line(20) = {26,27};
Line(21) = {27,28}; Line(22) = {28,29}; Line(23) = {29,30}; Line(24) = {30,31}; Line(25) = {31,32};
Line(26) = {32,33}; Line(27) = {33,34}; Line(28) = {34,35}; Line(29) = {35,36}; Line(30) = {36,37};
Line(31) = {37,38}; Line(32) = {38,39};

// Electrode Gauche (V-)
Line(33) = {40,41}; Line(34) = {41,42}; Line(35) = {42,43}; Line(36) = {43,44}; Line(37) = {44,45};
Line(38) = {45,46}; Line(39) = {46,47}; Line(40) = {47,48}; Line(41) = {48,49}; Line(42) = {49,50};
Line(43) = {50,51}; Line(44) = {51,52}; Line(45) = {52,53}; Line(46) = {53,54}; Line(47) = {54,55};
Line(48) = {55,56}; Line(49) = {56,57}; Line(50) = {57,58}; Line(51) = {58,59}; Line(52) = {59,60};
Line(53) = {60,61}; Line(54) = {61,62}; Line(55) = {62,63}; Line(56) = {63,64}; Line(57) = {64,65};
Line(58) = {65,1}; Line(59) = {1,2}; Line(60) = {2,3};

// Fermetures et Interface
Line(61) = {39, 40}; // Milieu Haut
Line(62) = {6, 7};   // Coin Haut Droit
Line(63) = {101, 6}; // Vertical Haut Droit
Line(64) = {100, 3}; // Vertical Haut Gauche (Sens: Bas -> Haut pour boucle anti-horaire)
Line(100) = {100, 101}; // INTERFACE

Line(65) = {4, 5};   // Bas horizontal
Line(66) = {100, 4}; // Vertical Bas Gauche
Line(67) = {5, 101}; // Vertical Bas Droit

// ---- 3. SURFACES ----
// PDMS (Haut)
// On part de 100 -> 63 -> 62 -> 1:32 -> 61 -> 33:60 -> -64
Curve Loop(1) = {100, 63, 62, 1:32, 61, 33:60, -64};
Plane Surface(1) = {1};

// OIL (Bas)
Curve Loop(2) = {-100, 66, 65, 67};
Plane Surface(2) = {2};

// ---- 4. PHYSIQUES ----
Physical Surface("BEM_domain_1", 1) = {1};
Physical Surface("BEM_domain_2", 2) = {2};

Physical Curve("positiveElectrode", 301) = {1:32};
Physical Curve("negativeElectrode", 302) = {33:60};
Physical Curve("BEM_interface_1_2_jux", 400) = {100};
Physical Curve("neutral_boundary", 303) = {61, 62, 63, 64, 65, 66, 67};

// Paramètres BEM
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * eps0);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL * eps0);
SetNumber("Boundary Conditions/positiveElectrode/BEM_domain_1/dirichlet", 2000.0);
SetNumber("Boundary Conditions/negativeElectrode/BEM_domain_1/dirichlet", 0.0);
SetNumber("Boundary Conditions/neutral_boundary/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/neutral_boundary/BEM_domain_2/neumann", 0.0);
