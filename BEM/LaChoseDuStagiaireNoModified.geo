scale = 1e-6;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL = 9.62;

// Outside
Lx_out = 1000.0*scale;
Ly_out = 200.0*scale;
nx_out = 160;
ny_out = 40;
flux_out = 0.;

// Electrode 1
Lx1 = 100.0*scale;
Ly1 = 60.0*scale;
cx1 = 350.0*scale;
cy1 = 0.0*scale;
nx1 = 12;
ny1 = 8;
phi1 = 2000.;

// Electrode 2
Lx2 = 100.0*scale;
Ly2 = 60.0*scale;
cx2 = -350.0*scale;
cy2 = 0.0*scale;
nx2 = 12;
ny2 = 8;
phi2 = 0.;

// Channel
Lx3 = 300.0*scale;
Ly3 = 60.0*scale;
cx3 = 0.0*scale;
cy3 = 0.0*scale;
nx3 = 12;
ny3 = 12;
flux3 = 1000.0;

// Outside
Point(1) = {-Lx_out/2, -Ly_out/2, 0, 1.0};
Point(2) = {Lx_out/2, -Ly_out/2, 0, 1.0};
Point(3) = {Lx_out/2, Ly_out/2, 0, 1.0};
Point(4) = {-Lx_out/2, Ly_out/2, 0, 1.0};

// Electrode 1
Point(5) = {cx1 - Lx1/2, cy1 - Ly1/2, 0, 1.0};
Point(6) = {cx1 - Lx1/2, cy1 + Ly1/2, 0, 1.0};
Point(7) = {cx1 + Lx1/2, cy1 + Ly1/2, 0, 1.0};
Point(8) = {cx1 + Lx1/2, cy1 - Ly1/2, 0, 1.0};

// Electrode 2
Point(9) = {cx2 - Lx2/2, cy2 - Ly2/2, 0, 1.0};
Point(10) = {cx2 - Lx2/2, cy2 + Ly2/2, 0, 1.0};
Point(11) = {cx2 + Lx2/2, cy2 + Ly2/2, 0, 1.0};
Point(12) = {cx2 + Lx2/2, cy2 - Ly2/2, 0, 1.0};

// Channel
Point(13) = {cx3 - Lx3/2, cy3 - Ly3/2, 0, 1.0};
Point(14) = {cx3 - Lx3/2, cy3, 0, 1.0};
Point(15) = {cx3 - 40*scale, 0, 0, 1.0};
Point(16) = {cx3 - 40*scale, cy3 + Ly3/2, 0, 1.0};
Point(17) = {cx3 + 40*scale, cy3 + Ly3/2, 0, 1.0};
Point(18) = {cx3 + 40*scale, cy3, 0, 1.0};
Point(19) = {cx3 + Lx3/2, cy3, 0, 1.0};
Point(20) = {cx3 + Lx3/2, cy3 - Ly3/2, 0, 1.0};

// Outside
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Electrode 1
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Electrode 2
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// Channel
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 13};

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {8, 5, 6, 7};
Curve Loop(3) = {12, 9, 10, 11};
Curve Loop(4) = {13, 14, 15, 16, 17, 18, 19, 20};

Transfinite Curve {1} = nx_out + 1 Using Progression 1;
Transfinite Curve {2} = ny_out + 1 Using Progression 1;
Transfinite Curve {3} = nx_out + 1 Using Progression 1;
Transfinite Curve {4} = ny_out + 1 Using Progression 1;

Transfinite Curve {5} = ny1 + 1 Using Progression 1;
Transfinite Curve {6} = nx1 + 1 Using Progression 1;
Transfinite Curve {7} = ny1 + 1 Using Progression 1;
Transfinite Curve {8} = nx1 + 1 Using Progression 1;

Transfinite Curve {9} = ny2 + 1 Using Progression 1;
Transfinite Curve {10} = nx2 + 1 Using Progression 1;
Transfinite Curve {11} = ny2 + 1 Using Progression 1;
Transfinite Curve {12} = nx2 + 1 Using Progression 1;

Transfinite Curve {13} = ny3 + 1 Using Progression 1;
Transfinite Curve {14} = nx3 + 1 Using Progression 1;
Transfinite Curve {15} = ny3 + 1 Using Progression 1;
Transfinite Curve {16} = nx3 + 1 Using Progression 1;
Transfinite Curve {17} = ny3 + 1 Using Progression 1;
Transfinite Curve {18} = nx3 + 1 Using Progression 1;
Transfinite Curve {19} = ny3 + 1 Using Progression 1;
Transfinite Curve {20} = 3*nx3 + 1 Using Progression 1;

//Recombine Surface {1};

// ===== AJOUTS/REMPLACEMENTS =====

// Matériaux
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS*8.8541878128e-12);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL*8.8541878128e-12);

// Surfaces de domaine
Plane Surface(1) = {1, 2, 3, 4};        // PDMS (comme avant)
Plane Surface(2) = {4};                  // NOUVEAU : huile (l’intérieur du channel)

Physical Surface("BEM_domain_1", 9) = {1};   // PDMS
Physical Surface("BEM_domain_2", 10) = {2};  // Huile

// Frontières
Physical Curve("outside", 5) = {1, 2, 3, 4};
Physical Curve("hole1", 6)    = {5, 6, 7, 8};
Physical Curve("hole2", 7)    = {9, 10, 11, 12};
Physical Curve("channel", 8)  = {13,14,15,16,17,18,19,20};

SetNumber("Boundary Conditions/outside/BEM_domain_1/neumann", flux_out);
SetNumber("Boundary Conditions/hole1/BEM_domain_1/dirichlet", phi1);
SetNumber("Boundary Conditions/hole2/BEM_domain_1/dirichlet", phi2);
SetNumber("Boundary Conditions/channel/BEM_domain_1/neumann", flux3);
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS*8.8541878128e-12); // dielectric permittivity
