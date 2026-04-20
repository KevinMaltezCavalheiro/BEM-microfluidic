// ============================================================
//  LaChoseDuStagiaire.geo  – 2-domain BEM example
//  Domain 1 : PDMS (outside rectangle with electrode holes)
//  Domain 2 : Oil  (inside the channel)
//  Interface : the channel boundary (shared between the two domains)
// ============================================================

scale = 1e-6;
epsilon_r_PDMS = 2.8;
epsilon_r_OIL  = 5.8;

// ---- Geometry parameters ----
Lx_out = 1000.0*scale;   Ly_out = 200.0*scale;
nx_out = 16;             ny_out = 4;
flux_out = 0.;

Lx1 = 100.0*scale;  Ly1 = 60.0*scale;
cx1 = 350.0*scale;  cy1 = 0.0*scale;
nx1 = 2;           ny1 = 2;
phi1 = 2000.;

Lx2 = 100.0*scale;  Ly2 = 60.0*scale;
cx2 = -350.0*scale; cy2 = 0.0*scale;
nx2 = 2;           ny2 = 2;
phi2 = 0.;

Lx3 = 300.0*scale;  Ly3 = 60.0*scale;
cx3 = 0.0*scale;    cy3 = 0.0*scale;
nx3 = 2;           ny3 = 2;

// ---- Points ----
// Outside rectangle
Point(1) = {-Lx_out/2, -Ly_out/2, 0, 1.0};
Point(2) = { Lx_out/2, -Ly_out/2, 0, 1.0};
Point(3) = { Lx_out/2,  Ly_out/2, 0, 1.0};
Point(4) = {-Lx_out/2,  Ly_out/2, 0, 1.0};

// Electrode 1
Point(5) = {cx1 - Lx1/2, cy1 - Ly1/2, 0, 1.0};
Point(6) = {cx1 - Lx1/2, cy1 + Ly1/2, 0, 1.0};
Point(7) = {cx1 + Lx1/2, cy1 + Ly1/2, 0, 1.0};
Point(8) = {cx1 + Lx1/2, cy1 - Ly1/2, 0, 1.0};

// Electrode 2
Point(9)  = {cx2 - Lx2/2, cy2 - Ly2/2, 0, 1.0};
Point(10) = {cx2 - Lx2/2, cy2 + Ly2/2, 0, 1.0};
Point(11) = {cx2 + Lx2/2, cy2 + Ly2/2, 0, 1.0};
Point(12) = {cx2 + Lx2/2, cy2 - Ly2/2, 0, 1.0};

// Channel (same points as before)
Point(13) = {cx3 - Lx3/2, cy3 - Ly3/2, 0, 1.0};
Point(14) = {cx3 - Lx3/2, cy3,          0, 1.0};
Point(15) = {cx3 - 40*scale, 0,          0, 1.0};
Point(16) = {cx3 - 40*scale, cy3+Ly3/2,  0, 1.0};
Point(17) = {cx3 + 40*scale, cy3+Ly3/2,  0, 1.0};
Point(18) = {cx3 + 40*scale, cy3,         0, 1.0};
Point(19) = {cx3 + Lx3/2,    cy3,         0, 1.0};
Point(20) = {cx3 + Lx3/2,    cy3-Ly3/2,   0, 1.0};

// ---- Lines ----
// Outside rectangle
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
Line(9)  = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// Channel boundary
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 13};

// ---- Transfinite mesh ----
Transfinite Curve {1} = nx_out + 1 Using Progression 1;
Transfinite Curve {2} = ny_out + 1 Using Progression 1;
Transfinite Curve {3} = nx_out + 1 Using Progression 1;
Transfinite Curve {4} = ny_out + 1 Using Progression 1;

Transfinite Curve {5} = ny1 + 1 Using Progression 1;
Transfinite Curve {6} = nx1 + 1 Using Progression 1;
Transfinite Curve {7} = ny1 + 1 Using Progression 1;
Transfinite Curve {8} = nx1 + 1 Using Progression 1;

Transfinite Curve {9}  = ny2 + 1 Using Progression 1;
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

// ---- Curve loops ----
Curve Loop(1) = {4, 1, 2, 3};        // outer rectangle
Curve Loop(2) = {8, 5, 6, 7};        // electrode 1 hole
Curve Loop(3) = {12, 9, 10, 11};     // electrode 2 hole
Curve Loop(4) = {13,14,15,16,17,18,19,20}; // channel

// ---- Surfaces ----
// Domain 1 (PDMS): outer rectangle with electrode holes AND channel hole.
Plane Surface(1) = {1, 2, 3, 4};

// Domain 2 (Oil): interior of the channel.
Plane Surface(2) = {4};

// ---- Physical groups ----

// Surfaces → BEM domains
Physical Surface("BEM_domain_1", 9)  = {1};   // PDMS
Physical Surface("BEM_domain_2", 10) = {2};   // Oil (channel)

// External boundary of domain 1
Physical Curve("outside", 5) = {1, 2, 3, 4};
Physical Curve("hole1",   6) = {5, 6, 7, 8};
Physical Curve("hole2",   7) = {9, 10, 11, 12};

// Interface: channel boundary is SHARED between domain 1 and domain 2.
// Convention: "BEM_interface_<domainA_name>_<domainB_name>"
// The name tokens must match the domain suffixes above (1 and 2 here).
Physical Curve("BEM_interface_1_2", 8) = {13,14,15,16,17,18,19,20};

// ---- Material properties ----
SetNumber("Materials/BEM_domain_1/Epsilon", epsilon_r_PDMS * 8.8541878128e-12);
SetNumber("Materials/BEM_domain_2/Epsilon", epsilon_r_OIL  * 8.8541878128e-12);

// ---- Boundary conditions (only on external boundaries) ----
// Domain 1 external BCs:
SetNumber("Boundary Conditions/outside/BEM_domain_1/neumann",   flux_out);
SetNumber("Boundary Conditions/hole1/BEM_domain_1/dirichlet",   phi1);
SetNumber("Boundary Conditions/hole2/BEM_domain_1/dirichlet",   phi2);

// Domain 2 has no external boundary → no BCs to set.
// (Its only boundary is the interface with domain 1, which is handled
//  automatically by the continuity conditions.)

// ---- Notes ----
// The channel curve loop (4) appears in BOTH Surface(1) (as a hole) and
// Surface(2) (as the outer boundary of the oil domain).  gmsh will generate
// the SAME nodes on both sides, which is exactly what the BEM subdomain
// method requires: a conformal mesh at the interface.
//
// If you add more interfaces in the future, just add:
//   Physical Curve("BEM_interface_<A>_<B>", <tag>) = {...};
// and set the corresponding Material/Epsilon values.  The C++ solver will
// detect them automatically.
