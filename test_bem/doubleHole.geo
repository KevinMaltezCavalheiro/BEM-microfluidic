// Outside
Lx_out = 2.0;
Ly_out = 1.5;
nx_out = 40;
ny_out = 30;
phi_out = 100.;

// Hole 1
Lx1 = 0.25;
Ly1 = 0.25;
cx1 = -0.5;
cy1 = 0.35;
nx1 = 8;
ny1 = 8;
phi1 = 215.;

// Hole 2
Lx2 = 0.25;
Ly2 = 0.25;
cx2 = 0.5;
cy2 = -0.35;
nx2 = 8;
ny2 = 8;
phi2 = 250.;

// Outside
Point(1) = {-Lx_out/2, -Ly_out/2, 0, 1.0};
Point(2) = {Lx_out/2, -Ly_out/2, 0, 1.0};
Point(3) = {Lx_out/2, Ly_out/2, 0, 1.0};
Point(4) = {-Lx_out/2, Ly_out/2, 0, 1.0};

// Hole 1
Point(5) = {cx1 - Lx1/2, cy1 - Ly1/2, 0, 1.0};
Point(6) = {cx1 - Lx1/2, cy1 + Ly1/2, 0, 1.0};
Point(7) = {cx1 + Lx1/2, cy1 + Ly1/2, 0, 1.0};
Point(8) = {cx1 + Lx1/2, cy1 - Ly1/2, 0, 1.0};

// Hole 2
Point(9) = {cx2 - Lx2/2, cy2 - Ly2/2, 0, 1.0};
Point(10) = {cx2 - Lx2/2, cy2 + Ly2/2, 0, 1.0};
Point(11) = {cx2 + Lx2/2, cy2 + Ly2/2, 0, 1.0};
Point(12) = {cx2 + Lx2/2, cy2 - Ly2/2, 0, 1.0};

// Outside
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Hole 1
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Hole 2
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {8, 5, 6, 7};
Curve Loop(3) = {12, 9, 10, 11};

Plane Surface(1) = {1, 2, 3};

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

//Recombine Surface {1};

Physical Curve("outside", 5) = {1, 2, 3, 4};
Physical Curve("hole1", 6) = {5, 6, 7, 8};
Physical Curve("hole2", 7) = {9, 10, 11, 12};

SetNumber("Boundary Conditions/outside/dirichlet", phi_out);
SetNumber("Boundary Conditions/hole1/dirichlet", phi1);
SetNumber("Boundary Conditions/hole2/dirichlet", phi2);