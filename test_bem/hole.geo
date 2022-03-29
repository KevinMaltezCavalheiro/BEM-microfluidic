Lx_in = 1.0;
Ly_in = 0.1;
cx = 0.0;
cy = -0.6;

Lx_out = 2.0;
Ly_out = 1.5;

phi_in = 200.;
phi_out = 100.;

nx_in = 40;
ny_in = 4;

nx_out = 80;
ny_out = 60;

Point(1) = {-Lx_out/2, -Ly_out/2, 0, 1.0};
Point(2) = {Lx_out/2, -Ly_out/2, 0, 1.0};
Point(3) = {Lx_out/2, Ly_out/2, 0, 1.0};
Point(4) = {-Lx_out/2, Ly_out/2, 0, 1.0};


Point(5) = {-Lx_in/2 + cx, -Ly_in/2 + cy, 0, 1.0};
Point(6) = {-Lx_in/2 + cx, Ly_in/2 + cy, 0, 1.0};
Point(7) = {Lx_in/2 + cx, Ly_in/2 + cy, 0, 1.0};
Point(8) = {Lx_in/2 + cx, -Ly_in/2 + cy, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {8, 5, 6, 7};

Plane Surface(1) = {1, 2};

Transfinite Curve {1} = nx_out + 1 Using Progression 1;
Transfinite Curve {2} = ny_out + 1 Using Progression 1;
Transfinite Curve {3} = nx_out + 1 Using Progression 1;
Transfinite Curve {4} = ny_out + 1 Using Progression 1;

Transfinite Curve {5} = ny_in + 1 Using Progression 1;
Transfinite Curve {6} = nx_in + 1 Using Progression 1;
Transfinite Curve {7} = ny_in + 1 Using Progression 1;
Transfinite Curve {8} = nx_in + 1 Using Progression 1;

//Recombine Surface {1};

Physical Curve("inside", 5) = {5, 6, 7, 8};
Physical Curve("outside", 6) = {1, 2, 3, 4};

SetNumber("Boundary Conditions/inside/dirichlet", phi_in);
SetNumber("Boundary Conditions/outside/dirichlet", phi_out);