w_inlet = 0.5;
w_outlet = 1.0;
L_inlet = 1.0;
L_outlet = 1.5;

phi_inlet = 200.;
phi_outlet = 100.;

n_trans_inlet = 10;
n_long_inlet = 10;
n_trans_outlet = 10;
n_long_outlet = 15;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {w_inlet, 0, 0, 1.0};
Point(3) = {w_inlet + L_outlet, 0, 0, 1.0};
Point(4) = {w_inlet + L_outlet, w_outlet, 0, 1.0};
Point(5) = {w_inlet, w_outlet, 0, 1.0};
Point(6) = {w_inlet, w_outlet + L_inlet, 0, 1.0};
Point(7) = {0, w_outlet + L_inlet, 0, 1.0};
Point(8) = {0, w_outlet, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Curve Loop(1) = {8, 1, 2, 3, 4, 5, 6, 7};

Plane Surface(1) = {1};

Transfinite Curve {1} = n_trans_inlet + 1 Using Progression 1;
Transfinite Curve {2} = n_long_outlet + 1 Using Progression 1;
Transfinite Curve {3} = n_trans_outlet + 1 Using Progression 1;
Transfinite Curve {4} = n_long_outlet + 1 Using Progression 1;
Transfinite Curve {5} = n_long_inlet + 1 Using Progression 1;
Transfinite Curve {6} = n_trans_inlet + 1 Using Progression 1;
Transfinite Curve {7} = n_long_inlet + 1 Using Progression 1;
Transfinite Curve {8} = n_trans_outlet + 1 Using Progression 1;

//Recombine Surface {1};

Physical Curve("sides", 5) = {1, 2, 4, 5, 7, 8};
Physical Curve("inlet", 6) = {3};
Physical Curve("outlet", 7) = {6};

SetNumber("Boundary Conditions/sides/neumann", 0.);
SetNumber("Boundary Conditions/inlet/dirichlet", phi_inlet);
SetNumber("Boundary Conditions/outlet/dirichlet", phi_outlet);