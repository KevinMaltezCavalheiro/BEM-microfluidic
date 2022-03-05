Lx = 10;
Ly = 10;
nx = 5;
ny = 5;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
Transfinite Curve {3, 1} = nx+1 Using Progression 1;
Transfinite Curve {4, 2} = ny+1 Using Progression 1;
Transfinite Surface {1};

Recombine Surface {1}; // quads instead of triangles

Physical Curve("bottom_edge", 1) = {1};
Physical Curve("left_edge", 2) = {2};
Physical Surface("domain", 1) = {5};
Physical Curve("top_edge",3) = {3};
Physical Curve("right_edge", 4) = {4};
Physical Point("top_left", 1) = {6};
Physical Point("top_right", 2) = {7};
Physical Point("bottom_right", 3) = {8};
Physical Point("bottom_left", 4) = {9};


// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.);
SetNumber("Boundary Conditions/left_edge/uy", 0.5);
SetNumber("Materials/domain/Young", 210e3);
SetNumber("Materials/domain/Poisson", 0.3);
SetNumber("Boundary Conditions/top_edge/tx", 0.);
SetNumber("Boundary Conditions/top_edge/ty", -100.);