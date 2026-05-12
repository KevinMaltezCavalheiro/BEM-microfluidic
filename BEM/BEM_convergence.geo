l = 1.0;
h = 1.0;

phiLeft = 100.0;
phiRight = 200.0;

n = 20;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {l, 0, 0, 1.0};
Point(3) = {l, h, 0, 1.0};
Point(4) = {0, h, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
Transfinite Curve {3, 1} = n+1 Using Progression 1;
Transfinite Curve {4, 2} = n+1 Using Progression 1;
Transfinite Surface {1};

Recombine Surface {1}; // quads instead of triangles

Physical Curve("left_edge", 5) = {4};
Physical Curve("bottom_edge", 6) = {1};
Physical Curve("right_edge", 7) = {2};
Physical Curve("top_edge", 8) = {3};

Physical Curve("BEM_boundary", 9) = {1,2,3,4};

Physical Surface("BEM_domain_1", 10) = {1};

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/BEM_domain_1/dirichlet", phiLeft);
SetNumber("Boundary Conditions/bottom_edge/BEM_domain_1/neumann", 0.0);
SetNumber("Boundary Conditions/right_edge/BEM_domain_1/dirichlet", phiRight);
SetNumber("Boundary Conditions/top_edge/BEM_domain_1/neumann", 0.0);