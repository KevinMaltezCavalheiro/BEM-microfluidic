Lx = 5;
Ly = 2;
nx = 50;
ny = 20;

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

Physical Curve("left_edge", 5) = {4};
Physical Surface("FEM_domain", 6) = {1};
Physical Curve("top_edge", 7) = {3};
Physical Curve("right_edge", 8) = {2};
Physical Point("fixed_node", 9) = {1};

Mesh.ElementOrder = 1;

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); 
//SetNumber("Boundary Conditions/left_edge/uy", 0.);
SetNumber("Materials/FEM_domain/Young", 210e3);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800); 
SetNumber("Boundary Conditions/top_edge/tx", 0.);
SetNumber("Boundary Conditions/top_edge/ty", 0.);
SetNumber("Boundary Conditions/right_edge/tx", 21e3); // for simple tension conditions
SetNumber("Boundary Conditions/right_edge/ty", 0.);
SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.);
SetNumber("Boundary Conditions/fixed_node/uy",2.);

SetNumber("Non_linear_solver",0);