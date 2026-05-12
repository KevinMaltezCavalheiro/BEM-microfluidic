Lx_poutre = 5;
Ly_poutre = 2;
h_poutre = 2;
h_tot = 10;
width = 10;
nx = 10;
ny = 4;

// THIS .geo FILE IS NOT USED IN THE REPORT.

Point(1) = {0, h_poutre, 0, 2};
Point(2) = {Lx_poutre, h_poutre, 0, 2};
Point(3) = {Lx_poutre, h_poutre + Ly_poutre, 0, 2};
Point(4) = {0, h_poutre + Ly_poutre, 0, 2};

Point(5) = {0, 0, 0, 2};
Point(6) = {width, 0, 0, 2};
Point(7) = {width, h_tot, 0, 2};
Point(8) = {0, h_tot, 0, 2};

Point(9) = {Lx_poutre/2, h_poutre, 0, 2};
Point(10) = {Lx_poutre/2, h_poutre + Ly_poutre, 0, 2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 1};
Line(6) = {4, 8};
Line(7) = {8, 7};
Line(8) = {7, 6};
Line(9) = {6, 5};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 1, 2, 3, 6, 7, 8, 9};
Plane Surface(2) = {2};

Transfinite Curve {3, 1} = nx+1 Using Progression 1;
Transfinite Curve {2, 4} = ny+1 Using Progression 1;
Transfinite Surface {1};

Transfinite Curve {9,7} = 2*nx+1 Using Progression 1;
Transfinite Curve {8} = 5*ny+1 Using Progression 1;
Transfinite Curve {5} = ny+1 Using Progression 1;
Transfinite Curve {6} = 3*ny+1 Using Progression 1;
//Transfinite Surface {2};

Recombine Surface {1}; // quads instead of triangles
Recombine Surface {2};

Mesh.ElementOrder = 1;

Physical Curve("left_edge", 5) = {4};
Physical Surface("FEM_domain", 6) = {1};
Physical Curve("top_edge", 7) = {3};
Physical Curve("right_edge", 8) = {2};
Physical Point("fixed_node", 9) = {1};

Physical Surface("BEM_domain", 10) = {2};

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); 
//SetNumber("Boundary Conditions/left_edge/uy", 2.);
SetNumber("Materials/FEM_domain/Young", 210e3);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800); 
SetNumber("Boundary Conditions/top_edge/tx", 0.); 
SetNumber("Boundary Conditions/top_edge/ty", 0.); 
SetNumber("Boundary Conditions/right_edge/tx", 21e3);
SetNumber("Boundary Conditions/right_edge/ty", 0.);
SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.); //set to -9.81 for gravity
SetNumber("Boundary Conditions/fixed_node/uy",2.);
