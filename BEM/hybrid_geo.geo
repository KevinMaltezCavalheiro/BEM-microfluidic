scale = 4e-6;
Lx_poutre = 5*scale;
Ly_poutre = 2*scale;
h_poutre = 2*scale;
h_tot = 10*scale;
width = 10*scale;
nx = 20; // prend beaucoup de temps àpd de 200x40
ny = 8;

Point(1) = {0, h_poutre, 0, 2};
Point(2) = {Lx_poutre, h_poutre, 0, 2};
Point(3) = {Lx_poutre, h_poutre + Ly_poutre, 0, 2};
Point(4) = {0, h_poutre + Ly_poutre, 0, 2};

Point(5) = {0, 0, 0, 2};
Point(6) = {Lx_poutre, 0, 0, 2};
Point(7) = {width, 0, 0, 2};
Point(8) = {width, h_tot, 0, 2};
Point(9) = {0, h_tot, 0, 2};

Line(1) = {2, 1};
Line(2) = {3, 2};
Line(3) = {4, 3};
Line(4) = {1, 4};

Line(5) = {1, 5};
Line(6) = {9, 4};
Line(7) = {8, 9};
Line(8) = {7, 8};
Line(9) = {6, 7};
Line(10) = {5, 6};

Curve Loop(1) = {4, 3, 2, 1};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 10, 9, 8, 7, 6, 3, 2, 1};
Plane Surface(2) = {2};

Transfinite Curve {3, 1, 10, 9} = nx+1 Using Progression 1;
Transfinite Curve {7} = 2*nx+1 Using Progression 1;
Transfinite Curve {2, 4, 5} = ny+1 Using Progression 1;
Transfinite Surface {1};

Transfinite Curve {6} = 3*ny+1 Using Progression 1;
Transfinite Curve {8} = 5*ny+1 Using Progression 1;
//Transfinite Surface {2};

Recombine Surface {1}; // quads instead of triangles
//Recombine Surface {2};

Mesh.ElementOrder = 1;

Physical Curve("left_edge", 5) = {4};
Physical Surface("FEM_domain", 6) = {1};
Physical Curve("top_edge", 7) = {3};
Physical Curve("right_edge", 8) = {2};
Physical Point("fixed_node", 9) = {1};

Physical Surface("BEM_domain_1", 10) = {2};

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); // HERE YOU DO NOT HAVE TO IMPOSE BOTH ux and uy simultaneously ! (permet aussi de simuler appuis à roulettes)
SetNumber("Boundary Conditions/left_edge/uy", 0);
SetNumber("Materials/domain/Young", 210e3);
SetNumber("Materials/domain/Poisson", 0.3);
SetNumber("Materials/domain/rho",7800); //volumic mass of acier
SetNumber("Boundary Conditions/top_edge/tx", 0.); // ALWAYS NEED TO IMPOSE BOTH tx AND ty ON A GIVEN EDGE (realiste, OK) !
SetNumber("Boundary Conditions/top_edge/ty", 0.); //set to some non-zero value to induce vertical deflection
//SetNumber("Boundary Conditions/right_edge/tx", 21e3); // for simple tension conditions
SetNumber("Boundary Conditions/right_edge/ty", 0.);
SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.); //set to -9.81 for gravity
//SetNumber("Boundary Conditions/fixed_node/uy",2.);

// BEM geometry
Physical Curve("BEM_FEM_boundary", 11) = {1,2,3};
Physical Curve("mass", 12) = {10}; //lower electrode
Physical Curve("rest_of_outside", 13) = {5,6,7,8,9};
phi_top = 95; 
SetNumber("Boundary Conditions/mass/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/BEM_FEM_boundary/BEM_domain_1/dirichlet", phi_top);
SetNumber("Boundary Conditions/rest_of_outside/BEM_domain_1/neumann", 0);
//Physical Curve("BEM_boundary", 14) = {1,3,5,6,7,8,9,10,2};

Physical Curve("Electrode", 15) = {1};
//SetNumber("Boundary Conditions/bottom_edge/tx", 0);
//SetNumber("Boundary Conditions/bottom_edge/ty", -25000);