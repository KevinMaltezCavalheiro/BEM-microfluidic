scale = 1e-6;
Lx_poutre = 25;
Ly_poutre = 1;
h_poutre = 1;
h_tot = 5;
width = 35;

n = 10; // FEM MESH DENSITY

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); 
SetNumber("Boundary Conditions/left_edge/uy", 0);
SetNumber("Materials/FEM_domain/Young", 150e9);
SetNumber("Materials/FEM_domain/Poisson", 0.27);
SetNumber("Materials/FEM_domain/rho",2300);
SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.);

phi_top = 90;
SetNumber("Boundary Conditions/mass/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/BEM_FEM_boundary/BEM_domain_1/dirichlet", phi_top);
SetNumber("Boundary Conditions/rest_of_outside/BEM_domain_1/neumann", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity

SetNumber("Non_linear_solver",1);

Point(1) = {0, h_poutre*scale, 0, 2};
Point(2) = {Lx_poutre*scale, h_poutre*scale, 0, 2};
Point(3) = {Lx_poutre*scale, h_poutre*scale + Ly_poutre*scale, 0, 2};
Point(4) = {0, h_poutre*scale + Ly_poutre*scale, 0, 2};

Point(5) = {0, 0, 0, 2};
Point(6) = {Lx_poutre*scale, 0, 0, 2};
Point(7) = {width*scale, 0, 0, 2};
Point(8) = {width*scale, h_tot*scale, 0, 2};
Point(9) = {0, h_tot*scale, 0, 2};

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

Transfinite Curve {3, 1, 10} = Lx_poutre*n+1 Using Progression 1;
Transfinite Curve {7} = width*n+1 Using Progression 1;
Transfinite Curve {8} = h_tot*n+1 Using Progression 1;
Transfinite Curve {2, 4} = Ly_poutre*n+1 Using Progression 1;
Transfinite Curve {5} = h_poutre*n+1 Using Progression 1;
Transfinite Curve {9} = (width-Lx_poutre)*n+1 Using Progression 1;
Transfinite Curve {6} = (h_tot-h_poutre-Ly_poutre)*n+1 Using Progression 1;
Transfinite Surface {1};

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

// BEM geometry
Physical Curve("BEM_FEM_boundary", 11) = {1,2,3};
Physical Curve("mass", 12) = {10}; //lower electrode
Physical Curve("rest_of_outside", 13) = {5,6,7,8,9};
Physical Curve("Electrode", 15) = {1};
//SetNumber("Boundary Conditions/bottom_edge/tx", 0);
//SetNumber("Boundary Conditions/bottom_edge/ty", -25000);
