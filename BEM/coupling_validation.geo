scale = 1e-6;

Lx = 5;
Ly = 2;
l = 3;

n = 25;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {Lx*scale, 0, 0, 0.1};
Point(3) = {Lx*scale, Ly*scale, 0, 0.1};
Point(4) = {0, Ly*scale, 0, 0.2};

Point(5) = {Lx*scale + l*scale, 0, 0, 0.2};
Point(6) = {Lx*scale + l*scale, Ly*scale, 0, 0.2};

Line(1) = {2, 1};
Line(2) = {3, 2};
Line(3) = {4, 3};
Line(4) = {1, 4};
Curve Loop(1) = {4, 3, 2, 1};
Plane Surface(1) = {1};

Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Curve Loop(2) = {2, 5, 6, 7};
Plane Surface(2) = {2};


Transfinite Curve {3, 1} = n*Lx+1 Using Progression 1;
Transfinite Curve {2, 4, 6} = n*Ly+1 Using Progression 1;
Transfinite Curve {5, 7} = n*l+1 Using Progression 1;
Transfinite Surface {1};
Transfinite Surface {2};

Recombine Surface {1}; // quads instead of triangles
Recombine Surface {2};

Mesh.ElementOrder = 1;

Physical Curve("left_edge", 1) = {4};
Physical Surface("FEM_domain", 2) = {1};
Physical Point("fixed_node", 3) = {1};

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); 
//SetNumber("Boundary Conditions/left_edge/uy", 0.);
SetNumber("Boundary Conditions/fixed_node/uy", 0.);
SetNumber("Materials/FEM_domain/Young", 210e3);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800); 

Physical Curve("BEM_FEM_boundary", 4) = {2};
Physical Surface("BEM_domain_1", 5) = {2};
Physical Curve("right_BEM", 6) = {6};
Physical Curve("homogeneous_field", 7) = {5, 7};

phi = 90;
SetNumber("Boundary Conditions/right_BEM/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/BEM_FEM_boundary/BEM_domain_1/dirichlet", phi);
SetNumber("Boundary Conditions/homogeneous_field/BEM_domain_1/neumann", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity

SetNumber("Non_linear_solver", 0);