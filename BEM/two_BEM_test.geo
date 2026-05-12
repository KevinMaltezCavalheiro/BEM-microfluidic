scale = 1;

h_tot = 10;
w_tot = 5;
h_FEM = 4;

nx = 80;
ny = 40;

Point(1) = {0, 0, 0, 2};
Point(2) = {w_tot, 0, 0, 2};
Point(3) = {w_tot, h_FEM, 0, 2};
Point(4) = {w_tot, h_tot-h_FEM, 0, 2};
Point(5) = {w_tot, h_tot, 0, 2};
Point(6) = {0, h_tot, 0, 2};
Point(7) = {0, h_tot-h_FEM, 0, 2};
Point(8) = {0, h_FEM, 0, 2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {8, 7};
Line(8) = {8, 1};
Line(9) = {3, 8};
Line(10) = {7, 4};


Curve Loop(1) = {1, 2, 9, 8};
Plane Surface(1) = {1};

Curve Loop(2) = {10, 4, 5, 6};
Plane Surface(2) = {2};

Curve Loop(3) = {9, 7, 10, 3};
Plane Surface(3) = {3};


Transfinite Curve {8, 2, 6, 4} = 2*ny;
Transfinite Curve {7, 3} = ny;
Transfinite Curve {1, 9, 10, 5} = nx;

Recombine Surface {1};

Mesh.ElementOrder = 1;

Physical Curve("left_FEM", 1) = {7};
Physical Curve("right_FEM", 2) = {3};

Physical Surface("FEM_domain",3) = {3};
Physical Surface("BEM_domain_1", 4) = {1};
Physical Surface("BEM_domain_2", 5) = {2};

SetNumber("Boundary Conditions/left_FEM/ux", 0.); 
SetNumber("Boundary Conditions/left_FEM/uy", 0.); 
SetNumber("Boundary Conditions/right_FEM/ux", 0);
SetNumber("Boundary Conditions/right_FEM/uy", 0);

SetNumber("Materials/FEM_domain/Young", 210e3);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800);

SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.);

Physical Curve("BEM_FEM_boundary_total", 6) = {9,10};
Physical Curve("BEM_FEM_boundary_1", 7) = {9};
Physical Curve("BEM_FEM_boundary_2", 8) = {10};
Physical Curve("Electrode1", 9) = {9};
Physical Curve("mass1", 10) = {1}; 

Physical Curve("Electrode2", 11) = {6};
Physical Curve("mass2", 12) = {4};
Physical Curve("rest_of_outside", 13) = {2, 3, 5, 7, 8};

phi1 = 100;
SetNumber("Boundary Conditions/mass1/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/Electrode1/BEM_domain_1/dirichlet", phi1);

SetNumber("Boundary Conditions/mass2/BEM_domain_2/dirichlet", 150);
SetNumber("Boundary Conditions/Electrode2/BEM_domain_2/dirichlet", 200);

Physical Curve("BEM_boundary_1", 14) = {8, 1, 2, 9};
Physical Curve("BEM_boundary_2", 15) = {10, 4, 5, 6};

Physical Curve("left1", 16) = {8};
Physical Curve("right1", 17) = {2};

Physical Curve("bottom2", 18) = {10};
Physical Curve("top2", 19) = {5};
SetNumber("Boundary Conditions/left1/BEM_domain_1/neumann", 0);
SetNumber("Boundary Conditions/right1/BEM_domain_1/neumann", 0);
SetNumber("Boundary Conditions/bottom2/BEM_domain_2/neumann", 0);
SetNumber("Boundary Conditions/top2/BEM_domain_2/neumann", 0);
