scale = 4e-6;

h_tot = 5*scale;
w_tot = 10*scale;
h_foot = 4*scale;
h_diaph = 4.4*scale;
h_end_diaph = 4.6*scale;
Lx_diaph = 1*scale;

nx = 10; // prend beaucoup de temps àpd de 200x40
ny = 4;

Point(1) = {0, 0, 0, 2};
Point(2) = {w_tot, 0, 0, 2};
Point(3) = {w_tot, h_foot, 0, 2};
Point(4) = {w_tot, h_tot, 0, 2};
Point(5) = {w_tot-Lx_diaph, h_end_diaph, 0, 2};
Point(6) = {Lx_diaph, h_end_diaph, 0, 2};
Point(7) = {0, h_tot, 0, 2};
Point(8) = {0, h_foot, 0, 2};
Point(9) = {Lx_diaph, h_diaph, 0, 2};
Point(10) = {w_tot-Lx_diaph, h_diaph, 0, 2};

Line(1) = {9, 8};
Line(2) = {10, 9};
Line(3) = {3, 10};
Line(4) = {4, 3};
Line(5) = {5, 4};
Line(6) = {6, 5};
Line(7) = {7, 6};
Line(8) = {8, 7};

Line(9) = {8, 1};
Line(10) = {1, 2};
Line(11) = {2 ,3};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

Curve Loop(2) = {9, 10, 11, 3, 2, 1};
Plane Surface(2) = {2}; 

Transfinite Curve {10} = 10*nx+1 Using Progression 1;
Transfinite Curve {1, 2, 3} = 3*nx+1 Using Progression 1;
Transfinite Curve {7, 6, 5} = 3*nx+1 Using Progression 1;
Transfinite Curve {8, 4} = 2*ny+1 Using Progression 1;
Transfinite Curve {7, 1, 5, 3} = 1*ny+1 Using Progression 1;
Transfinite Curve {9, 11}  = 4*ny+1 Using Progression 1;

//Transfinite Surface {1};

Recombine Surface {1};

Mesh.ElementOrder = 1;

Physical Curve ("left_edge_foot", 5) = {8};
Physical Curve ("right_edge_foot", 6) = {4};
Physical Curve ("left_low_foot", 7) = {1};
Physical Curve ("left_up_foot", 8) = {7};
Physical Curve ("right_low_foot", 9) = {3};
Physical Curve ("right_up_foot", 10) = {5};
Physical Curve ("up_diaph", 11) = {6};
Physical Curve ("low_diaph", 12) = {2};

Physical Surface("FEM_domain", 13) = {1};
Physical Surface("BEM_domain", 14) = {2};


SetNumber("Boundary Conditions/left_edge_foot/ux", 0.); 
SetNumber("Boundary Conditions/left_edge_foot/uy", 0);
SetNumber("Boundary Conditions/right_edge_foot/ux", 0.);
SetNumber("Boundary Conditions/right_edge_foot/uy", 0);
SetNumber("Boundary Conditions/left_low_foot/ux", 0.);
SetNumber("Boundary Conditions/left_low_foot/uy", 0);
SetNumber("Boundary Conditions/left_up_foot/ux", 0.);
SetNumber("Boundary Conditions/left_up_foot/uy", 0);
SetNumber("Boundary Conditions/right_low_foot/ux", 0.);
SetNumber("Boundary Conditions/right_low_foot/uy", 0);
SetNumber("Boundary Conditions/right_up_foot/ux", 0.);
SetNumber("Boundary Conditions/right_up_foot/uy", 0);

SetNumber("Materials/FEM_domain/Young", 210e3);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800);

SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.);

//BEM Geometry
Physical Curve("BEM_FEM_boundary", 15) = {3,2,1};
Physical Curve("Electrode", 18) = {2};
Physical Curve("mass", 16) = {10}; 
Physical Curve("rest_of_outside", 17) = {9,11,1,3,8,4,7,5,6};

phi_top = 20;
SetNumber("Boundary Conditions/mass/dirichlet", 0);
SetNumber("Boundary Conditions/Electrode/dirichlet", phi_top);
SetNumber("Boundary Conditions/rest_of_outside/neumann", 0);

Physical Curve("BEM_boundary", 19) = {1,2,3,11,10,9};
