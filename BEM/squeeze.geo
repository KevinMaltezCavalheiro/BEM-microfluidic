scale = 4e-6;

h_tot0 = 5;
w_tot0 = 10;
h_foot0 = 4;
h_diaph0 = 4.4;
h_end_diaph0 = 4.6;
Lx_diaph0 = 1;

unit_l = 0.2; // serves as a reference

h_tot = h_tot0*scale;
w_tot = w_tot0*scale;
h_foot = h_foot0*scale;
h_diaph = h_diaph0*scale;
h_end_diaph = h_end_diaph0*scale;
Lx_diaph = Lx_diaph0*scale;

n = 8; // FEM elements density
nBEM = 0.25; // BEM elements density

// THIS .geo FILE IS NOT USED IN THE REPORT

Point(1) = {0, 0, 0, 0.8*nBEM*scale};
Point(2) = {w_tot, 0, 0, 0.8*nBEM*scale};
Point(3) = {w_tot, h_foot, 0, 0.4*nBEM*scale};
Point(4) = {w_tot, h_tot, 0, 2};
Point(5) = {w_tot-Lx_diaph, h_end_diaph, 0, 2};
Point(6) = {Lx_diaph, h_end_diaph, 0, 2};
Point(7) = {0, h_tot, 0, 2};
Point(8) = {0, h_foot, 0, 0.4*nBEM*scale};
Point(9) = {Lx_diaph, h_diaph, 0, 0.1*nBEM*scale};
Point(10) = {w_tot-Lx_diaph, h_diaph, 0, 0.1*nBEM*scale};

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

Line(13) = {9, 6};
Line(14) = {5, 10};

Curve Loop(1) = {-13, 1, 8, 7};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 13, 6, 14};
Plane Surface(2) = {2};

Curve Loop(3) = {-14, 5, 4, 3};
Plane Surface(3) = {3};

Curve Loop(4) = {9, 10, 11, 3, 2, 1};
Plane Surface(4) = {4}; 

Transfinite Curve {13, 14} = (h_end_diaph0 - h_diaph0)/unit_l*n+1 Using Progression 1;
Transfinite Curve {2, 6} = (w_tot0 - 2*Lx_diaph0)/unit_l*n+1 Using Progression 1;
Transfinite Curve {4, 8} = (h_tot0-h_foot0)/unit_l*n+1 Using Progression 1;
Transfinite Curve {1, 3, 5, 7} = (h_tot0-h_foot0)/unit_l*n+1 Using Progression 1;

Transfinite Surface {2};

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};

Mesh.ElementOrder = 1;

Physical Curve ("left_edge_foot", 5) = {8};
Physical Curve ("right_edge_foot", 6) = {4};
Physical Curve ("left_low_foot", 7) = {1};
Physical Curve ("left_up_foot", 8) = {7};
Physical Curve ("right_low_foot", 9) = {3};
Physical Curve ("right_up_foot", 10) = {5};
Physical Curve ("up_diaph", 11) = {6};
Physical Curve ("low_diaph", 12) = {2};

Physical Surface("FEM_domain", 13) = {1,2,3};
Physical Surface("BEM_domain_1", 14) = {4};


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

SetNumber("Materials/FEM_domain/Young", 150e9);
SetNumber("Materials/FEM_domain/Poisson", 0.27);
SetNumber("Materials/FEM_domain/rho",2300);

SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.);

//BEM Geometry
Physical Curve("BEM_FEM_boundary", 15) = {3,2,1};
Physical Curve("Electrode", 18) = {2};
Physical Curve("mass", 16) = {10}; 
Physical Curve("rest_of_outside", 17) = {9,11,1,3};

phi_top = 20;
SetNumber("Boundary Conditions/mass/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/Electrode/BEM_domain_1/dirichlet", phi_top);
SetNumber("Boundary Conditions/rest_of_outside/BEM_domain_1/neumann", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity
