scale = 4e-6;
h_tot = 4*scale;
w_tot = 10*scale;

h_mems = 4*scale;
h_diaph = 3.5*scale;
Lx_foot = 1*scale;
Lx_diaph = 2*scale;
nx = 20;
ny = 50;

// monter à 500+ steps en FEM non linéaire car converge très très lentement

Point(1) = {0, 0, 0, 2};
Point(2) = {Lx_foot, 0, 0, 2};
Point(3) = {Lx_diaph, h_diaph, 0, 2};
Point(4) = {w_tot - Lx_diaph, h_diaph, 0, 2};
Point(5) = {w_tot - Lx_foot, 0, 0, 2};
Point(6) = {w_tot, 0, 0, 2};
Point(7) = {w_tot, h_mems, 0, 2};
Point(8) = {0, h_mems, 0, 2};

Point(9) = {Lx_diaph, 0, 0, 2};
Point(10) = {w_tot - Lx_diaph, 0, 0, 2};
Point(11) = {w_tot, h_tot, 0, 2};
Point(12) = {0, h_tot, 0, 2};

Line(1) = {2, 1};
Line(2) = {3, 2};
Line(3) = {4, 3};
Line(4) = {5, 4};
Line(5) = {6, 5};
Line(6) = {7, 6};
Line(7) = {8, 7};
Line(8) = {1, 8};

Line(9) = {2, 9};
Line(10) = {9, 10};
Line(11) = {10, 5};
//Line(12) = {11, 7};
//Line(13) = {12, 11};
//Line(14) = {8, 12};

Curve Loop(1) = {8, 7, 6, 5, 4, 3, 2, 1};
Plane Surface(1) = {1};

//Curve Loop(2) = {-9, -10, -11, 4, 3, 2};
Curve Loop(2) = {2, 9, 10, 11, 4, 3};
Plane Surface(2) = {2};

Transfinite Curve{7} = 10*nx+1 Using Progression 1;
Transfinite Curve{1, 5} = nx+1 Using Progression 1;
Transfinite Curve{8, 6} = ny+1 Using Progression 1;
Transfinite Curve{10, 3} = 6*nx+1 Using Progression 1;
Transfinite Curve{2, 4} = 2*ny+1 Using Progression 1;
Transfinite Curve{9, 11} = nx+1 Using Progression 1;

Recombine Surface {1};

Mesh.ElementOrder = 1;

Physical Curve("diaphragme", 5) = {3};
Physical Surface("FEM_domain", 6) = {1};
Physical Curve("left_edge", 7) = {8};
Physical Curve("top_edge", 8) = {7};
Physical Curve("right_edge", 9) = {6};
Physical Curve("right_foot", 10) = {5};
Physical Curve("left_foot", 11) = {1};
Physical Curve("left_height", 12) = {2};
Physical Curve("right_height", 13) = {4};
Physical Point("fixed_node", 19) = {8};

Physical Surface("BEM_domain", 14) = {2};

SetNumber("Boundary Conditions/left_edge/ux", 0.); 
SetNumber("Boundary Conditions/left_edge/uy", 0);
SetNumber("Boundary Conditions/left_foot/ux", 0.); 
SetNumber("Boundary Conditions/left_foot/uy", 0);
SetNumber("Boundary Conditions/right_edge/ux", 0.); 
SetNumber("Boundary Conditions/right_edge/uy", 0);
SetNumber("Boundary Conditions/right_foot/ux", 0.); 
SetNumber("Boundary Conditions/right_foot/uy", 0);
SetNumber("Boundary Conditions/fixed_node/ux",0);
SetNumber("Boundary Conditions/fixed_node/uy",0);
SetNumber("Boundary Conditions/right_height/ux", 0.); 
SetNumber("Boundary Conditions/right_height/uy", 0);
SetNumber("Boundary Conditions/left_height/ux", 0.); 
SetNumber("Boundary Conditions/left_height/uy", 0);
SetNumber("Materials/FEM_domain/Young", 210e3);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800);

SetNumber("Boundary Conditions/top_edge/tx", 0.); // ALWAYS NEED TO IMPOSE BOTH tx AND ty ON A GIVEN EDGE (realiste, OK) !
SetNumber("Boundary Conditions/top_edge/ty", 0.);
SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.);

//BEM Geometry
Physical Curve("BEM_FEM_boundary", 15) = {2, 3, 4};
Physical Curve("mass", 16) = {10};//Lower electrode
Physical Curve("rest_of_outside", 17) = {9, 11, 2, 4};
Physical Curve("dirichlet", 20) = {3};
phi_top = 100;
SetNumber("Boundary Conditions/mass/dirichlet", 0);
SetNumber("Boundary Conditions/dirichlet/dirichlet", phi_top);
SetNumber("Boundary Conditions/rest_of_outside/neumann", 0);

Physical Curve("BEM_boundary", 18) = {9, 10, 11, 4, 3, 2};