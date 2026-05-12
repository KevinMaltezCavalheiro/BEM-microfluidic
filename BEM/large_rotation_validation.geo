
h = 1;
H = 10;

n = 4;

Point(1) = {0, 0, 0, 0.5};
Point(2) = {0.9*H, 0, 0, 0.5};
Point(3) = {0.9*H, h, 0, 0.5};
Point(4) = {0, h, 0, 0.5};

Point(5) = {0.95*H, 0, 0, 0.5};
Point(6) = {0.95*H, h, 0, 0.5};

Point(7) = {0.95*H, H, 0, 0.5};
Point(8) = {0.9*H, H, 0, 0.5};

Point(9) = {H, 0, 0, 0.5};
Point(10) = {H, h, 0, 0.5};
Point(11) = {H, H, 0, 0.5};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Curve Loop(2) = {5, 6, 7, -2};
Plane Surface(2) = {2};

Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 3};
Curve Loop(3) = {-7, 8, 9, 10};
Plane Surface(3) = {3};

Line(11) = {5, 9};
Line(12) = {9, 10};
Line(13) = {10, 6};
Curve Loop(4) = {11, 12, 13, -6};
Plane Surface(4) = {4};

Line(14) = {10, 11};
Line(15) = {11, 7};
Curve Loop(5) = {14, 15, -8, -13};
Plane Surface(5) = {5};

Transfinite Curve {1, 3, 8, 10, 14} = 18*n+1 Using Progression 1;
Transfinite Curve {2, 4, 6, 12} = 2*n+1 Using Progression 1;
Transfinite Curve {5, 7, 9, 11, 13, 15} = n+1 Using Progression 1;
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};
Recombine Surface {5};

Physical Curve("left_edge", 1) = {4};
Physical Surface("FEM_domain", 2) = {1, 2, 3, 4, 5}; // the trick is to include both plane surfaces in one single domain
Physical Curve("top_edge", 3) = {9, 15};

Physical Point("point_A", 5) = {7};

F = 40000;

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); // ALWAYS NEED TO IMPOSE BOTH ux AND uy ON A GIVEN EDGE !! (pas très réaliste, faut y réfléchir)
SetNumber("Boundary Conditions/left_edge/uy", 0.);
SetNumber("Materials/FEM_domain/Young", 3e7);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800); //volumic mass of acier
SetNumber("Boundary Conditions/top_edge/tx", F); // ALWAYS NEED TO IMPOSE BOTH tx AND ty ON A GIVEN EDGE (realiste, OK) !
SetNumber("Boundary Conditions/top_edge/ty", 0); //set to some other value for vertical deflection
SetNumber("Volumic Forces/FEM_domain/bx",0.);
SetNumber("Volumic Forces/FEM_domain/by",0.); //set to -9.81 for gravity

SetNumber("Non_linear_solver",1); // activate non linear solver