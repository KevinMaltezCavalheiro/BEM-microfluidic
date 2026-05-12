
h = 1;
H = 2;

n = 10;

// THIS .geo FILE IS NOT USED IN THE REPORT.

Point(1) = {0, 0, 0, 0.5};
Point(2) = {(H-h), 0, 0, 0.5};
Point(3) = {H-h, h, 0, 0.5};
Point(4) = {0, h, 0, 0.5};

Point(5) = {H, 0, 0, 0.5};
Point(6) = {H, h, 0, 0.5};

Point(7) = {H, H, 0, 0.5};
Point(8) = {H-h, H, 0, 0.5};

Line(1) = {1, 2};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {2, 5};
Line(6) = {5, 6};

Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 3};

Curve Loop(1) = {1, 5, 6, 8, 9, 10, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1, 3, 8, 10} = (H-h)*n+1 Using Progression 1;
Transfinite Curve {4, 5, 6, 9} = n+1 Using Progression 1;

Recombine Surface {1};

Physical Curve("left_edge", 1) = {4};
Physical Surface("FEM_domain", 2) = {1}; // the trick is to include both plane surfaces in one single domain
Physical Curve("top_edge", 3) = {9};

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

Physical Curve("BEM_FEM_boundary", 4) = {4};//+

Mesh.Algorithm = 8;
