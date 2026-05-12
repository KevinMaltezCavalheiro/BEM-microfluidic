//THIS CODE DEFINES THE GEOMETRY OF A CLAMPED BAR WITH TWO TYPES OF MESHES, TRIANGULAR AND QUADRANGULAR

Lx = 5;
Ly = 2;
nx = 50;
ny = 20;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {Lx, 0, 0, 0.2};
Point(3) = {Lx, Ly, 0, 0.2};
Point(4) = {0, Ly, 0, 0.1};
Point(5) = {Lx/2, 0, 0, 0.15}; 
Point(6) = {Lx/2, Ly, 0, 0.15}; 
Line(1) = {1, 5};
Line(2) = {2, 5};
Line(3) = {3, 2};
Line(4) = {6, 3};
Line(5) = {6, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};
Curve Loop(1) = {1, 7, 5, 6};
Curve Loop(2) = {7, 4, 3, 2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
//Transfinite Curve {1, 5} = nx+1 Using Progression 1;
//Transfinite Curve {7, 6} = ny+1 Using Progression 1;
//Transfinite Curve {2, 4} = nx+1 Using Progression 1;
//Transfinite Curve {3, 7} = ny+1 Using Progression 1;
//Transfinite Surface {1};
//Transfinite Surface {2};


Mesh.ElementOrder = 2;

Recombine Surface {1}; // quads instead of triangles

Physical Curve("left_edge", 5) = {6};
Physical Surface("domain", 6) = {1, 2}; // the trick is to include both plane surfaces in one single domain
Physical Curve("top_edge", 7) = {4, 5};
Physical Curve("right_edge", 8) = {3};
Physical Point("fixed_node", 9) = {1};

// additional parameters given to the solver
SetNumber("Boundary Conditions/left_edge/ux", 0.); // ALWAYS NEED TO IMPOSE BOTH ux AND uy ON A GIVEN EDGE !! (pas très réaliste, faut y réfléchir)
//SetNumber("Boundary Conditions/left_edge/uy", 0.);
SetNumber("Materials/domain/Young", 210e3);
SetNumber("Materials/domain/Poisson", 0.3);
SetNumber("Materials/domain/rho",7800); //volumic mass of acier
SetNumber("Boundary Conditions/top_edge/tx", 0.); // ALWAYS NEED TO IMPOSE BOTH tx AND ty ON A GIVEN EDGE (realiste, OK) !
SetNumber("Boundary Conditions/top_edge/ty", 0.); //set to some other value for vertical deflection
SetNumber("Boundary Conditions/right_edge/tx", 21e3); //for simple tension conditions
SetNumber("Boundary Conditions/right_edge/ty", 0.);
SetNumber("Volumic Forces/domain/bx",0.);
SetNumber("Volumic Forces/domain/by",0.); //set to -9.81 for gravity
SetNumber("Boundary Conditions/fixed_node/uy",0.);

Physical Curve("BEM_FEM_boundary", 10) = {1,2,3,4,5,6};
