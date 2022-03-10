Lx = 5;
Ly = 2;
nx = 3;
ny = 2;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {4, 1, 2, 3};
//Plane Surface(1) = {1};
Transfinite Curve {3, 1} = nx+1 Using Progression 1;
Transfinite Curve {4, 2} = ny+1 Using Progression 1;
//Transfinite Surface {1};

//Recombine Surface {1}; // quads instead of triangles

//à gauche c'est le tag qu'on donne, il doit être différent pour chaque curve ou surface et à droite c'est l'élément définit correspondant 4 courbes et 1 surface
Physical Curve("bottom_edge", 1) = {1};
Physical Curve("right_edge", 2) = {2};
Physical Curve("top_edge", 3) = {3};
Physical Curve("left_edge", 4) = {4};
//Physical Curve("boundary", 5) = {1,2,3,4};
//Physical Surface("domain", 6) = {1};

// additional parameters given to the solver

// BC domain
//SetNumber("Materials/domain/Young", 210e3);
//SetNumber("Materials/domain/Poisson", 0.3);
//SetNumber("Materials/domain/rho",7800); //volumic mass of acier
//SetNumber("Volumic Forces/domain/bx",0.);
//SetNumber("Volumic Forces/domain/by",0.); //set to -9.81 for gravity

//BC top edge
SetNumber("Boundary Conditions/top_edge/tx", 0.); // ALWAYS NEED TO IMPOSE BOTH tx AND ty ON A GIVEN EDGE (realiste, OK) !
SetNumber("Boundary Conditions/top_edge/ty", 0.); //set to some non-zero value to induce vertical deflection

//BC right edge
SetNumber("Boundary Conditions/right_edge/tx", 21e3); // for simple tension conditions
SetNumber("Boundary Conditions/right_edge/ty", 0.);

//BC bottom edge
SetNumber("Boundary Conditions/bottom_edge/tx", 21e3);
SetNumber("Boundary Conditions/bottom_edge/ty", 0.);

// BC left edge
SetNumber("Boundary Conditions/left_edge/ux", 0.); // HERE YOU DO NOT HAVE TO IMPOSE BOTH ux and uy simultaneously ! (permet aussi de simuler appuis à roulettes)
//SetNumber("Boundary Conditions/left_edge/uy", 0.);
//+
Show "*";
