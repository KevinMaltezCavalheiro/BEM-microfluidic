//THIS CODE DEFINES THE GEOMETRY OF A CLAMPED BAR WITH TWO TYPES OF MESHES, TRIANGULAR AND QUADRANGULAR

Lx = 4;
Ly = 2;
nx = 20;
ny = 20;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {Lx, 0, 0, 0.2};
Point(3) = {Lx, Ly, 0, 0.2};
Point(4) = {0, Ly, 0, 0.1};
Point(5) = {Lx/2, 0, 0, 0.15};
Point(6) = {Lx/2, Ly, 0, 0.15};

Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};
Line(8) = {6, 5};

Curve Loop(1) = {1, 7, 5, 6};
Curve Loop(2) = {2, 3, 4, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Curve {1, 5} = nx+1 Using Progression 1;
Transfinite Curve {7, 6} = ny+1 Using Progression 1;
Transfinite Curve {2, 4} = nx+1 Using Progression 1;
Transfinite Curve {3, 8} = ny+1 Using Progression 1;
//Transfinite Surface {1};
//Transfinite Surface {2};

Mesh.ElementOrder = 1;

Recombine Surface {1}; // quads instead of triangles

Physical Curve("left_edge", 5) = {6};
Physical Curve("right_edge", 6) = {3};
Physical Curve("bottom_edge", 7) = {1, 2};
Physical Curve("top_edge", 8) = {4, 5};
Physical Surface("BEM_domain", 6) = {1, 2};
Physical Curve("between", 9) = {7, 8};

SetNumber("Boundary Conditions/left_edge/dirichlet", 200.);
SetNumber("Boundary Conditions/right_edge/dirichlet", 100.);
SetNumber("Boundary Conditions/bottom_edge/neumann", 0.);
SetNumber("Boundary Conditions/top_edge/neumann", 0.);
