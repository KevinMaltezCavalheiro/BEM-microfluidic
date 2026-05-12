Include "options01_data.geo";
h = 1/N;

// Points 
Point(1) = {1, 0, 0, h};
Point(2) = {1, 1, 0, h};
Point(3) = {-1, 1, 0, h};
Point(4) = {-1, 0, 0, h};
Point(5) = {0, 0, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line Loop(11) = {1, 2, 3, 4, 5};
Plane Surface(14) = {11};
Physical Line(1) = {4};
Physical Line(2) = {5};
Physical Line(3) = {1};
Physical Line(4) = {2};
Physical Line(5) = {3};
Physical Surface(1) = {14};
// Mesh.Algorithm=6;
