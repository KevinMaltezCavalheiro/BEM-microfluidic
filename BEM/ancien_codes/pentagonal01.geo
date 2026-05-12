
DefineConstant[
  N = {10, Name "Input/1Points "}
];

h = 1/N;

// Points 
Point(1) = {-1, 1, 0, h};
Point(2) = {-1, 0, 0, h};
Point(3) = {0, -1, 0, h};
Point(4) = {2, -1, 0, h};
Point(5) = {2, 1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line Loop(6) = {5, 1, 2, 3, 4};
Plane Surface(7) = {6};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Surface(1) = {7};
