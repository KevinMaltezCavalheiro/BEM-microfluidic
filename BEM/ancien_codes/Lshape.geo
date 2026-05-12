DefineConstant[
  N = {10, Name "Input/1Points "}
];
h = 1/N;
Point(1) = {0, 0, 0, h}; 
Point(2) = {2, 0, 0, h};
Point(3) = {2, 1, 0, h};
Point(4) = {1,1, 0, h};
Point(5) = {1,2, 0, h};
Point(6) = {0,2, 0, h};
Line(1) = {6, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line Loop(7) = {6, 1, 2, 3, 4, 5};
Plane Surface(8) = {7};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Surface(1) = {8};
