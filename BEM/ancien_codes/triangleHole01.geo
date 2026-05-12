Include "options01_data.geo";
h = 1/N;
R=Sqrt(3.);
r=0.05;

// Center Points of circle
Point(1) = {0, R/3., 0, h}; // Point A
Point(2) = {-0.5, -R/6., 0, h};
Point(3) = {0.5,-R/6., 0, h};
Point(4) = {0,0, 0, h};
Point(5) = {r,0, 0, h};
Point(6) = {-r,0, 0, h};
Line(1) = {2, 3};
Line(2) = {3, 1};
Line(3) = {1, 2};
Circle(4) = {5, 4, 6};
Circle(5) = {6, 4, 5};
Line Loop(6) = {3, 1, 2};
Line Loop(7) = {4, 5};
Plane Surface(8) = {6, 7};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(10) = {4, 5};
Physical Surface(1) = {8};
