N=20;
h = 1/N;
R=1.;
theta=Pi/2;

// Center Points of circle
Point(1) = {R*Cos(theta), R*Sin(theta), 0, h}; 
Point(2) = {R*Cos(theta+2*Pi/3),R*Sin(theta+2*Pi/3), 0, h};
Point(3) = {R*Cos(theta+4*Pi/3),R*Sin(theta+4*Pi/3), 0, h};
Line(1) = {2, 3};
Line(2) = {3, 1};
Line(3) = {1, 2};
Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};
Physical Surface(1) = {5};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};