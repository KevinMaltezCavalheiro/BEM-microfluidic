scale = 1e-5;

h_tot = 20*scale;
h_base = 3*scale;
h_pin = 2.5*scale;
h_space = 1*scale;

l_tot = 20*scale;
l_base = 1*scale;
l_pin = 0.8*scale;
l_space = 1.2*scale;

t_electrode = 1*scale; // width of fixed electrodes
// FEM domain
Point(1) = {-l_tot/2, -h_tot/2, 0, 0.2*scale};
Point(2) = {-l_tot/2, h_tot/2, 0, 0.2*scale};
Point(3) = {-(l_tot/2 - l_base), h_tot/2, 0, 0.2*scale};
Point(4) = {-(l_tot/2 - l_base), h_base/2, 0, 0.2*scale};
Point(5) = {-l_pin/2, h_base/2, 0, 0.2*scale};
Point(6) = {-l_pin/2, h_base/2 + h_pin, 0, 0.2*scale};
Point(7) = {l_pin/2, h_base/2 + h_pin, 0, 0.2*scale};
Point(8) = {l_pin/2, h_base/2, 0, 0.2*scale};
Point(9) = {(l_tot/2 - l_base), h_base/2, 0, 0.2*scale};
Point(10) = {(l_tot/2 - l_base), h_tot/2, 0, 0.2*scale};
Point(11) = {l_tot/2, h_tot/2, 0, 0.2*scale};
Point(12) = {l_tot/2, -h_tot/2, 0, 0.2*scale};
Point(13) = {(l_tot/2 - l_base), -h_tot/2, 0, 0.2*scale};
Point(14) = {(l_tot/2 - l_base), -h_base/2, 0, 0.2*scale};
Point(15) = {l_pin/2, -h_base/2, 0, 0.2*scale};
Point(16) = {l_pin/2, -(h_base/2 + h_pin), 0, 0.2*scale};
Point(17) = {-l_pin/2, -(h_base/2 + h_pin), 0, 0.2*scale};
Point(18) = {-l_pin/2, -h_base/2, 0, 0.2*scale};
Point(19) = {-(l_tot/2 - l_base), -h_base/2, 0, 0.2*scale};
Point(20) = {-(l_tot/2 - l_base), -h_tot/2, 0, 0.2*scale};

// FEM domain
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
Plane Surface(1) = {1};

Recombine Surface {1}; // quads instead of triangles

Physical Curve("bottom_left", 1) = {20};
Physical Curve("top_left", 2) = {2};
Physical Curve("top_right", 3) = {10};
Physical Curve("bottom_right", 4) = {12};

Physical Surface("FEM_domain", 5) = {1};
Physical Curve("BEM_FEM_boundary", 6) = {3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 19};

// mechanical properties and boundary conditions
SetNumber("Boundary Conditions/bottom_left/ux", 0.); // encastrement
SetNumber("Boundary Conditions/bottom_left/uy", 0);
SetNumber("Boundary Conditions/top_left/ux", 0.); // encastrement
SetNumber("Boundary Conditions/top_left/uy", 0);
SetNumber("Boundary Conditions/top_right/ux", 0.); // encastrement
SetNumber("Boundary Conditions/top_right/uy", 0);
SetNumber("Boundary Conditions/bottom_right/ux", 0.); // encastrement
SetNumber("Boundary Conditions/bottom_right/uy", 0);
SetNumber("Materials/domain/Young", 210e3); // A DETERMINER PRECISEMENT
SetNumber("Materials/domain/Poisson", 0.3);
SetNumber("Materials/domain/rho",7800); //volumic mass of acier
SetNumber("Volumic Forces/FEM_domain/bx",50); // acceleration of accelerometer
SetNumber("Volumic Forces/FEM_domain/by",0.); 

// BEM domain 1
Point(21) = {-(l_pin/2 + l_space + t_electrode), -(h_base/2 + h_pin + h_space + t_electrode), 0, 0.2*scale}; // first fixed electrode
Point(22) = {-(l_pin/2 + l_space + t_electrode), -(h_base/2 + h_space), 0, 0.2*scale};
Point(23) = {-(l_pin/2 + l_space), -(h_base/2 + h_space), 0, 0.2*scale};
Point(24) = {-(l_pin/2 + l_space), -(h_base/2 + h_pin + h_space), 0, 0.2*scale};
Point(25) = {l_pin/2 + l_space, -(h_base/2 + h_pin + h_space), 0, 0.2*scale};
Point(26) = {l_pin/2 + l_space, -(h_base/2 + h_space), 0, 0.2*scale};
Point(27) = {l_pin/2 + l_space + t_electrode, -(h_base/2 + h_space), 0, 0.2*scale};
Point(28) = {l_pin/2 + l_space + t_electrode, -(h_base/2 + h_pin + h_space + t_electrode), 0, 0.2*scale};

Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 28};
Line(28) = {28, 21};

Line(37) = {20, 13};

Curve Loop(2) = {37, 13, 14, 15, 16, 17, 18, 19};
Curve Loop(3) = {21, 22, 23, 24, 25, 26, 27, 28};
Plane Surface(2) = {2, 3};

Physical Surface("BEM_domain_1", 7) = {2};
Physical Curve("BEM_FEM_boundary_1", 8) = {13, 14, 15, 16, 17, 18, 19};
Physical Curve("electrode_1", 9) = {21, 22, 23, 24, 25, 26, 27, 28};
Physical Curve("outside_1", 10) = {37};

phi_1 = 30;
SetNumber("Boundary Conditions/BEM_FEM_boundary_1/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/electrode_1/BEM_domain_1/dirichlet", phi_1);
SetNumber("Boundary Conditions/outside_1/BEM_domain_1/neumann", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity

// BEM domain 2
Point(29) = {-(l_pin/2 + l_space + t_electrode), h_base/2 + h_space, 0, 0.2*scale}; // second fixed electrode
Point(30) = {-(l_pin/2 + l_space + t_electrode), h_base/2 + h_pin + h_space + t_electrode, 0, 0.2*scale}; 
Point(31) = {l_pin/2 + l_space + t_electrode, h_base/2 + h_pin + h_space + t_electrode, 0, 0.2*scale};
Point(32) = {l_pin/2 + l_space + t_electrode, h_base/2 + h_space, 0, 0.2*scale};
Point(33) = {l_pin/2 + l_space, h_base/2 + h_space, 0, 0.2*scale};
Point(34) = {l_pin/2 + l_space, h_base/2 + h_pin + h_space, 0, 0.2*scale};
Point(35) = {-(l_pin/2 + l_space), h_base/2 + h_pin + h_space, 0, 0.2*scale};
Point(36) = {-(l_pin/2 + l_space), h_base/2 + h_space, 0, 0.2*scale};


Line(29) = {29, 30};
Line(30) = {30, 31};
Line(31) = {31, 32};
Line(32) = {32, 33};
Line(33) = {33, 34};
Line(34) = {34, 35};
Line(35) = {35, 36};
Line(36) = {36, 29};

Line(38) = {10, 3};

Curve Loop(4) = {38, 3, 4, 5, 6, 7, 8, 9};
Curve Loop(5) = {29, 30, 31, 32, 33, 34, 35, 36};
Plane Surface(3) = {4, 5};

Physical Surface("BEM_domain_2", 11) = {3};
Physical Curve("BEM_FEM_boundary_2", 12) = {3, 4, 5, 6, 7, 8, 9};
Physical Curve("electrode_2", 13) = {29, 30, 31, 32, 33, 34, 35, 36};
Physical Curve("outside_2", 14) = {38};

phi_2 = 30;
SetNumber("Boundary Conditions/BEM_FEM_boundary_2/BEM_domain_2/dirichlet", 0);
SetNumber("Boundary Conditions/electrode_2/BEM_domain_2/dirichlet", phi_2);
SetNumber("Boundary Conditions/outside_2/BEM_domain_2/neumann", 0);
SetNumber("Materials/BEM_domain_2/Epsilon", 8.8541878128e-12); // dielectric permittivity