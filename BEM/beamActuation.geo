scale = 1e-6;

// THIS .geo FILE IS NOT USED IN THE REPORT

bt = 1; // Beam tickness
bh = 1; // Beam height (position)

eL = 3; // Electrode base length
el = 2.5; // Electrode tip length
et = 1; // Electrode tickness
es = 1.5; // Electrode spacing

W = 3*es + 3*eL; // Width of the space
H = bh + 3*(bt+et); // Height of the space

nEs = 9;
nEL = 12;
nBt = 6;
nEt = 3;
nBh = 3;

bt = bt * scale;
bh = bh * scale;
eL = eL * scale;
el = el * scale;
et = et * scale;
es = es * scale;
W = W * scale;
H = H * scale;

e = (eL - el)/2;    // Electrode constant

/*
 *  Points
 */

// Space
Point(1) = {0, 0, 0, 1};
Point(2) = {W, 0, 0, 1};
Point(3) = {W, H, 0, 1};
Point(4) = {0, H, 0, 1};

// Beam
Point(5) = {0, bh+bt, 0, 1};
Point(6) = {es, bh+bt, 0, 1};
Point(7) = {es+eL, bh+bt, 0, 1};
Point(8) = {2*es+eL, bh+bt, 0, 1};
Point(9) = {2*es+2*eL, bh+bt, 0, 1};
Point(10) = {3*es+2*eL, bh+bt, 0, 1};

Point(11) = {3*es+2*eL, bh, 0, 1};
Point(12) = {2*es+2*eL, bh, 0, 1};
Point(13) = {2*es+eL, bh, 0, 1};
Point(14) = {es+eL, bh, 0, 1};
Point(15) = {es, bh, 0, 1};
Point(16) = {0, bh, 0, 1};

// Electrodes
Point(17) = {es+e, bh+bt+et, 0, 1};
Point(18) = {es+e+el, bh+bt+et, 0, 1};
Point(19) = {2*es+eL+e, bh+bt+et, 0, 1};
Point(20) = {2*es+eL+e+el, bh+bt+et, 0, 1};

/*
 *  Lines
 */

 // Space
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {16, 1};

// Beam
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 5};
Line(18) = {15, 6};
Line(19) = {14, 7};
Line(20) = {13, 8};
Line(21) = {12, 9};

// Electrode
Line(22) = {6, 17};
Line(23) = {17, 18};
Line(24) = {18, 7};
Line(25) = {8, 19};
Line(26) = {19, 20};
Line(27) = {20, 9};

/*
 *  Transfinite Curve
 */
Transfinite Curve{6, 8, 10, 12, 14, 16} = nEs + 1 Using Progression 1;
Transfinite Curve{7, 9, 13, 15, 23, 26} = nEL + 1 Using Progression 1;
Transfinite Curve{17, 18, 19, 20, 21, 11} = nBt + 1 Using Progression 1;
Transfinite Curve{22, 24, 25, 27} = nEt + 1 Using Progression 1;

Transfinite Curve{1, 3} = 3*nEs + 3*nEL + 1 Using Progression 1;
Transfinite Curve(2) = 4 * (nBt + nEt) + 1 Using Progression 1;
Transfinite Curve(4) = 3*nEt + 2*nBt + 1 Using Progression 1;
Transfinite Curve(5) = nBh + 1 Using Progression 1;

/*
 *  Curve Loop
 */
Curve Loop(1) = {6, -18, 16, 17};
Curve Loop(2) = {7, -19, 15, 18};
Curve Loop(3) = {8, -20, 14, 19};
Curve Loop(4) = {9, -21, 13, 20};
Curve Loop(5) = {10, 11, 12, 21};
Curve Loop(6) = {23, 24, -7, 22};
Curve Loop(7) = {26, 27, -9, 25};
Curve Loop(8) = {1, 2, 3, 4, 6, 22, 23, 24, 8, 25, 26, 27, 10, 11, 12, 13, 14, 15, 16, 5};

/*
 *  Plane Surface
 */
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

/*
 *  Surface Manipulation
 */
Transfinite Surface{1, 2, 3, 4, 5, 6, 7};
Recombine Surface{1, 2, 3, 4, 5, 6, 7};

/*
 *  Physical Curve
 */
Physical Curve("clamp", 1) = {17};
Physical Curve("outside", 2) = {1, 2, 3, 4, 5};
Physical Curve("left_electrode", 3) = {24};
Physical Curve("right_electrode", 4) = {25};
Physical Curve("beam", 5) = {6, 22, 23, 8, 26, 10, 11, 12, 13, 14, 15, 16};
Physical Curve("BEM_FEM_boundary", 6) = {6, 22, 23, 24, 8, 26, 26, 27, 10, 11, 12, 13, 14, 15, 16};

/*
 *  Physical Surface
 */
Physical Surface("FEM_domain", 1) = {1, 2, 3, 4, 5, 6, 7};
Physical Surface("BEM_domain_1", 2) = {8};

/*
 *  Parameters
 */
SetNumber("Boundary Conditions/clamp/ux", 0.0);
SetNumber("Boundary Conditions/clamp/uy", 0.0);

SetNumber("Materials/FEM_domain/Young", 3e7);
SetNumber("Materials/FEM_domain/Poisson", 0.3);
SetNumber("Materials/FEM_domain/rho",7800);

SetNumber("Volumic Forces/FEM_domain/bx", 0.0);
SetNumber("Volumic Forces/FEM_domain/by", 0.0);

SetNumber("Boundary Conditions/outside/BEM_domain_1/neumann", 0);
SetNumber("Boundary Conditions/beam/BEM_domain_1/neumann", 0);
SetNumber("Boundary Conditions/left_electrode/BEM_domain_1/dirichlet", 200);
SetNumber("Boundary Conditions/right_electrode/BEM_domain_1/dirichlet", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity

SetNumber("Non_linear_solver", 1);