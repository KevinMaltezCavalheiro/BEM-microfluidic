scale = 1e-6;

bt = 0.2; // Beam tickness
bh = 1; // Beam height (position)

eL = 3; // Electrode base length
el = 1.5; // Electrode tip length
et = 1; // Electrode tickness
es = 1.5; // Electrode spacing

nbElectrodes = 10;

W = nbElectrodes*(es + eL) + 2*es;    // Width of the space
H = W;  // Height of the space

nEs = 12;
nEL = 15;
nBt = 8;
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

For i In {0:nbElectrodes-1}
    Point(6*i + 1) = {(es+eL)*i, bh + bt, 0, 1};
    Point(6*i + 2) = {(es+eL)*i + es, bh + bt, 0, 1};
    Point(6*i + 3) = {(es+eL)*i + es+e, bh + bt + et, 0, 1};
    Point(6*i + 4) = {(es+eL)*i + es+e+el, bh + bt + et, 0, 1};
    Point(6*i + 5) = {(es+eL)*i + es, bh, 0, 1};
    Point(6*i + 6) = {(es+eL)*i, bh, 0, 1};
EndFor
Point(6*nbElectrodes + 1) = {nbElectrodes*(es+eL), bh + bt, 0, 1};
Point(6*nbElectrodes + 2) = {nbElectrodes*(es+eL) + es, bh + bt, 0, 1};
Point(6*nbElectrodes + 3) = {nbElectrodes*(es+eL) + es, bh, 0, 1};
Point(6*nbElectrodes + 6) = {nbElectrodes*(es+eL), bh, 0, 1};

Point(6*nbElectrodes + 4) = {0, 0, 0, 1};
Point(6*nbElectrodes + 5) = {W, 0, 0, 1};
Point(6*nbElectrodes + 7) = {W, H, 0, 1};
Point(6*nbElectrodes + 8) = {0, H, 0, 1};

For i In {0:nbElectrodes-1}
    Point(6*nbElectrodes + 8 + i + 1) = {i*(es+eL) + es + e + el/2, bh + bt + et - ((el*e)/(2*et)), 0, 1};
EndFor

/*
 *  Lines
 */
Line(9*nbElectrodes + 1) = {6*nbElectrodes + 4, 6*nbElectrodes + 5};
Line(9*nbElectrodes + 2) = {6*nbElectrodes + 5, 6*nbElectrodes + 7};
Line(9*nbElectrodes + 3) = {6*nbElectrodes + 7, 6*nbElectrodes + 8};
Line(9*nbElectrodes + 4) = {6*nbElectrodes + 8, 1};
Line(9*nbElectrodes + 5) = {6, 6*nbElectrodes + 4};

Line(9*nbElectrodes + 9) = {6*nbElectrodes + 1, 6*nbElectrodes + 2};
Line(9*nbElectrodes + 6) = {6*nbElectrodes + 2, 6*nbElectrodes + 3};
Line(9*nbElectrodes + 7) = {6*nbElectrodes + 3, 6*nbElectrodes + 6};
Line(9*nbElectrodes + 8) = {6*nbElectrodes + 6, 6*nbElectrodes + 1};
For i In {0:nbElectrodes-1}
    Line(9*i + 1) = {6*i + 1, 6*i + 2};
    Line(9*i + 2) = {6*i + 2, 6*i + 3};
    // Line(9*i + 3) = {6*i + 3, 6*i + 4};
    Circle(9*i + 3) = {6*i + 3, 6*nbElectrodes + 8 + i + 1, 6*i + 4};
    Line(9*i + 4) = {6*i + 4, 6*(i+1) + 1};
    Line(9*i + 5) = {6*(i+1) + 1, 6*i + 2};
    Line(9*i + 6) = {6*i + 5, 6*i + 2};
    Line(9*i + 7) = {6*i + 5, 6*i + 6};
    Line(9*i + 8) = {6*i + 6, 6*i + 1};
    Line(9*i + 9) = {6*(i+1) + 6, 6*i + 5};
EndFor

/*
 *  Transfinite Curve
 */
For i In {0:nbElectrodes-1}
    Transfinite Curve{9*i + 1, 9*i + 7} = nEs + 1 Using Progression 1;
    Transfinite Curve{9*i + 2, 9*i + 4} = nEt + 1 Using Progression 1;
    Transfinite Curve{9*i + 3, 9*i + 5, 9*i + 9} = nEL + 1 Using Progression 1;
    Transfinite Curve{9*i + 8, 9*i + 6} = nBt + 1 Using Progression 1;
EndFor
Transfinite Curve{9*nbElectrodes + 9, 9*nbElectrodes + 7} = nEs + 1 Using Progression 1;
Transfinite Curve{9*nbElectrodes + 8, 9*nbElectrodes + 6} = nBt + 1 Using Progression 1;

Transfinite Curve{9*nbElectrodes + 1, 9*nbElectrodes + 3} = nbElectrodes*(nEs + nEL) + 2*nEs + 1 Using Progression 1;
Transfinite Curve{9*nbElectrodes + 2} = nbElectrodes*(nEs + nEL) + 2*nEs + 1 Using Progression 1;
Transfinite Curve{9*nbElectrodes + 4} = nbElectrodes*(nEs + nEL) + 2*nEs - nBt - nBh + 1 Using Progression 1;
Transfinite Curve{9*nbElectrodes + 5} = nBh + 1 Using Progression 1;

/*
 *  Curve Loop
 */

For i In {0:nbElectrodes-1}
    Curve Loop(3*i + 1) = {9*i + 1, -(9*i + 6), 9*i + 7, 9*i + 8};
    Curve Loop(3*i + 2) = {9*i + 3, 9*i + 4, 9*i + 5, 9*i + 2};
    Curve Loop(3*i + 3) = {-(9*i + 5), -(9*(i+1) + 8), 9*i + 9, 9*i + 6};
EndFor
Curve Loop(3*nbElectrodes + 1) = {9*nbElectrodes + 9, 9*nbElectrodes + 6, 9*nbElectrodes + 7, 9*nbElectrodes + 8};
// How to fill a vector in a loop ?

bemBoundary = {1:6*nbElectrodes + 8};
bemBoundaryWithoutE = {1:4*nbElectrodes + 8 + 2};
bemFemBoundary = {1:6*nbElectrodes + 3};
leftElectrodesFull = {1:nbElectrodes};
rightElectrodesFull = {1:nbElectrodes};
For i In {0:nbElectrodes-1}
    bemBoundary[4*i+1 - 1] = 9*i + 1;
    bemBoundary[4*i+2 - 1] = 9*i + 2;
    bemBoundary[4*i+3 - 1] = 9*i + 3;
    bemBoundary[4*i+4 - 1] = 9*i + 4;

    bemBoundaryWithoutE[2*i+1 - 1] = 9*i + 1;
    bemBoundaryWithoutE[2*i+2 - 1] = 9*i + 3;

    bemFemBoundary[4*i+1 - 1] = 9*i + 1;
    bemFemBoundary[4*i+2 - 1] = 9*i + 2;
    bemFemBoundary[4*i+3 - 1] = 9*i + 3;
    bemFemBoundary[4*i+4 - 1] = 9*i + 4;

    leftElectrodesFull[i] = 9*i + 2;
    rightElectrodesFull[i] = 9*i + 4;
EndFor
bemBoundary[4*nbElectrodes + 1 - 1] = 9*nbElectrodes + 9;
bemBoundary[4*nbElectrodes + 2 - 1] = 9*nbElectrodes + 6;
bemBoundary[4*nbElectrodes + 3 - 1] = 9*nbElectrodes + 7;

bemBoundaryWithoutE[2*nbElectrodes + 1 - 1] = 9*nbElectrodes + 9;
bemBoundaryWithoutE[2*nbElectrodes + 2 - 1] = 9*nbElectrodes + 6;
bemBoundaryWithoutE[2*nbElectrodes + 3 - 1] = 9*nbElectrodes + 7;

bemFemBoundary[4*nbElectrodes + 1 - 1] = 9*nbElectrodes + 9;
bemFemBoundary[4*nbElectrodes + 2 - 1] = 9*nbElectrodes + 6;
bemFemBoundary[4*nbElectrodes + 3 - 1] = 9*nbElectrodes + 7;
For i In {0:nbElectrodes-1}
    bemBoundary[4*nbElectrodes + 3 + 2*i+1 - 1] = 9*i + 9;
    bemBoundary[4*nbElectrodes + 3 + 2*i+2 - 1] = 9*i + 7;

    bemBoundaryWithoutE[2*nbElectrodes + 3 + 2*i+1 - 1] = 9*i + 9;
    bemBoundaryWithoutE[2*nbElectrodes + 3 + 2*i+2 - 1] = 9*i + 7;

    bemFemBoundary[4*nbElectrodes + 3 + 2*i+1 - 1] = 9*i + 9;
    bemFemBoundary[4*nbElectrodes + 3 + 2*i+2 - 1] = 9*i + 7;
EndFor
bemBoundary[6*nbElectrodes + 4 - 1] = 9*nbElectrodes + 5;
bemBoundary[6*nbElectrodes + 5 - 1] = 9*nbElectrodes + 1;
bemBoundary[6*nbElectrodes + 6 - 1] = 9*nbElectrodes + 2;
bemBoundary[6*nbElectrodes + 7 - 1] = 9*nbElectrodes + 3;
bemBoundary[6*nbElectrodes + 8 - 1] = 9*nbElectrodes + 4;

bemBoundaryWithoutE[4*nbElectrodes + 4 - 1] = 9*nbElectrodes + 5;
bemBoundaryWithoutE[4*nbElectrodes + 5 - 1] = 9*nbElectrodes + 1;
bemBoundaryWithoutE[4*nbElectrodes + 6 - 1] = 9*nbElectrodes + 2;
bemBoundaryWithoutE[4*nbElectrodes + 7 - 1] = 9*nbElectrodes + 3;
bemBoundaryWithoutE[4*nbElectrodes + 8 - 1] = 9*nbElectrodes + 4;

leftElectrodes = {1:nbElectrodes-1};
rightElectrodes = {1:nbElectrodes-1};
For i In {0:nbElectrodes-2}
    leftElectrodes[i] = leftElectrodesFull[i+1];
    rightElectrodes[i] = rightElectrodesFull[i];
EndFor

bemBoundaryWithoutE[4*nbElectrodes + 9 - 1] = 2;
bemBoundaryWithoutE[4*nbElectrodes + 10 - 1] = 9*nbElectrodes + 4;

Curve Loop(3*nbElectrodes + 2) = bemBoundary[];

/*
 *  Plane Surface
 */
For i In {0:nbElectrodes-1}
    Plane Surface(3*i + 1) = {3*i + 1};
    Plane Surface(3*i + 2) = {3*i + 2};
    Plane Surface(3*i + 3) = {3*i + 3};
EndFor
Plane Surface(3*nbElectrodes + 1) = {3*nbElectrodes + 1};
Plane Surface(3*nbElectrodes + 2) = {3*nbElectrodes + 2};

/*
 *  Surface operations
 */
Recombine Surface{1:3*nbElectrodes + 1};
Transfinite Surface{1:3*nbElectrodes + 1};

/*
 *  Physical Curve
 */
Physical Curve("clamp", 1) = {8};
Physical Curve("BEM_boundary_without_e", 2) = bemBoundaryWithoutE[];
Physical Curve("BEM_FEM_boundary", 3) = bemFemBoundary[];
Physical Curve("left_electrodes", 4) = leftElectrodes[];
Physical Curve("right_electrodes", 5) = rightElectrodes[];

Physical Surface("FEM_domain") = {1:3*nbElectrodes + 1};
Physical Surface("BEM_domain_1") = {3*nbElectrodes + 2};

/*
 *  Parameters
 */
 SetNumber("Boundary Conditions/clamp/ux", 0.0);
 SetNumber("Boundary Conditions/clamp/uy", 0.0);
 
//  SetNumber("Materials/FEM_domain/Young", 3e7);
SetNumber("Materials/FEM_domain/Young", 1e7);
 SetNumber("Materials/FEM_domain/Poisson", 0.3);
 SetNumber("Materials/FEM_domain/rho",7800);
 
 SetNumber("Volumic Forces/FEM_domain/bx", 0.0);
 SetNumber("Volumic Forces/FEM_domain/by", 0.0);
 
 SetNumber("Boundary Conditions/BEM_boundary_without_e/BEM_domain_1/neumann", 0);
 SetNumber("Boundary Conditions/left_electrodes/BEM_domain_1/dirichlet", 8);
 SetNumber("Boundary Conditions/right_electrodes/BEM_domain_1/dirichlet", 0);
 SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity
 
 SetNumber("Non_linear_solver", 1);