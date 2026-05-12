scale = 1e-6;

nFEM = 2; // FEM element density

nBEM = 5/nFEM; // BEM element density (mettre ça dans quatrième coordonnée)

phi = 110;
SetNumber("Boundary Conditions/BEM_FEM_boundary/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/top_electrode/BEM_domain_1/dirichlet", phi);
SetNumber("Boundary Conditions/rest_of_outside/BEM_domain_1/neumann", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity

// mechanical properties and boundary conditions
SetNumber("Boundary Conditions/anchor/ux", 0.); // encastrement
SetNumber("Boundary Conditions/anchor/uy", 0);
SetNumber("Boundary Conditions/left_edge/ux", 0);
SetNumber("Materials/FEM_domain/Young", 160e9); // A DETERMINER PRECISEMENT
SetNumber("Materials/FEM_domain/Poisson", 0.22);
SetNumber("Materials/FEM_domain/rho",2300);
SetNumber("Volumic Forces/FEM_domain/bx",0);
SetNumber("Volumic Forces/FEM_domain/by",0);
SetNumber("Non_linear_solver",1);

Wc = 4;
g = 8;
Lg = 20;
Ls = 280;
Ht = 40;
Wt = 16;
Ws = 2;
Lc = 30;
y0 = 20;
Hbottom = 50;
Htop = 180;
Wtot = 400;

// scaled quantities
Wcs = Wc*scale;
gs = g*scale;
Lgs = Lg*scale;
Lss = Ls*scale;
Hts = Ht*scale;
Wts = Wt*scale;
Wss = Ws*scale;
Lcs = Lc*scale;
y0s = y0*scale;
Hbottoms = Hbottom*scale;
Htops = Htop*scale;
Wtots = Wtot*scale;

Ncombs = 4; // ne pas changer, pas généralisé

// fixed bottom anchor
y = -(Lg-2-4)*scale;
x = (5*g + 4.5*Wc)*scale;
Point(1) = {x, y-Wss-2*scale, 0, 10*scale}; // limite avec suite
Point(2) = {x, y-Wts, 0, 10*scale};
Point(3) = {x-Wts, y-Wts, 0, 10*scale};
Point(4) = {x-Wts, y, 0, 10*scale};
Point(5) = {x, y, 0, 10*scale};
Point(6) = {x, y-2*scale, 0, 10*scale}; // limite avec suite
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1}; // limite avec suite
Transfinite Curve{1} = (Wt - 2 - Ws)*nFEM+1 Using Progression 1;
Transfinite Curve{2, 3, 4} = Wt*nFEM+1 Using Progression 1;
Transfinite Curve{5} = 2*nFEM+1 Using Progression 1;
Transfinite Curve{6} = Ws*nFEM+1 Using Progression 1;
Curve Loop(1) = {1:6};
Plane Surface(1) = {1};
Transfinite Surface{1} = {2:5};
Recombine Surface{1};

// bottom long beam
Point(7) = {x+Lss, y-2*scale, 0, 10*scale}; // limite avec suite
Point(8) = {x+Lss, y-Wss-2*scale, 0, 10*scale}; // limite avec suite
Line(7) = {6, 7};
Line(8) = {7, 8}; // limite avec suite
Line(9) = {8, 1};
Transfinite Curve{7, 9} = Ls*nFEM+1 Using Progression 1;
Transfinite Curve{8} = Ws*nFEM+1 Using Progression 1;
Curve Loop(2) = {-6, 7:9};
Plane Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2};

//truss
x = x + Lss;
y = y - 10*scale - Wss;
Point(9) = {x, y + 8*scale + Wss + Lgs, 0 , 10*scale}; // limite avec suite
Point(10) = {x, y + 8*scale + 2*Wss + Lgs, 0 , 10*scale}; // limite avec suite
Point(11) = {x, y + 16*scale + 2*Wss + Lgs, 0 , 10*scale};
Point(12) = {x + Wts, y + 16*scale + 2*Wss + Lgs, 0 , 10*scale};
Point(13) = {x + Wts, y, 0 , 10*scale};
Point(14) = {x, y, 0 , 10*scale};
Line(10) = {7, 9};
Line(11) = {9, 10}; // limite avec suite
Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 13};
Line(15) = {13, 14};
Line(16) = {14, 8};
Transfinite Curve{10} = Lg*nFEM+1 Using Progression 1;
Transfinite Curve{11} = Ws*nFEM+1 Using Progression 1;
Transfinite Curve{12, 16} = 8*nFEM+1 Using Progression 1;
Transfinite Curve{13, 15} = Wt*nFEM+1 Using Progression 1;
Transfinite Curve{14} = (16 + 2*Ws + Lg)*nFEM+1 Using Progression 1;
Curve Loop(3) = {-8, 10:16};
Plane Surface(3) = {3};
Transfinite Surface{3} = {11:14};
Recombine Surface{3};

//top long beam
x = (5*g + 4.5*Wc)*scale;
Point(15) = {x, 4*scale, 0 , 10*scale}; // limite avec suite
Point(16) = {x, 4*scale + Wss, 0 , 10*scale}; // limite avec suite
Line(17) = {9, 15};
Line(18) = {15, 16}; // limite avec suite
Line(19) = {16, 10};
Transfinite Curve{17, 19} = Ls*nFEM+1 Using Progression 1;
Transfinite Curve{18} = Ws*nFEM+1 Using Progression 1;
Curve Loop(4) = {-11, 17:19};
Plane Surface(4) = {4};
Transfinite Surface{4};
Recombine Surface{4};

//bottom base moveable combs
Point(17) = {x, 0, 0 , 10*scale};
Point(18) = {0, 0, 0 , 10*scale}; // roulements
Point(19) = {0, 14*scale + Wss, 0 , 10*scale}; // roulements + limite avec suite
Point(20) = {x, 14*scale + Wss, 0 , 10*scale}; // limite avec suite
Line(20) = {15, 17};
Line(21) = {17, 18};
Line(22) = {18, 19}; // roulements
Line(23) = {19, 20}; // limite avec suite
Line(24) = {20, 16};
Transfinite Curve{20} = 4*nFEM+1 Using Progression 1;
Transfinite Curve{21, 23} = (5*g + 4.5*Wc)*nFEM+1 Using Progression 1;
Transfinite Curve{22} = (14+Ws)*nFEM+1 Using Progression 1;
Transfinite Curve{24} = 10*nFEM+1 Using Progression 1;
Curve Loop(5) = {-18, 20:24};
Plane Surface(5) = {5};
Transfinite Surface{5} = {17:20};
Recombine Surface{5};

//top base moveable combs
y = 26*scale+Wss;
Point(21) = {0, y, 0, nBEM*0.1*scale}; // roulements
Point(22) = {Wcs/2+gs, y, 0, nBEM*0.1*scale};
Point(23) = {Wcs/2+gs, y+Lcs, 0, nBEM*0.1*scale}; // top comb 1
Point(24) = {3*Wcs/2+gs, y+Lcs, 0, nBEM*0.1*scale}; // top comb 1
Point(25) = {3*Wcs/2+gs, y, 0, nBEM*0.1*scale};
Point(26) = {5*Wcs/2+3*gs, y, 0, nBEM*0.1*scale};
Point(27) = {5*Wcs/2+3*gs, y+Lcs, 0, nBEM*0.1*scale}; // top comb 2
Point(28) = {7*Wcs/2+3*gs, y+Lcs, 0, nBEM*0.1*scale}; // top comb 2
Point(29) = {7*Wcs/2+3*gs, y, 0, nBEM*0.1*scale};
Point(30) = {9*Wcs/2+5*gs, y, 0, nBEM*0.1*scale};

For i In {1:Ncombs-3}
    Point(30+4*i-3) = {9*Wcs/2+5*gs + (i-1)*(2*Wcs+2*gs), y+Lcs, 0, nBEM*0.1*scale}; // top comb
    Point(30+4*i-2) = {9*Wcs/2+5*gs + (i-1)*(2*Wcs+2*gs) + Wcs, y+Lcs, 0, nBEM*0.1*scale}; // top comb
    Point(30+4*i-1) = {9*Wcs/2+5*gs + (i-1)*(2*Wcs+2*gs) + Wcs, y, 0, nBEM*0.1*scale};
    Point(30+4*i) = {9*Wcs/2+5*gs + i*(2*Wcs+2*gs), y, 0, nBEM*0.1*scale};
EndFor

nbPts1 = 30 + (Ncombs-3)*4;
x = (5*g + 4.5*Wc)*scale + 2*(Wcs+gs);
Point(nbPts1 + 1) = {x, y+Lcs, 0, nBEM*0.1*scale}; // top last comb
Point(nbPts1 + 2) = {x+Wcs, y+Lcs, 0, nBEM*0.1*scale}; // top last comb
Point(nbPts1 + 3) = {x+Wcs, y, 0, nBEM*0.1*scale}; 
Point(nbPts1 + 4) = {x+Wcs+gs, y, 0, nBEM*0.1*scale}; // last point right
Point(nbPts1 + 5) = {x+Wcs+gs, y-12*scale, 0, nBEM*0.4*scale}; // last point right
nbPts2 = nbPts1 + 5;


Line(25) = {19, 21}; // roulements
Line(26) = {21, 22};
Line(27) = {22, 23};
Line(28) = {23, 24}; // top comb 1
Line(29) = {24, 25};
Line(30) = {25, 26};
Line(31) = {26, 27};
Line(32) = {27, 28}; // top comb 2
Line(33) = {28, 29}; 
Line(34) = {29, 30};

For i In {1:Ncombs-3}
    Line(34+4*i-3) = {30+4*i-4, 30+4*i-3};
    Line(34+4*i-2) = {30+4*i-3, 30+4*i-2}; // top comb
    Line(34+4*i-1) = {30+4*i-2, 30+4*i-1};
    Line(34+4*i) = {30+4*i-1, 30+4*i}; 
    Transfinite Curve{34+4*i-3, 34+4*i-1} = Lc*nFEM+1 Using Progression 1;
    Transfinite Curve{34+4*i-2} = Wc*nFEM+1 Using Progression 1;
    Transfinite Curve{34+4*i} = (Wc+2*g)*nFEM+1 Using Progression 1;
EndFor
nbLines1 = 34 + (Ncombs-3)*4;
Line(nbLines1 + 1) = {nbPts1, nbPts1 + 1};
Line(nbLines1 + 2) = {nbPts1 + 1, nbPts1 + 2}; //top last comb
Line(nbLines1 + 3) = {nbPts1 + 2, nbPts1 + 3};
Line(nbLines1 + 4) = {nbPts1 + 3, nbPts1 + 4};
Line(nbLines1 + 5) = {nbPts1 + 4, nbPts1 + 5}; // last line right
Line(nbLines1 + 6) = {nbPts1 + 5, 20};
nbLines2 = nbLines1 + 6;

Transfinite Curve{25, nbLines1 + 5} = 12*nFEM+1 Using Progression 1;
Transfinite Curve{26} = (Wc/2+g)*nFEM+1 Using Progression 1;
Transfinite Curve{27, 29, 31, 33, nbLines1 + 1, nbLines1 + 3} = Lc*nFEM+1 Using Progression 1;
Transfinite Curve{28, 32, nbLines1 + 2} = Wc*nFEM+1 Using Progression 1;
Transfinite Curve{30, 34} = (Wc+2*g)*nFEM+1 Using Progression 1;
Transfinite Curve{nbLines1 + 4} = g*nFEM+1 Using Progression 1;
Transfinite Curve{nbLines1 + 6} = ((Ncombs-2)*(2*Wc+2*g)-g-Wc)*nFEM+1 Using Progression 1;

// bottom lines
Line(nbLines2 + 1) = {25, 22}; // bottom comb 1
Line(nbLines2 + 2) = {29, 26}; // bottom comb 2
For i In {1:Ncombs-3}
    Line(nbLines2 + 2 + i) = {30+4*i-1, 30+4*i-4}; // bottom comb
EndFor
nbLines3 = nbLines2 + 2 + Ncombs-2;
Line(nbLines3) = {nbPts1 + 3, nbPts1}; //bottom last comb
Transfinite Curve{nbLines2 + 1: nbLines3} = Wc*nFEM+1 Using Progression 1;

// top base combs
Curve Loop(6) = {-23, 25, 26, -45,30,-46,34, -47, 38, -48, nbLines1+4:nbLines1+6}; // hardcoded
Plane Surface(6) = {6};
Transfinite Surface{6} = {19,21,nbPts2-1,nbPts2};
Recombine Surface{6};
// hardcoder le nombre de combs à partir d'ici
// four combs
Curve Loop(7) = {27, 28, 29, 45};
Plane Surface(7) = {7};
Curve Loop(8) = {31, 32, 33, 46};
Plane Surface(8) = {8};
Curve Loop(9) = {35, 36, 37, 47};
Plane Surface(9) = {9};
Curve Loop(10) = {39, 40, 41, 48};
Plane Surface(10) = {10};
Transfinite Surface{7:10};
Recombine Surface{7:10};

// mechanical physical surfaces
Physical Curve("anchor", 1) = {1:6};

Physical Surface("FEM_domain", 2) = {1:10};
Physical Curve("BEM_FEM_boundary", 3) = {26:nbLines2, 24, 19, 12:13};
Physical Curve("left_edge", 4) = {22,25};

//top electrode
y = y + Lcs + (Lcs-y0s) + 12*scale;
xtot = Ncombs*(2*gs+2*Wcs) + Wcs/2;
Point(nbPts2 + 1) = {0, y, 0, 10*scale};
Point(nbPts2 + 2) = {xtot, y, 0, 10*scale};
For i In {1:Ncombs}
    Point(nbPts2 + 2 + 4*i - 3) = {xtot-(i-1)*(2*Wcs+2*gs), y-12*scale-Lcs, 0, 0.2*nBEM*scale};
    Point(nbPts2 + 2 + 4*i - 2) = {xtot-(i-1)*(2*Wcs+2*gs) - Wcs, y-12*scale - Lcs, 0, 0.2*nBEM*scale};
    Point(nbPts2 + 2 + 4*i - 1) = {xtot-(i-1)*(2*Wcs+2*gs)-Wcs, y-12*scale, 0, 0.2*nBEM*scale};
    Point(nbPts2 + 2 + 4*i) = {xtot-i*(2*Wcs+2*gs), y-12*scale, 0, 0.2*nBEM*scale};
EndFor
nbPts3 = nbPts2 + 2 + 4*Ncombs;
Point(nbPts3 + 1) = {Wcs/2, y-12*scale-Lcs, 0, 0.2*nBEM*scale};
Point(nbPts3 + 2) = {0, y-12*scale-Lcs, 0, 0.2*nBEM*scale};
nbPts4 = nbPts3 + 2;

For i In {nbPts2+1:nbPts4-1}
    Line(nbLines3 + i - nbPts2) = {i, i+1};
EndFor
Line(nbLines3 + nbPts4 - nbPts2) = {nbPts4, nbPts2+1};
nbLines4 = nbLines3 + nbPts4 - nbPts2;

Physical Curve("top_electrode", 5) = {nbLines3+1:nbLines4-1};

// rest of outside
Line(nbLines4 + 1) = {nbPts4, 21};
y = -(Lg-2-4)*scale - 10*scale - Wss + 16*scale + 2*Wss + Lgs;
Point(nbPts4 + 1) = {Wtots, y, 0, 20*scale};
Point(nbPts4 + 2) = {Wtots, Htops, 0, 20*scale};
Point(nbPts4 + 3) = {0, Htops, 0, 20*scale};

Line(nbLines4 + 2) = {12, nbPts4 + 1};
Line(nbLines4 + 3) = {nbPts4 + 1, nbPts4 + 2};
Line(nbLines4 + 4) = {nbPts4 + 2, nbPts4 + 3};
Line(nbLines4 + 5) = {nbPts4 + 3, nbPts2 + 1};

Curve Loop(11) = {nbLines4 + 2: nbLines4 + 5, nbLines3 + 1: nbLines4 - 1, nbLines4 + 1, 26:nbLines2, 24, 19, 12:13};
Plane Surface(11) = {11};
Physical Surface("BEM_domain_1", 6) = {11};

Physical Curve("rest_of_outside", 7) = {nbLines4+1: nbLines4+5};