scale = 2e-6;

// USE WITH MINIMUM 2 FINS, else use longitudinal_comb.geo
N_fins = 4; // number of fins on one side of the comb // MAX 12 si on change ps la longueur

// WARNING: when using more fins the pull-in voltage decreases

n = 6; // FEM elements density
nBEM = 1; // BEM elements density

a = 0; // acceleration verticale

// mechanical properties and boundary conditions
SetNumber("Boundary Conditions/left/ux", 0.); // encastrement
SetNumber("Boundary Conditions/left/uy", 0);
SetNumber("Boundary Conditions/right/ux", 0.); // encastrement
SetNumber("Boundary Conditions/right/uy", 0);
SetNumber("Materials/FEM_domain/Young", 150e9);
SetNumber("Materials/FEM_domain/Poisson", 0.27);
SetNumber("Materials/FEM_domain/rho",2300); 
SetNumber("Volumic Forces/FEM_domain/bx",0);
SetNumber("Volumic Forces/FEM_domain/by",a); // acceleration of accelerometer

SetNumber("Non_linear_solver", 1);

// BEM properties in bottom domain
phi_1 = 20;
SetNumber("Boundary Conditions/BEM_FEM_boundary_1/BEM_domain_1/dirichlet", 0);
SetNumber("Boundary Conditions/electrode_1/BEM_domain_1/dirichlet", phi_1);
SetNumber("Boundary Conditions/outside_1/BEM_domain_1/neumann", 0);
SetNumber("Materials/BEM_domain_1/Epsilon", 8.8541878128e-12); // dielectric permittivity

// BEM properties in top domain
phi_2 = 20;
SetNumber("Boundary Conditions/BEM_FEM_boundary_2/BEM_domain_2/dirichlet", 0);
SetNumber("Boundary Conditions/electrode_2/BEM_domain_2/dirichlet", phi_2);
SetNumber("Boundary Conditions/outside_2/BEM_domain_2/neumann", 0);
SetNumber("Materials/BEM_domain_2/Epsilon", 8.8541878128e-12); // dielectric permittivity

h_tot = 30*scale;
h_base = 0.4*scale;
h_fin = 2.8*scale;
h_space = 1*scale; // space between bout de l'electrode and base of the clamped beam of the comb
h_fin_elec = 2.4*scale; // length of the fins of the electrode, can be longer than the classical fins
// WARNING, abs(h_fin_elec - h_fin) must be lower than h_space

//l_bord = 10*scale;
//l_tot = 34.4*scale;
l_tot = 60*scale;
l_fin = 0.8*scale;
l_space = 0.8*scale;
t_electrode = 1*scale; // width of one electrode
l_periodic = l_fin + 2*l_space + t_electrode; // distance between two fins
//l_tot = 2*l_bord + l_fin + (N_fins-1)*l_periodic;
l_bord = (l_tot - l_fin - (N_fins-1)*l_periodic)/2;

unit_l = 0.2*scale; // reference for the transfinite curves

// définition des points du contour
Point(1) = {0, -h_base/2, 0, nBEM*scale};
Point(2) = {0, -h_tot/2, 0, nBEM*scale};
Point(3) = {l_tot, -h_tot/2, 0, nBEM*scale};
Point(4) = {l_tot, -h_base/2, 0, nBEM*scale};
Point(5) = {l_tot, h_base/2, 0, nBEM*scale};
Point(6) = {l_tot, h_tot/2, 0, nBEM*scale};
Point(7) = {0, h_tot/2, 0, nBEM*scale};
Point(8) = {0, h_base/2, 0, nBEM*scale};
// définition des lignes du contour
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Transfinite Curve {4,8} = (h_base/unit_l)*n+1 Using Progression 1;


// définition des points des fins du bas
offsetp1 = 8;
For i In {1:N_fins}
    Point(offsetp1+4*i-3) = {(i-1)*l_periodic + l_bord, -h_base/2, 0, 0.2*scale};
    Point(offsetp1+4*i-2) = {(i-1)*l_periodic + l_bord, -h_base/2 - h_fin, 0, 0.2*scale};
    Point(offsetp1+4*i-1) = {(i-1)*l_periodic + l_bord + l_fin, -h_base/2 - h_fin, 0, 0.2*scale};
    Point(offsetp1+4*i) = {(i-1)*l_periodic + l_bord + l_fin, -h_base/2, 0, 0.2*scale};
EndFor
// définition des lignes des fins du bas
Line(9) = {offsetp1+1, 1};
Transfinite Curve {9} = (l_bord/unit_l)*n+1 Using Progression 1;
offsetl1 = 9;
For i In {1:N_fins-1}
    Line(offsetl1+4*i-3) = {offsetp1+4*i-2, offsetp1+4*i-3};
    Line(offsetl1+4*i-2) = {offsetp1+4*i-1, offsetp1+4*i-2};
    Line(offsetl1+4*i-1) = {offsetp1+4*i, offsetp1+4*i-1};
    Line(offsetl1+4*i) = {offsetp1+4*i+1, offsetp1+4*i};
    Transfinite Curve {offsetl1+4*i-3, offsetl1+4*i-1} = (h_fin/unit_l)*n+1 Using Progression 1;
    Transfinite Curve {offsetl1+4*i-2} = (l_fin/unit_l)*n+1 Using Progression 1;
    Transfinite Curve {offsetl1+4*i} = ((l_periodic-l_fin)/unit_l)*n+1 Using Progression 1;
EndFor
Line(offsetl1+4*N_fins-3) = {offsetp1+4*N_fins-2, offsetp1+4*N_fins-3};
Line(offsetl1+4*N_fins-2) = {offsetp1+4*N_fins-1, offsetp1+4*N_fins-2};
Line(offsetl1+4*N_fins-1) = {offsetp1+4*N_fins, offsetp1+4*N_fins-1};
Line(offsetl1+4*N_fins) = {4, offsetp1+4*N_fins};
Transfinite Curve {offsetl1+4*N_fins-3, offsetl1+4*N_fins-1} = (h_fin/unit_l)*n+1 Using Progression 1;
Transfinite Curve {offsetl1+4*N_fins-2} = (l_fin/unit_l)*n+1 Using Progression 1;
Transfinite Curve {offsetl1+4*N_fins} = (l_bord/unit_l)*n+1 Using Progression 1;
// définition des lignes entre les fins du bas et la poutre principale
For i In{1:N_fins}
    Line(offsetl1+4*N_fins+i) = {offsetp1+4*i-3, offsetp1+4*i};
    Transfinite Curve {offsetl1+4*N_fins+i} = (l_fin/unit_l)*n+1 Using Progression 1;
EndFor

// définition des points des fins du haut
offsetp2 = offsetp1 + 4*N_fins;
For i In {1:N_fins}
    Point(offsetp2+4*i-3) = {(i-1)*l_periodic + l_bord, h_base/2, 0, 0.2*scale};
    Point(offsetp2+4*i-2) = {(i-1)*l_periodic + l_bord, h_base/2 + h_fin, 0, 0.2*scale};
    Point(offsetp2+4*i-1) = {(i-1)*l_periodic + l_bord + l_fin, +h_base/2 + h_fin, 0, 0.2*scale};
    Point(offsetp2+4*i) = {(i-1)*l_periodic + l_bord + l_fin, +h_base/2, 0, 0.2*scale};
EndFor
// définition des lignes des fins du haut
offsetl2 = offsetl1 + 5*N_fins + 1;
Line(offsetl2) = {8, offsetp2+1};
Transfinite Curve {offsetl2} = (l_bord/unit_l)*n+1 Using Progression 1;
For i In {1:N_fins-1}
    Line(offsetl2+4*i-3) = {offsetp2+4*i-3, offsetp2+4*i-2};
    Line(offsetl2+4*i-2) = {offsetp2+4*i-2, offsetp2+4*i-1};
    Line(offsetl2+4*i-1) = {offsetp2+4*i-1, offsetp2+4*i};
    Line(offsetl2+4*i) = {offsetp2+4*i, offsetp2+4*i+1};
    Transfinite Curve {offsetl2+4*i-3, offsetl2+4*i-1} = (h_fin/unit_l)*n+1 Using Progression 1;
    Transfinite Curve {offsetl2+4*i-2} = (l_fin/unit_l)*n+1 Using Progression 1;
    Transfinite Curve {offsetl2+4*i} = ((l_periodic-l_fin)/unit_l)*n+1 Using Progression 1;
EndFor
Line(offsetl2+4*N_fins-3) = {offsetp2+4*N_fins-3, offsetp2+4*N_fins-2};
Line(offsetl2+4*N_fins-2) = {offsetp2+4*N_fins-2, offsetp2+4*N_fins-1};
Line(offsetl2+4*N_fins-1) = {offsetp2+4*N_fins-1, offsetp2+4*N_fins};
Line(offsetl2+4*N_fins) = {offsetp2+4*N_fins, 5};
Transfinite Curve {offsetl2+4*N_fins-3, offsetl2+4*N_fins-1} = (h_fin/unit_l)*n+1 Using Progression 1;
Transfinite Curve {offsetl2+4*N_fins-2} = (l_fin/unit_l)*n+1 Using Progression 1;
Transfinite Curve {offsetl2+4*N_fins} = (l_bord/unit_l)*n+1 Using Progression 1;
// définition des lignes entre les fins du haut et la poutre principale
For i In{1:N_fins}
    Line(offsetl2+4*N_fins+i) = {offsetp2+4*i, offsetp2+4*i-3};
    Transfinite Curve {offsetl2+4*N_fins+i} = (l_fin/unit_l)*n+1 Using Progression 1;
EndFor

// définition des lignes dans la poutre centrale pour définir sous-domaines
offsetl3 = offsetl2 + 5*N_fins;
For i In {1:N_fins}
    Line(offsetl3+2*i-1) = {offsetp2+4*i-3, offsetp1+4*i-3};
    Line(offsetl3+2*i) = {offsetp1+4*i, offsetp2+4*i};
    Transfinite Curve {offsetl3+2*i-1, offsetl3+2*i} = (h_base/unit_l)*n+1 Using Progression 1;
EndFor

// définition des courbes et surfaces dans les fins du bas
For i In {1:N_fins}
    Curve Loop(i) = {offsetl1+4*i-1, offsetl1+4*i-2, offsetl1+4*i-3, offsetl1+4*N_fins+i};
    Plane Surface(i) = {i};
    Recombine Surface {i};
EndFor

// définition des courbes et surfaces dans les fins du haut
offsets1 = N_fins;
For i In {1:N_fins}
    Curve Loop(offsets1 + i) = {offsetl2+4*i-3, offsetl2+4*i-2, offsetl2+4*i-1, offsetl2+4*N_fins+i};
    Plane Surface(offsets1 + i) = {offsets1 + i};
    Recombine Surface {offsets1 + i};
EndFor

// définition des courbes et surfaces dans la poutre centrale
offsets2 = offsets1 + N_fins + 1;
Curve Loop(offsets2) = {-8, offsetl2, offsetl3+1, offsetl1};
Plane Surface(offsets2) = {offsets2};
Recombine Surface {offsets2};
For i In {1:N_fins-1}
    Curve Loop(offsets2 + 2*i - 1) = {offsetl3+2*i-1, offsetl1+4*N_fins+i, offsetl3+2*i, offsetl2+4*N_fins+i};
    Plane Surface(offsets2 + 2*i - 1) = {offsets2 + 2*i - 1};
    Recombine Surface {offsets2 + 2*i - 1};
    Curve Loop(offsets2 + 2*i) = {offsetl3+2*i, offsetl2+4*i, offsetl3+2*i+1, offsetl1+4*i};
    Plane Surface(offsets2 + 2*i) = {offsets2 + 2*i};
    Recombine Surface {offsets2 + 2*i};
EndFor
Curve Loop(offsets2 + 2*N_fins - 1) = {offsetl3+2*N_fins-1, offsetl1+4*N_fins+N_fins, offsetl3+2*N_fins, offsetl2+4*N_fins+N_fins};
Plane Surface(offsets2 + 2*N_fins - 1) = {offsets2 + 2*N_fins - 1};
Recombine Surface {offsets2 + 2*N_fins - 1};
Curve Loop(offsets2 + 2*N_fins) = {offsetl3+2*N_fins, offsetl2+4*N_fins, -4, offsetl1+4*N_fins};
Plane Surface(offsets2 + 2*N_fins) = {offsets2 + 2*N_fins};
Recombine Surface {offsets2 + 2*N_fins};


// mechanical physical surfaces
Physical Curve("left", 1) = {8};
Physical Curve("right", 2) = {4};

Physical Surface("FEM_domain", 3) = {1:offsets2 + 2*N_fins};
Transfinite Surface {1:offsets2 + 2*N_fins};
Physical Curve("BEM_FEM_boundary", 4) = {offsetl1:offsetl1+4*N_fins, offsetl2:offsetl2+4*N_fins};

// définition des points de la bottom electrode
x0_e = l_bord - l_space - t_electrode; // pour faciliter le bazar
y0_e = h_base/2 + h_space;
offsetp3 = offsetp2 + 4*N_fins;
Point(offsetp3 + 1) = {l_tot - x0_e - t_electrode, -y0_e - h_fin_elec, 0, nBEM*0.2*scale};
Point(offsetp3 + 2) = {l_tot - x0_e - t_electrode, -y0_e, 0, nBEM*0.2*scale};
Point(offsetp3 + 3) = {l_tot - x0_e, -y0_e, 0, nBEM*0.2*scale};
Point(offsetp3 + 4) = {l_tot - x0_e, -y0_e - h_fin_elec - t_electrode, 0, nBEM*0.2*scale};
Point(offsetp3 + 5) = {x0_e, -y0_e - h_fin_elec - t_electrode, 0, nBEM*0.2*scale};
Point(offsetp3 + 6) = {x0_e, -y0_e, 0, nBEM*0.2*scale};
Point(offsetp3 + 7) = {x0_e + t_electrode, -y0_e, 0, nBEM*0.2*scale};
Point(offsetp3 + 8) = {x0_e + t_electrode, -y0_e - h_fin_elec, 0, nBEM*0.2*scale};

offsetp4 = offsetp3 + 8;
For i In {1:N_fins-1}
    Point(offsetp4 + 4*i - 3) = {x0_e + l_periodic*i, -y0_e - h_fin_elec, 0, nBEM*0.2*scale};
    Point(offsetp4 + 4*i - 2) = {x0_e + l_periodic*i, -y0_e, 0, nBEM*0.2*scale};
    Point(offsetp4 + 4*i - 1) = {x0_e + l_periodic*i + t_electrode, -y0_e, 0, nBEM*0.2*scale};
    Point(offsetp4 + 4*i) = {x0_e + l_periodic*i + t_electrode, -y0_e - h_fin_elec, 0, nBEM*0.2*scale};
EndFor

// définition des lignes générales (indépendantes de N_fins) de la bottom electrode
offsetl4 = offsetl3 + N_fins*2;
Line(offsetl4 + 1) = {offsetp3 + 1, offsetp3 + 2};
Line(offsetl4 + 2) = {offsetp3 + 2, offsetp3 + 3};
Line(offsetl4 + 3) = {offsetp3 + 3, offsetp3 + 4};
Line(offsetl4 + 4) = {offsetp3 + 4, offsetp3 + 5};
Line(offsetl4 + 5) = {offsetp3 + 5, offsetp3 + 6};
Line(offsetl4 + 6) = {offsetp3 + 6, offsetp3 + 7};
Line(offsetl4 + 7) = {offsetp3 + 7, offsetp3 + 8};
Line(offsetl4 + 8) = {offsetp3 + 8, offsetp4 + 1};
offsetl5 = offsetl4 + 8;
// définition des lignes fins electrode du bas
For i In {1:N_fins-2}
    Line(offsetl5 + 4*i - 3) = {offsetp4 + 4*i - 3, offsetp4 + 4*i - 2};
    Line(offsetl5 + 4*i - 2) = {offsetp4 + 4*i - 2, offsetp4 + 4*i - 1};
    Line(offsetl5 + 4*i - 1) = {offsetp4 + 4*i - 1, offsetp4 + 4*i};
    Line(offsetl5 + 4*i) = {offsetp4 + 4*i, offsetp4 + 4*i + 1};
EndFor
Line(offsetl5 + 4*(N_fins-1) - 3) = {offsetp4 + 4*(N_fins-1) - 3, offsetp4 + 4*(N_fins-1) - 2};
Line(offsetl5 + 4*(N_fins-1) - 2) = {offsetp4 + 4*(N_fins-1) - 2, offsetp4 + 4*(N_fins-1) - 1};
Line(offsetl5 + 4*(N_fins-1) - 1) = {offsetp4 + 4*(N_fins-1) - 1, offsetp4 + 4*(N_fins-1)};
Line(offsetl5 + 4*(N_fins-1)) = {offsetp4 + 4*(N_fins-1), offsetp3+1};

// définition des curve loop et surface du BEM_domain1
offsets3 = offsets2 + 2*N_fins;
Curve Loop(offsets3 + 1) = {offsetl1+4*N_fins:offsetl1:-1,1,2,3}; // exterior boundary of BEM_domain_1
Curve Loop(offsets3 + 2) = {offsetl4 + 1 : offsetl5 + 4*(N_fins-1)};
Plane Surface(offsets3 + 1) = {offsets3 + 1, offsets3 + 2};

Physical Surface("BEM_domain_1", 5) = {offsets3 + 1};
Physical Curve("BEM_FEM_boundary_1", 6) = {offsetl1+4*N_fins:offsetl1:-1};
Physical Curve("electrode_1", 7) = {offsetl4 + 1 : offsetl5 + 4*(N_fins-1)};
Physical Curve("outside_1", 8) = {1, 2, 3};

// définition des points de la top electrode
offsetp5 = offsetp4 + 4*(N_fins-1);
Point(offsetp5 + 1) = {x0_e + t_electrode, y0_e + h_fin_elec, 0, nBEM*0.2*scale};
Point(offsetp5 + 2) = {x0_e + t_electrode, y0_e, 0, nBEM*0.2*scale};
Point(offsetp5 + 3) = {x0_e, y0_e, 0, nBEM*0.2*scale};
Point(offsetp5 + 4) = {x0_e, y0_e + h_fin_elec + t_electrode, 0, nBEM*0.2*scale};
Point(offsetp5 + 5) = {l_tot - x0_e, y0_e + h_fin_elec + t_electrode, 0, nBEM*0.2*scale};
Point(offsetp5 + 6) = {l_tot - x0_e, y0_e, 0, nBEM*0.2*scale};
Point(offsetp5 + 7) = {l_tot - x0_e - t_electrode, y0_e, 0, nBEM*0.2*scale};
Point(offsetp5 + 8) = {l_tot - x0_e - t_electrode, y0_e + h_fin_elec, 0, nBEM*0.2*scale};

offsetp6 = offsetp5 + 8;
For i In {1:N_fins-1}
    Point(offsetp6 + 4*i - 3) = {x0_e + l_periodic*i, y0_e + h_fin_elec, 0, nBEM*0.2*scale};
    Point(offsetp6 + 4*i - 2) = {x0_e + l_periodic*i, y0_e, 0, nBEM*0.2*scale};
    Point(offsetp6 + 4*i - 1) = {x0_e + l_periodic*i + t_electrode, y0_e, 0, nBEM*0.2*scale};
    Point(offsetp6 + 4*i) = {x0_e + l_periodic*i + t_electrode, y0_e + h_fin_elec, 0, nBEM*0.2*scale};
EndFor

// définition des lignes générales (indépendantes de N_fins) de la top electrode
offsetl6 = offsetl5 + 4*(N_fins-1);
Line(offsetl6 + 1) = {offsetp5 + 1, offsetp5 + 2};
Line(offsetl6 + 2) = {offsetp5 + 2, offsetp5 + 3};
Line(offsetl6 + 3) = {offsetp5 + 3, offsetp5 + 4};
Line(offsetl6 + 4) = {offsetp5 + 4, offsetp5 + 5};
Line(offsetl6 + 5) = {offsetp5 + 5, offsetp5 + 6};
Line(offsetl6 + 6) = {offsetp5 + 6, offsetp5 + 7};
Line(offsetl6 + 7) = {offsetp5 + 7, offsetp5 + 8};
Line(offsetl6 + 8) = {offsetp5 + 8, offsetp6 + 4*(N_fins-1)};
offsetl7 = offsetl6 + 8;
// définition des lignes fins electrode du haut
Line(offsetl7 + 1) = {offsetp6 + 4 - 3, offsetp5 + 1};
Line(offsetl7 + 2) = {offsetp6 + 4 - 2, offsetp6 + 4 - 3};
Line(offsetl7 + 3) = {offsetp6 + 4 - 1, offsetp6 + 4 - 2};
Line(offsetl7 + 4) = {offsetp6 + 4, offsetp6 + 4 - 1};
For i In {2:N_fins-1}
    Line(offsetl7 + 4*i - 3) = {offsetp6 + 4*i - 3, offsetp6 + 4*i - 4};
    Line(offsetl7 + 4*i - 2) = {offsetp6 + 4*i - 2, offsetp6 + 4*i - 3};
    Line(offsetl7 + 4*i - 1) = {offsetp6 + 4*i - 1, offsetp6 + 4*i - 2};
    Line(offsetl7 + 4*i) = {offsetp6 + 4*i, offsetp6 + 4*i - 1};
EndFor

// définition des curve loop et surface du BEM_domain2
Curve Loop(offsets3 + 3) = {offsetl2:offsetl2+4*N_fins,5,6,7}; // exterior boundary of BEM_domain_2
Curve Loop(offsets3 + 4) = {offsetl6 + 1 : offsetl6 + 8, offsetl7 + 4*(N_fins-1):offsetl7+1:-1};
Plane Surface(offsets3 + 2) = {offsets3 + 3, offsets3 + 4};

Physical Surface("BEM_domain_2", 9) = {offsets3 + 2};
Physical Curve("BEM_FEM_boundary_2", 10) = {offsetl2:offsetl2+4*N_fins};
Physical Curve("electrode_2", 11) = {offsetl6 + 1 : offsetl6 + 8, offsetl7 + 4*(N_fins-1):offsetl7+1:-1};
Physical Curve("outside_2", 12) = {5,6,7};
