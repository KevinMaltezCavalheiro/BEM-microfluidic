/* -------------------------------------------------------------------
   File "microstrip.geo"

   This file is the geometrical description used by Gmsh to produce the mesh
   file "microstrip.msh", which is read by GetDP when analysing the file
   "microstrip.pro".
   ------------------------------------------------------------------- */

/* Definition of some parameters for geometrical dimensions, i.e.  h (height of
   'Diel1'), w (width of 'Line'), t (thickness of 'Line') xBox (width of the air
   box) and yBox (height of the air box) */

h = 1.e-3 ; w = 4.72e-3 ;  t = 0.035e-3 ;
xBox = w/2. * 6. ;  yBox = h * 12. ;

/* Definition of parameters for local mesh dimensions */

DefineConstant[ s = {1., Name "Parameters/Global mesh size factor"} ] ;

p0 = h / 10. * s ;
pLine0 = w/2. / 10. * s ;  pLine1 = w/2. / 50. * s ;
pxBox = xBox / 10. * s ;  pyBox = yBox / 8. * s ;

/* Definition of geometrical points */

Point(1) = { 0   , 0, 0, p0} ;
Point(2) = { xBox, 0, 0, pxBox} ;
Point(3) = { xBox, h, 0, pxBox} ;
Point(4) = { 0   , h, 0, pLine0} ;
Point(5) = { w/2., h, 0, pLine1} ;
Point(6) = { 0   , h+t, 0, pLine0} ;
Point(7) = { w/2., h+t, 0, pLine1} ;
Point(8) = { 0   , yBox, 0, pyBox} ;
Point(9) = { xBox, yBox, 0, pyBox} ;

/* Definition of gemetrical lines */

Line(1) = {1,2};   Line(2) = {2,3};  Line(3) = {3,9};
Line(4) = {9,8};   Line(5) = {8,6};  Line(7) = {4,1};
Line(8) = {5,3};   Line(9) = {4,5};  Line(10) = {6,7};
Line(11) = {5,7};

/* Definition of geometrical surfaces */

Curve Loop(12) = {1, 2, -8, -9, 7}; Plane Surface(13) = {12};
Curve Loop(14) = {10,-11,8,3,4,5};  Plane Surface(15) = {14};

/* Definition of Physical entities.  The definition of Physical entities
   (Surfaces and Curves) tells Gmsh the elements and associated region numbers
   that have to be saved in the mesh file 'microstrip.msh'.  For example, Region
   111 is made of the triangle elements of the geometric surface 13, whereas
   Region 121 is made of line elements of the geometric curves 9, 10 and 11. */

Physical Surface("Air", 101) = {15};
Physical Surface("Dielectric", 111) = {13};

Physical Curve("Ground", 120) = {1} ;
Physical Curve("Electrode", 121) = {9,10,11} ;
Physical Curve("Surface infinity", 130) = {2,3,4} ;
