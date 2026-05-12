/* -------------------------------------------------------------------
   Tutorial 1 : electrostatic field of a microstrip line

   Features:
   - Physical and abstract regions
   - Scalar FunctionSpace with Dirichlet constraint
   - Galerkin term for stiffness matrix

   To compute the solution in a terminal:
       getdp microstrip -solve EleSta_v
       getdp microstrip -pos Map
       getdp microstrip -pos Cut

   To compute the solution interactively from the Gmsh GUI:
       File > Open > microstrip.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/* In this first tutorial we consider the calculation of the electric field
   given a static distribution of electric potential. This corresponds to an
   "electrostatic" physical model, obtained by combining the time-invariant
   Faraday equation (curl e = 0, with e the electric field) with Gauss'
   law (div d = rho, with d the displacement field and rho the charge density)
   and the dielectric constitutive law (d = epsilon e, with epsilon the
   dielectric permittivity).

   Since curl e = 0, e can be derived from a scalar electric potential v, such
   that e = -grad v. Plugging this potential in Gauss' law and using the
   constitutive law leads to a scalar Poisson equation in terms of the electric
   potential: -div(epsilon grad v) = rho.

   We consider here the special case where rho = 0 to model a conducting
   microstrip line on top of a dielectric substrate. A Dirichlet boundary
   condition sets the potential to 1 mV on the boundary of the microstrip line
   (called "Electrode" below) and to 0 V on the ground. A homogeneous Neumann
   boundary condition (zero flux of the displacement field, i.e. n.d = 0) is
   imposed on the left boundary of the domain to account for the symmetry of the
   problem. The domain is truncated on the top and right with a homogeneous
   Dirichlet boundary condition (v = 0), assumed to be imposed sufficiently far
   away from the microstrip. */

Group {
  /* One starts by giving explicit meaningful names to the Physical regions
     defined in the "microstrip.msh" mesh file. This model comprises only
     2 volume regions and 2 surface regions. */

  Air = Region[101];
  Diel1 = Region[111];

  Ground = Region[120];
  Electrode = Region[121];
  SurfInf = Region[130];

  /* We now define abstract regions to be used below in the definition of the
     scalar electric potential formulation:

     Vol_Ele    : volume where -div(epsilon grad v) = 0 is solved
     Sur_Neu_Ele: surface where non homogeneous Neumann boundary conditions
                  (on n.d = -n . (epsilon grad v)) are imposed

     Vol_xxx groups contain only volume elements of the mesh (triangles here).
     Sur_xxx groups contain only surface elements of the mesh (lines here).

     Since there are no non-homogeneous Neumann conditions in this particular
     example, Sur_Neu_Ele is defined as empty.

     Note that volume elements are those that correspond to the higher dimension
     of the model at hand (2D elements here), surface elements correspond to the
     higher dimension of the model minus one (1D elements here).
     */

  Vol_Ele = Region[ {Air, Diel1} ];
  Sur_Neu_Ele = Region[ {} ];
}

Function {
  /* Material laws (here the dielectric permittivity) are defined as piecewise
     functions (note the square brackets), in terms of the above defined
     physical regions. */

  eps0 = 8.854187818e-12;
  epsilon[Air] = 1. * eps0;
  epsilon[Diel1] = 9.8 * eps0;
}

Constraint {
  /* The Dirichlet boundary condition is also defined piecewise.
     The constraint "Dirichlet_Ele" is invoked in the FunctionSpace
     below. */

  { Name Dirichlet_Ele; Type Assign;
    Case {
      { Region Ground; Value 0.; }
      { Region Electrode; Value 1.e-3; }
      { Region SurfInf; Value 0; }
    }
  }
}

Group{
  /* The domain of definition of a FunctionSpace lists all regions on which a
     field is defined.

     Contrary to Vol_xxx and Sur_xxx regions, which contain only volume or
     surface regions, resp., domains of definition Dom_xxx regions may contain
     both volume and surface regions. Hence the use of the prefixes Vol_, Sur_
     and Dom_ to avoid confusion.*/

  Dom_Hgrad_v_Ele =  Region[ {Vol_Ele, Sur_Neu_Ele} ];
}

FunctionSpace {
  /* The function space in which we pick the electric scalar potential "v"
     solution is defined by
     - a domain of definition (the "Support": "Dom_Hgrad_v_Ele")
     - a type ("Form0" means scalar field)
     - a set of basis functions ("BF_Node" means scalar nodal basis functions)
     - a set of entities to which the basis functions are associated ("Entity":
       here all the nodes of the domain of definition "NodesOf[All]")
     - a constraint (here the Dirichlet boundary conditions)

     The FE expansion of the unknown field "v" reads

     v(x,y) = Sum_k vn_k sn_k(x,y)

     where the "vn_k" coefficients are the nodal values (connectors) and
     "sn_k(x,y)" the nodal basis functions.  Not all connectors are unknowns of
     the FE problem, due to the "Constraint", which assigns particular values to
     the nodes of the Ground and Electrode regions. GetDP deals with that
     automatically on basis of the definition of the FunctionSpace. */

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Hgrad_v_Ele; Entity NodesOf[ All ]; }
      // using "NodesOf[All]" instead of "NodesOf[Dom_Hgrad_v_Ele]" is an
      // optimization, which avoids explicitly building the list of
      // all the nodes
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_Ele; }
    }
  }
}

Jacobian {
  /* Jacobians are used to specify the mapping between elements in the mesh and
     the reference elements (defined in standardized/reference unit cells) over which
     integration is performed. "Vol" represents the classical 1-to-1 mapping
     between identical spatial dimensions, i.e. in this case a reference
     triangle/quadrangle onto triangles/quadrangles in the z=0 plane (2D <->
     2D). "Sur" is used to map the reference triangle/quadrangle onto
     triangles/quadrangles in a 3D space (2D <-> 3D), or to map the reference
     line segment onto segments in 2D space (1D <-> 2D). "Lin" is used to map
     the reference line segment onto segments in 3D space (1D <-> 3D). */

  { Name Vol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
}

Integration {
  /* A Gauss quadrature rule with 4 points is used for all integrations below. */

  { Name Int ;
    Case { { Type Gauss ;
             Case { { GeoElement Line        ; NumberOfPoints  4 ; }
                    { GeoElement Triangle    ; NumberOfPoints  4 ; }
                    { GeoElement Quadrangle  ; NumberOfPoints  4 ; } }
      }
    }
  }
}

Formulation {
  /* The syntax of the Formulation section is a harder nut to crack. So let us
     deal with it carefully.

     A GetDP Formulation encodes a so-called weak formulation of the original
     partial differential equation, i.e. of -div(epsilon grad v) = 0. This weak
     formulation involves finding v such that

     (-div(epsilon grad v) , v')_Vol_Ele = 0

     holds for all so-called "test-functions" v', where (.,.)_D denotes an inner
     product over a domain D. If the test-functions v' are differentiable,
     integration by parts using Green's identity leads to finding v such that

     (epsilon grad v, grad v')_Vol_Ele - (n . (epsilon grad v), v')_Bnd_Vol_Ele = 0

     holds for all v', where Bnd_Vol_Ele is the boundary of Vol_Ele. In our
     microstrip example this surface term vanishes, as there is either no test
     function v' (on the Dirichlet boundary), or "n. (epsilon grad v)" is zero
     (on the homogeneous Neumann boundary). We are thus eventually looking for
     functions v in the function space Hgrad_v_Ele such that

     (epsilon grad v, grad v')_Vol_Ele = 0

     holds for all v'. Finally, our choice here is to use a Galerkin method,
     where the test functions v' are the same basis functions ("sn_k") as the
     ones used to interpolate the unknown potential v.

     The "Integral" statement in the Formulation is a symbolic representation of
     this weak formulation. It has 4 semicolon separated arguments:
     * the density [.,.] to be integrated (note the square brackets instead of
       the parentheses), with the test-functions (always) after the comma;
     * the domain of integration;
     * the Jacobian of the transformation reference element <-> real element;
     * the integration method to be used.

     In the density, braces around a quantity (such as "{v}") indicate that this
     quantity belongs to a FunctionSpace. Differential operators can be applied
     within braces (such as "{Grad v}"); in particular the symbol "d" represents
     the exterior derivative, and it is a synonym of "Grad" when applied to a
     scalar function (such as "{d v}").

     What can be a bit confusing is that the two comma-separated terms of the
     density are not interpreted exactly in the same way. Let us unravel this in
     detail.

     As the Galerkin method uses as test functions the same basis functions
     "sn_k" as for the unknown field "v", the second term in the density should
     be something like this [ ... , basis_functions_of {d v} ]. However, since
     the second term is always devoted to test functions, the operator
     "basis_functions_of" would always be there.  It can therefore be made
     implicit and, according to the GetDP syntax, omitted. So, one simply writes
     [ ... , {d v} ].

     On the other hand, the first term can contain a much wider variety of
     expressions than the second one. In our case it should be expressed in
     terms of the FE expansion of "v" at the present system solution, i.e. when
     the coefficients vn_k in the expansion of "v = Sum_k vn_k sn_k" are
     unknown. This is indicated by prefixing the braces with "Dof" (degrees of
     freedom), which leads to the following density:

     [ epsilon[] * Dof{d v} , {d v} ],

     a so-called bilinear term that will contribute to the stiffness matrix of
     the electrostatic problem at hand.

     Another option, which does not work here, is to evaluate the first
     argument with the last available computed solution, i.e. simply
     perform the interpolation with known coefficients vn_k. The Dof
     prefix is then omitted and we have:

     [ epsilon[] * {d v} , {d v} ],

     a so-called linear term that contributes to the right-hand side of the
     linear system.

     Both choices are commonly used in different contexts, and we shall often
     come back in subsequent tutorials. */

  { Name Electrostatics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      Integral { [ epsilon[] * Dof{d v} , {d v} ];
	    In Vol_Ele; Jacobian Vol; Integration Int; }

   /* Additional Integral terms can be added here. For example, the
      following term may account for non-homogeneous Neumann boundary
      conditions, provided that the function nd[] is defined:

      Integral { [ nd[] , {v} ];
	    In Sur_Neu_Ele; Jacobian Sur; Integration Int; }

      All the terms in the Equation environment are added,
      and an implicit "= 0" is considered at the end. */
    }
  }
}

/* In the Resolution environment we specify what to do with a weak formulation:
   here we simply generate a linear system, solve it and save the solution (.res file)
   to disk. */

Resolution {
  { Name EleSta_v;
    System {
      { Name Sys_Ele; NameOfFormulation Electrostatics_v; }
    }
    Operation {
      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
    }
  }
}

/* Post-processing is done in two parts.

   The first part defines, in terms of the Formulation, which itself refers to
   the FunctionSpace, a number of quantities that can be evaluated at the
   postprocessing stage. The three quantities defined here are:
   - the scalar electric potential;
   - the electric field;
   - the electric displacement.

   The second part consists in defining groups of post-processing operations,
   which can be invoked separately. The first group is invoked by default when
   Gmsh is run interactively. Each Operation specifies
   - a quantity to evaluate;
   - the region on which the evaluation is done;
   - the name of the output file.
   The generated post-processing files are automatically displayed by Gmsh if
   the "Merge result automatically" option is enabled (the default). */

PostProcessing {
  { Name EleSta_v; NameOfFormulation Electrostatics_v;
    Quantity {
      { Name v; Value {
          Term { [ {v} ]; In Dom_Hgrad_v_Ele; Jacobian Vol; }
        }
      }
      { Name e; Value {
          Term { [ -{d v} ]; In Dom_Hgrad_v_Ele; Jacobian Vol; }
        }
      }
      { Name d; Value {
          Term { [ -epsilon[] * {d v} ]; In Dom_Hgrad_v_Ele; Jacobian Vol; }
        }
      }
    }
  }
}

e = 1.e-7; // tolerance to ensure that the cut is inside the simulation domain
h = 2.e-3; // vertical position of the cut

PostOperation {
  { Name Map; NameOfPostProcessing EleSta_v;
     Operation {
       Print [ v, OnElementsOf Vol_Ele, File "mStrip_v.pos" ];
       Print [ e, OnElementsOf Vol_Ele, File "mStrip_e.pos" ];
       Print [ e, OnLine {{e,h,0}{14.e-3,h,0}}{60}, File "Cut_e.pos" ];
     }
  }
  { Name Cut; NameOfPostProcessing EleSta_v;
    // same cut as above, with more points and exported in raw text format
    Operation {
      Print [ e, OnLine {{e,e,0}{14.e-3,e,0}} {500}, Format TimeTable, File "Cut_e.txt" ];
    }
  }
}
