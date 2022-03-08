#include <gmsh.h>
#include <iostream>
#include <sstream>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Core> 
#include <Eigen/Sparse> 
#include <set>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <vector>
//#include <QImage>
#include "test.h"

using namespace std;
using namespace gmsh;
using namespace Eigen;

typedef Triplet<double> T;

void cheeseHole(double x, double y, double z, double r, double lc, vector<int> &shells, vector<int> &volumes)
{
 // This function will create a spherical hole in a volume. We don't specify
 // tags manually, and let the functions return them automatically:

    int p1 = model::geo::addPoint(x, y, z, lc);
    int p2 = model::geo::addPoint(x + r, y, z, lc);
    int p3 = model::geo::addPoint(x, y + r, z, lc);
    int p4 = model::geo::addPoint(x, y, z + r, lc);
    int p5 = model::geo::addPoint(x - r, y, z, lc);
    int p6 = model::geo::addPoint(x, y - r, z, lc);
    int p7 = model::geo::addPoint(x, y, z - r, lc);

    int c1 = model::geo::addCircleArc(p2, p1, p7);
    int c2 = model::geo::addCircleArc(p7, p1, p5);
    int c3 = model::geo::addCircleArc(p5, p1, p4);
    int c4 = model::geo::addCircleArc(p4, p1, p2);
    int c5 = model::geo::addCircleArc(p2, p1, p3);
    int c6 = model::geo::addCircleArc(p3, p1, p5);
    int c7 = model::geo::addCircleArc(p5, p1, p6);
    int c8 = model::geo::addCircleArc(p6, p1, p2);
    int c9 = model::geo::addCircleArc(p7, p1, p3);
    int c10 = model::geo::addCircleArc(p3, p1, p4);
    int c11 = model::geo::addCircleArc(p4, p1, p6);
    int c12 = model::geo::addCircleArc(p6, p1, p7);

    int l1 = model::geo::addCurveLoop({c5, c10, c4});
    int l2 = model::geo::addCurveLoop({c9, -c5, c1});
    int l3 = model::geo::addCurveLoop({c12, -c8, -c1});
    int l4 = model::geo::addCurveLoop({c8, -c4, c11});
    int l5 = model::geo::addCurveLoop({-c10, c6, c3});
    int l6 = model::geo::addCurveLoop({-c11, -c3, c7});
    int l7 = model::geo::addCurveLoop({-c2, -c7, -c12});
    int l8 = model::geo::addCurveLoop({-c6, -c9, c2});

    // We need non-plane surfaces to define the spherical holes. Here we use the
    // `model::geo::addSurfaceFilling()' function, which can be used for
    // surfaces with 3 or 4 curves on their boundary. With the he built-in kernel,
    // if the curves are circle arcs, ruled surfaces are created; otherwise
    // transfinite interpolation is used.
    //
    // With the OpenCASCADE kernel, `model::occ::addSurfaceFilling()' uses a
    // much more general generic surface filling algorithm, creating a BSpline
    // surface passing through an arbitrary number of boundary curves. The
    // `model::geo::addThruSections()' allows to create ruled surfaces (see
    // `t19.cpp').

    int s1 = model::geo::addSurfaceFilling({l1});
    int s2 = model::geo::addSurfaceFilling({l2});
    int s3 = model::geo::addSurfaceFilling({l3});
    int s4 = model::geo::addSurfaceFilling({l4});
    int s5 = model::geo::addSurfaceFilling({l5});
    int s6 = model::geo::addSurfaceFilling({l6});
    int s7 = model::geo::addSurfaceFilling({l7});
    int s8 = model::geo::addSurfaceFilling({l8});

    int sl = model::geo::addSurfaceLoop({s1, s2, s3, s4, s5, s6, s7, s8});
    int v = model::geo::addVolume({sl});
    shells.push_back(sl);
    volumes.push_back(v);
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 1
//
//  Geometry and mesh data
//
// -----------------------------------------------------------------------------

// The C++ API allows to do much more than what can be done in .geo files. These
// additional features are introduced gradually in the extended tutorials,
// starting with `x1.cpp'.

// In this first extended tutorial, we start by using the API to access basic
// geometrical and mesh data.
void x1(int argc, char **argv)
{
    // You can run this tutorial on any file that Gmsh can read, e.g. a mesh file
  // in the MSH format: `t1.exe file.msh'
  open(argv[1]);

  // Print the model name and dimension:
  string name;
  model::getCurrent(name);
  cout << "Model " << name << " (" << model::getDimension()
            << "D)\n";

  // Geometrical data is made of elementary model `entities', called `points'
  // (entities of dimension 0), `curves' (entities of dimension 1), `surfaces'
  // (entities of dimension 2) and `volumes' (entities of dimension 3). As we
  // have seen in the other C++ tutorials, elementary model entities are
  // identified by their dimension and by a `tag': a strictly positive
  // identification number. Model entities can be either CAD entities (from the
  // built-in `geo' kernel or from the OpenCASCADE `occ' kernel) or `discrete'
  // entities (defined by a mesh). `Physical groups' are collections of model
  // entities and are also identified by their dimension and by a tag.

  // Get all the elementary entities in the model, as a vector of (dimension,
  // tag) pairs:
  vector<pair<int, int> > entities;
  model::getEntities(entities);

  for(auto e : entities) {
    // Dimension and tag of the entity:
    int dim = e.first, tag = e.second;

    // Mesh data is made of `elements' (points, lines, triangles, ...), defined
    // by an ordered list of their `nodes'. Elements and nodes are identified by
    // `tags' as well (strictly positive identification numbers), and are stored
    // ("classified") in the model entity they discretize. Tags for elements and
    // nodes are globally unique (and not only per dimension, like entities).

    // A model entity of dimension 0 (a geometrical point) will contain a mesh
    // element of type point, as well as a mesh node. A model curve will contain
    // line elements as well as its interior nodes, while its boundary nodes
    // will be stored in the bounding model points. A model surface will contain
    // triangular and/or quadrangular elements and all the nodes not classified
    // on its boundary or on its embedded entities. A model volume will contain
    // tetrahedra, hexahedra, etc. and all the nodes not classified on its
    // boundary or on its embedded entities.

    // Get the mesh nodes for the entity (dim, tag):
    vector<size_t> nodeTags;
    vector<double> nodeCoords, nodeParams;
    model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

    // Get the mesh elements for the entity (dim, tag):
    vector<int> elemTypes;
    vector<vector<size_t> > elemTags, elemNodeTags;
    model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

    // Elements can also be obtained by type, by using `getElementTypes()'
    // followed by `getElementsByType()'.

    // Let's print a summary of the information available on the entity and its
    // mesh.

    // * Type of the entity:
    string type;
    model::getType(dim, tag, type);
    string name;
    model::getEntityName(dim, tag, name);
    if(name.size()) name += " ";
    cout << "Entity " << name << "(" << dim << "," << tag << ") of type "
              << type << "\n";

    // * Number of mesh nodes and elements:
    int numElem = 0;
    for(auto &tags : elemTags) numElem += tags.size();
    cout << " - Mesh has " << nodeTags.size() << " nodes and " << numElem
              << " elements\n";

    // * Upward and downward adjacencies:
    vector<int> up, down;
    model::getAdjacencies(dim, tag, up, down);
    if(up.size()) {
      cout << " - Upward adjacencies: ";
      for(auto e : up) cout << e << " ";
      cout << "\n";
    }
    if(down.size()) {
      cout << " - Downward adjacencies: ";
      for(auto e : down) cout << e << " ";
      cout << "\n";
    }

    // * Does the entity belong to physical groups?
    vector<int> physicalTags;
    model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
    if(physicalTags.size()) {
      cout << " - Physical group: ";
      for(auto physTag : physicalTags) {
        string n;
        model::getPhysicalName(dim, physTag, n);
        if(n.size()) n += " ";
        cout << n << "(" << dim << ", " << physTag << ") ";
      }
      cout << "\n";
    }

    // * Is the entity a partition entity? If so, what is its parent entity?
    vector<int> partitions;
    model::getPartitions(dim, tag, partitions);
    if(partitions.size()) {
      cout << " - Partition tags:";
      for(auto part : partitions) cout << " " << part;
      int parentDim, parentTag;
      model::getParent(dim, tag, parentDim, parentTag);
      cout << " - parent entity (" << parentDim << "," << parentTag
                << ")\n";
    }

    // * List all types of elements making up the mesh of the entity:
    for(auto elemType : elemTypes) {
      string name;
      int d, order, numv, numpv;
      vector<double> param;
      model::mesh::getElementProperties(elemType, name, d, order, numv,
                                              param, numpv);
      cout << " - Element type: " << name << ", order " << order << "\n";
      cout << "   with " << numv << " nodes in param coord: (";
      for(auto p : param) cout << p << " ";
      cout << ")\n";
    }
  }

  // We can use this to clear all the model data:
  clear();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 2
//
//  Mesh import, discrete entities, hybrid models, terrain meshing
//
// -----------------------------------------------------------------------------
// The API can be used to import a mesh without reading it from a file, by
// creating nodes and elements on the fly and storing them in model
// entities. These model entities can be existing CAD entities, or can be
// discrete entities, entirely defined by the mesh.
//
// Discrete entities can be reparametrized (see `t13.py') so that they can be
// remeshed later on; and they can also be combined with built-in CAD entities
// to produce hybrid models.
//
// We combine all these features in this tutorial to perform terrain meshing,
// where the terrain is described by a discrete surface (that we then
// reparametrize) combined with a CAD representation of the underground.
void x2(int argc, char **argv)
{
    model::add("x2");

  // We will create the terrain surface mesh from N x N input data points:
  int N = 100;

  // Helper function to return a node tag given two indices i and j:
  auto tag = [N](size_t i, size_t j) { return (N + 1) * i + j + 1; };

  // The x, y, z coordinates of all the points:
  vector<double> coords;

  // The tags of the corresponding nodes:
  vector<size_t> nodes;

  // The connectivities of the triangle elements (3 node tags per triangle) on
  // the terrain surface:
  vector<size_t> tris;

  // The connectivities of the line elements on the 4 boundaries (2 node tags
  // for each line element):
  vector<size_t> lin[4];

  // The connectivities of the point elements on the 4 corners (1 node tag for
  // each point element):
  size_t pnt[4] = {tag(0, 0), tag(N, 0), tag(N, N), tag(0, N)};

  for(size_t i = 0; i < N + 1; i++) {
    for(size_t j = 0; j < N + 1; j++) {
      nodes.push_back(tag(i, j));
      coords.insert(coords.end(), {(double)i / N, (double)j / N,
                                   0.05 * sin(10 * (double)(i + j) / N)});
      if(i > 0 && j > 0) {
        tris.insert(tris.end(),
                    {tag(i - 1, j - 1), tag(i, j - 1), tag(i - 1, j)});
        tris.insert(tris.end(), {tag(i, j - 1), tag(i, j), tag(i - 1, j)});
      }
      if((i == 0 || i == N) && j > 0) {
        size_t s = (i == 0) ? 3 : 1;
        lin[s].insert(lin[s].end(), {tag(i, j - 1), tag(i, j)});
      }
      if((j == 0 || j == N) && i > 0) {
        size_t s = (j == 0) ? 0 : 2;
        lin[s].insert(lin[s].end(), {tag(i - 1, j), tag(i, j)});
      }
    }
  }

  // Create 4 discrete points for the 4 corners of the terrain surface:
  for(int i = 0; i < 4; i++) model::addDiscreteEntity(0, i + 1);
  model::setCoordinates(1, 0, 0, coords[3 * tag(0, 0) - 1]);
  model::setCoordinates(2, 1, 0, coords[3 * tag(N, 0) - 1]);
  model::setCoordinates(3, 1, 1, coords[3 * tag(N, N) - 1]);
  model::setCoordinates(4, 0, 1, coords[3 * tag(0, N) - 1]);

  // Create 4 discrete bounding curves, with their boundary points:
  for(int i = 0; i < 4; i++)
    model::addDiscreteEntity(1, i + 1, {i + 1, (i < 3) ? (i + 2) : 1});

  // Create one discrete surface, with its bounding curves:
  model::addDiscreteEntity(2, 1, {1, 2, -3, -4});

  // Add all the nodes on the surface (for simplicity... see below):
  model::mesh::addNodes(2, 1, nodes, coords);

  // Add point elements on the 4 points, line elements on the 4 curves, and
  // triangle elements on the surface:
  for(int i = 0; i < 4; i++) {
    // Type 15 for point elements:
    model::mesh::addElementsByType(i + 1, 15, {}, {pnt[i]});
    // Type 1 for 2-node line elements:
    model::mesh::addElementsByType(i + 1, 1, {}, lin[i]);
  }
  // Type 2 for 3-node triangle elements:
  model::mesh::addElementsByType(1, 2, {}, tris);

  // Reclassify the nodes on the curves and the points (since we put them all on
  // the surface before with `addNodes' for simplicity)
  model::mesh::reclassifyNodes();

  // Create a geometry for the discrete curves and surfaces, so that we can
  // remesh them later on:
  model::mesh::createGeometry();

  // Note that for more complicated meshes, e.g. for on input unstructured STL
  // mesh, we could use `classifySurfaces()' to automatically create the
  // discrete entities and the topology; but we would then have to extract the
  // boundaries afterwards.

  // Create other CAD entities to form one volume below the terrain surface.
  // Beware that only built-in CAD entities can be hybrid, i.e. have discrete
  // entities on their boundary: OpenCASCADE does not support this feature.
  int p1 = model::geo::addPoint(0, 0, -0.5);
  int p2 = model::geo::addPoint(1, 0, -0.5);
  int p3 = model::geo::addPoint(1, 1, -0.5);
  int p4 = model::geo::addPoint(0, 1, -0.5);
  int c1 = model::geo::addLine(p1, p2);
  int c2 = model::geo::addLine(p2, p3);
  int c3 = model::geo::addLine(p3, p4);
  int c4 = model::geo::addLine(p4, p1);
  int c10 = model::geo::addLine(p1, 1);
  int c11 = model::geo::addLine(p2, 2);
  int c12 = model::geo::addLine(p3, 3);
  int c13 = model::geo::addLine(p4, 4);
  int ll1 = model::geo::addCurveLoop({c1, c2, c3, c4});
  int s1 = model::geo::addPlaneSurface({ll1});
  int ll3 = model::geo::addCurveLoop({c1, c11, -1, -c10});
  int s3 = model::geo::addPlaneSurface({ll3});
  int ll4 = model::geo::addCurveLoop({c2, c12, -2, -c11});
  int s4 = model::geo::addPlaneSurface({ll4});
  int ll5 = model::geo::addCurveLoop({c3, c13, 3, -c12});
  int s5 = model::geo::addPlaneSurface({ll5});
  int ll6 = model::geo::addCurveLoop({c4, c10, 4, -c13});
  int s6 = model::geo::addPlaneSurface({ll6});
  int sl1 = model::geo::addSurfaceLoop({s1, s3, s4, s5, s6, 1});
  int v1 = model::geo::addVolume({sl1});
  model::geo::synchronize();

  // Set this to true to build a fully hex mesh:
  bool transfinite = false;
  bool transfiniteAuto = false;

  if(transfinite) {
    int NN = 30;
    vector<pair<int, int> > tmp;
    model::getEntities(tmp, 1);
    for(auto c : tmp) { model::mesh::setTransfiniteCurve(c.second, NN); }
    model::getEntities(tmp, 2);
    for(auto s : tmp) {
      model::mesh::setTransfiniteSurface(s.second);
      model::mesh::setRecombine(s.first, s.second);
      model::mesh::setSmoothing(s.first, s.second, 100);
    }
    model::mesh::setTransfiniteVolume(v1);
  }
  else if(transfiniteAuto) {
    option::setNumber("Mesh.MeshSizeMin", 0.5);
    option::setNumber("Mesh.MeshSizeMax", 0.5);
    // setTransfiniteAutomatic() uses the sizing constraints to set the number
    // of points
    model::mesh::setTransfiniteAutomatic();
  }
  else {
    option::setNumber("Mesh.MeshSizeMin", 0.05);
    option::setNumber("Mesh.MeshSizeMax", 0.05);
  }

  model::mesh::generate(3);
  write("x2.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 3
//
//  Post-processing data import: list-based
//
// -----------------------------------------------------------------------------
void x3(int argc, char **argv)
{
    // Gmsh supports two types of post-processing data: "list-based" and
  // "model-based". Both types of data are handled through the `view' interface.

  // List-based views are completely independent from any model and any mesh:
  // they are self-contained and simply contain lists of coordinates and values,
  // element by element, for 3 types of fields (scalar "S", vector "V" and
  // tensor "T") and several types of element shapes (point "P", line "L",
  // triangle "T", quadrangle "Q", tetrahedron "S", hexahedron "H", prism "I"
  // and pyramid "Y"). (See `x4.cpp' for a tutorial on model-based views.)

  // To create a list-based view one should first create a view:
  int t1 = view::add("A list-based view");

  // List-based data is then added by specifying the type as a 2 character
  // string that combines a field type and an element shape (e.g. "ST" for a
  // scalar field on triangles), the number of elements to be added, and the
  // concatenated list of coordinates (e.g. 3 "x" coordinates, 3 "y"
  // coordinates, 3 "z" coordinates for first order triangles) and values for
  // each element (e.g. 3 values for first order scalar triangles, repeated for
  // each step if there are several time steps).

  // Let's create two triangles...
  vector<double> triangle1 = {
    0., 1., 1., // x coordinates of the 3 triangle nodes
    0., 0., 1., // y coordinates of the 3 triangle nodes
    0., 0., 0.}; // z coordinates of the 3 triangle nodes
  vector<double> triangle2 = {0., 1., 0., 0., 1., 1., 0., 0., 0.};

  // ... and append values for 10 time steps
  for(int step = 0; step < 10; step++) {
    // 3 node values for each step
    triangle1.insert(triangle1.end(), {10., 11. - step, 12.});
    triangle2.insert(triangle2.end(), {11., 12., 13. + step});
  }

  // List-based data is just added by concatenating the data for all the
  // triangles:
  vector<double> triangles(triangle1);
  triangles.insert(triangles.end(), triangle2.begin(), triangle2.end());
  view::addListData(t1, "ST", 2, triangles);

  // Internally, post-processing views parsed by the .geo file parser create
  // such list-based data (see e.g. `t7.cpp', `t8.cpp' and `t9.cpp'),
  // independently of any mesh.

  // Vector or tensor fields can be imported in the same way, the only
  // difference beeing the type (starting with "V" for vector fields and "T" for
  // tensor fields) and the number of components. For example a vector field on
  // a line element can be added as follows:
  vector<double> line = {0.,  1., // x coordinate of the 2 line nodes
                              1.2, 1.2, // y coordinate of the 2 line nodes
                              0.,  0.}; // z coordinate of the 2 line nodes

  for(int step = 0; step < 10; step++) {
    // 3 vector components for each node (2 nodes here), for each step
    line.insert(line.end(), {10. + step, 0., 0., 10. + step, 0., 0.});
  }
  view::addListData(t1, "VL", 1, line);

  // List-based data can also hold 2D (in window coordinates) and 3D (in model
  // coordinates) strings (see `t4.cpp'). Here we add a 2D string located on the
  // bottom-left of the window (with a 20 pixels offset), as well as a 3D string
  // located at model coordinates (0.5, 0.5. 0):
  view::addListDataString(t1, {20., -20.}, {"Created with Gmsh"});
  view::addListDataString(t1, {0.5, 1.5, 0.},
                                {"A multi-step list-based view"},
                                {"Align", "Center", "Font", "Helvetica"});

  // The various attributes of the view can be queried and changed using the
  // option interface. Beware that the option interface uses view indices
  // instead of view tags; so to change the current time step and the intervals
  // type, and to retrieve the total number of steps, one would do:
  string v1 = "View[" + to_string(view::getIndex(t1)) + "]";
  option::setNumber(v1 + ".TimeStep", 5);
  option::setNumber(v1 + ".IntervalsType", 3);
  double ns;
  option::getNumber(v1 + ".NbTimeStep", ns);
  cout << v1 << " with tag " << t1 << " has " << ns << " time steps\n";

  // Views can be queried and modified in various ways using plugins (see
  // `t9.py'), or probed directly using `view::probe()' - here at point
  // (0.9, 0.1, 0):
  vector<double> val;
  double distance;
  view::probe(t1, 0.9, 0.1, 0, val, distance);
  cout << "Value at (0.9, 0.1, 0):";
  for(auto v : val) cout << " " << v;
  cout << "\n";

  // Views can be saved to disk using `view::write()':
  view::write(t1, "x3.pos");

  // High-order datasets can be provided by setting the interpolation matrices
  // explicitly. Let's create a second view with second order interpolation on a
  // 4-node quadrangle.

  // Add a new view:
  int t2 = view::add("Second order quad");

  // Set the node coordinates:
  vector<double> quad = {
    0.,   1.,   1.,   0., // x coordinates of the 4 quadrangle nodes
    -1.2, -1.2, -0.2, -0.2, // y coordinates of the 4 quadrangle nodes
    0.,   0.,   0.,   0.}; // z coordinates of the 4 quadrangle nodes

  // Add nine values that will be interpolated by second order basis functions
  quad.insert(quad.end(), {1., 1., 1., 1., 3., 3., 3., 3., -3.});

  // Set the two interpolation matrices c[i][j] and e[i][j] defining the d = 9
  // basis functions: f[i](u, v, w) = sum_(j = 0, ..., d - 1) c[i][j] u^e[j][0]
  // v^e[j][1] w^e[j][2], i = 0, ..., d-1, with u, v, w the coordinates in the
  // reference element:
  view::setInterpolationMatrices(
    t2, "Quadrangle", 9,
    {0,    0,   0.25,  0,    0,   -0.25, -0.25, 0,    0.25, 0,     0,   0.25,
     0,    0,   -0.25, 0.25, 0,   -0.25, 0,     0,    0.25, 0,     0,   0.25,
     0.25, 0,   0.25,  0,    0,   0.25,  0,     0,    0.25, -0.25, 0,   -0.25,
     0,    0,   -0.5,  0.5,  0,   0.5,   0,     -0.5, 0,    0,     0.5, -0.5,
     0,    0.5, 0,     -0.5, 0,   0,     0,     0,    -0.5, 0.5,   0,   -0.5,
     0,    0.5, 0,     0,    0.5, -0.5,  0,     -0.5, 0,    0.5,   0,   0,
     1,    -1,  1,     -1,   0,   0,     0,     0,    0},
    {0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 1, 0,
     0, 2, 1, 0, 1, 2, 0, 0, 1, 0, 1, 1, 0});

  // Note that two additional interpolation matrices could also be provided to
  // interpolate the geometry, i.e. to interpolate curved elements.

  // Add the data to the view:
  view::addListData(t2, "SQ", 1, quad);

  // In order to visualize the high-order field, one must activate adaptive
  // visualization, set a visualization error threshold and a maximum
  // subdivision level (Gmsh does automatic mesh refinement to visualize the
  // high-order field with the requested accuracy):
  string v2 = "View[" + to_string(view::getIndex(t2)) + "]";
  option::setNumber(v2 + ".AdaptVisualizationGrid", 1);
  option::setNumber(v2 + ".TargetError", 1e-2);
  option::setNumber(v2 + ".MaxRecursionLevel", 5);

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 4
//
//  Post-processing data import: model-based
//
// -----------------------------------------------------------------------------
void x4(int argc, char **argv)
{
    // Contrary to list-based view (see `x3.cpp'), model-based views are based on
  // one or more meshes. Compared to list-based views, they are thus linked to
  // one model (per step). Post-processing data stored in MSH files create such
  // model-based views.

  // Let's create a first model-based view using a simple mesh contructed by
  // hand. We create a model with a discrete surface
  model::add("simple model");
  int surf = model::addDiscreteEntity(2);

  // We add 4 nodes and 2 3-node triangles (element type "2")
  model::mesh::addNodes(2, surf, {1, 2, 3, 4},
                              {0., 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 0.});
  model::mesh::addElementsByType(surf, 2, {1, 2}, {1, 2, 3, 1, 3, 4});

  // We can now create a new model-based view, to which we add 10 steps of
  // node-based data:
  int t1 = view::add("A model-based view");
  for(int step = 0; step < 10; step++) {
    view::addHomogeneousModelData(
      t1, step, "simple model", "NodeData", {1, 2, 3, 4}, // tags of nodes
      {10., 10., 12. + step, 13. + step}); // data, per node
  }

  // Besided node-based data, which result in continuous fields, one can also
  // add general discontinous fields defined at the nodes of each element, using
  // "ElementNodeData":
  int t2 = view::add("A discontinuous model-based view");
  for(int step = 0; step < 10; step++) {
    view::addHomogeneousModelData(
      t2, step, "simple model", "ElementNodeData", {1, 2}, // tags of elements
      {10., 10., 12. + step, 14., 15., 13. + step}); // data per element nodes
  }

  // Constant per element datasets can also be created using "ElementData". Note
  // that a more general function `addModelData' to add data for hybrid meshes
  // (when data is not homogeneous, i.e. when the number of nodes changes
  // between elements) is also available.

  // Each step of a model-based view can be defined on a different model,
  // i.e. on a different mesh. Let's define a second model and mesh it
  model::add("another model");
  model::occ::addBox(0, 0, 0, 1, 1, 1);
  model::occ::synchronize();
  model::mesh::generate(3);

  // We can add other steps to view "t" based on this new mesh:
  vector<size_t> nodes;
  vector<double> coord, coordParam;
  model::mesh::getNodes(nodes, coord, coordParam);
  for(int step = 11; step < 20; step++) {
    vector<double> val;
    for(size_t i = 0; i < coord.size(); i += 3)
      val.push_back(step * coord[i]);
    view::addHomogeneousModelData(t1, step, "another model", "NodeData",
                                        nodes, val);
  }

  // This feature allows to create seamless animations for time-dependent
  // datasets on deforming or remeshed models.

  // High-order node-based datasets are supported without needing to supply the
  // interpolation matrices (iso-parametric Lagrange elements). Arbitrary
  // high-order datasets can be specified as "ElementNodeData", with the
  // interpolation matrices specified in the same as as for list-based views
  // (see `x3.cpp').

  // Model-based views can be saved to disk using `view::write()'; note
  // that saving a view based on multiple meshes (like the view `t1') will
  // automatically create several files. If the `PostProcessing.SaveMesh' option
  // is not set, `view::write()' will only save the view data, without the
  // mesh (which could be saved independently with `write()').
  view::write(t1, "x4_t1.msh");
  view::write(t2, "x4_t2.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 5
//
//  Additional geometrical data: parametrizations, normals, curvatures
//
// -----------------------------------------------------------------------------
void x5(int argc, char **argv)
{
    // The API provides access to geometrical data in a CAD kernel agnostic
  // manner.

  // Let's create a simple CAD model by fusing a sphere and a cube, then mesh
  // the surfaces:
  model::add("x5");
  int s = model::occ::addSphere(0, 0, 0, 1);
  int b = model::occ::addBox(0.5, 0, 0, 1.3, 2, 3);
  vector<pair<int, int> > ov;
  vector<vector<pair<int, int> > > ovv;
  model::occ::fuse({{3, s}}, {{3, b}}, ov, ovv);
  model::occ::synchronize();
  model::mesh::generate(2);

  // We can for example retrieve the exact normals and the curvature at all the
  // mesh nodes (i.e. not normals and curvatures computed from the mesh, but
  // directly evaluated on the geometry), by querying the CAD kernels at the
  // corresponding parametric coordinates.
  vector<double> normals, curvatures;

  // For each surface in the model:
  vector<pair<int, int> > entities;
  model::getEntities(entities, 2);
  for(auto e : entities) {
    // Retrieve the surface tag
    int s = e.second;

    // Get the mesh nodes on the surface, including those on the boundary
    // (contrary to internal nodes, which store their parametric coordinates,
    // boundary nodes will be reparametrized on the surface in order to compute
    // their parametric coordinates, the result being different when
    // reparametrized on another adjacent surface)
    vector<size_t> tags;
    vector<double> coord, param;
    model::mesh::getNodes(tags, coord, param, 2, s, true);

    // Get the surface normals on all the points on the surface corresponding to
    // the parametric coordinates of the nodes
    vector<double> norm;
    model::getNormal(s, param, norm);

    // In the same way, get the curvature
    vector<double> curv;
    model::getCurvature(2, s, param, curv);

    // Store the normals and the curvatures so that we can display them as
    // list-based post-processing views
    for(size_t i = 0; i < coord.size(); i += 3) {
      normals.push_back(coord[i]);
      normals.push_back(coord[i + 1]);
      normals.push_back(coord[i + 2]);
      normals.push_back(norm[i]);
      normals.push_back(norm[i + 1]);
      normals.push_back(norm[i + 2]);
      curvatures.push_back(coord[i]);
      curvatures.push_back(coord[i + 1]);
      curvatures.push_back(coord[i + 2]);
      curvatures.push_back(curv[i / 3]);
    }
  }

  // Create a list-based vector view on points to display the normals, and a
  // scalar view on points to display the curvatures
  int vn = view::add("normals");
  view::addListData(vn, "VP", normals.size() / 6, normals);
  int vc = view::add("curvatures");
  view::addListData(vc, "SP", curvatures.size() / 4, curvatures);

  // We can also retrieve the parametrization bounds of model entities, e.g. of
  // curve 5, and evaluate the parametrization for several parameter values:
  vector<double> bounds[2], t, xyz1;
  model::getParametrizationBounds(1, 5, bounds[0], bounds[1]);
  int N = 20;
  for(int i = 0; i < N; i++)
    t.push_back(bounds[0][0] + i * (bounds[1][0] - bounds[0][0]) / N);
  model::getValue(1, 5, t, xyz1);

  // We can also reparametrize curve 5 on surface 1, and evaluate the points in
  // the parametric plane of the surface:
  vector<double> uv, xyz2;
  model::reparametrizeOnSurface(1, 5, t, 1, uv);
  model::getValue(2, 1, uv, xyz2);

  // Hopefully we get the same x, y, z coordinates!
  vector<double> diff(xyz1.size());
  transform(xyz1.begin(), xyz1.end(), xyz2.begin(), diff.begin(),
                 [](double a, double b){ return abs(a - b); });
  if(*max_element(diff.begin(), diff.end()) < 1e-12)
    logger::write("Evaluation on curve and surface match!");
  else
    logger::write("Evaluation on curve and surface do not match!",
                        "error");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 6
//
//  Additional mesh data: integration points, Jacobians and basis functions
//
// -----------------------------------------------------------------------------
void x6(int argc, char **argv)
{
    model::add("x6");

  // The API provides access to all the elementary building blocks required to
  // implement finite-element-type numerical methods. Let's create a simple 2D
  // model and mesh it:
  model::occ::addRectangle(0, 0, 0, 1, 0.1);
  model::occ::synchronize();
  model::mesh::setTransfiniteAutomatic();
  model::mesh::generate(2);

  // Set the element order and the desired interpolation order:
  int elementOrder = 1, interpolationOrder = 2;
  model::mesh::setOrder(elementOrder);

  auto pp = [](const string &label, const vector<double> &v, int mult)
  {
    cout << " * " << v.size() / mult << " " << label << ": ";
    for(auto c : v) cout << c << " ";
    cout << "\n";
  };

  // Iterate over all the element types present in the mesh:
  vector<int> elementTypes;
  model::mesh::getElementTypes(elementTypes);
  for(auto t : elementTypes) {
    // Retrieve properties for the given element type
    string elementName;
    int dim, order, numNodes, numPrimNodes;
    vector<double> localNodeCoord;
    model::mesh::getElementProperties
      (t, elementName, dim, order, numNodes, localNodeCoord, numPrimNodes);
    cout << "\n** " << elementName << " **\n\n";

    // Retrieve integration points for that element type, enabling exact
    // integration of polynomials of order "interpolationOrder". The "Gauss"
    // integration family returns the "economical" Gauss points if available,
    // and defaults to the "CompositeGauss" (tensor product) rule if not.
    vector<double> localCoords, weights;
    model::mesh::getIntegrationPoints
      (t, "Gauss" + to_string(interpolationOrder), localCoords, weights);
    pp("integration points to integrate order " +
       to_string(interpolationOrder) + " polynomials", localCoords, 3);

    // Return the basis functions evaluated at the integration points. Selecting
    // "Lagrange" and "GradLagrange" returns the isoparamtric basis functions
    // and their gradient (in the reference space of the given element type). A
    // specific interpolation order can be requested using "LagrangeN" and
    // "GradLagrangeN" with N = 1, 2, ... Other supported function spaces
    // include "H1LegendreN", "GradH1LegendreN", "HcurlLegendreN",
    // "CurlHcurlLegendreN".
    int numComponents, numOrientations;
    vector<double> basisFunctions;
    model::mesh::getBasisFunctions
      (t, localCoords, "Lagrange", numComponents, basisFunctions,
       numOrientations);
    pp("basis functions at integration points", basisFunctions, 1);
    model::mesh::getBasisFunctions
      (t, localCoords, "GradLagrange", numComponents, basisFunctions,
       numOrientations);
    pp("basis function gradients at integration points", basisFunctions, 3);

    // Compute the Jacobians (and their determinants) at the integration points
    // for all the elements of the given type in the mesh. Beware that the
    // Jacobians are returned "by column": see the API documentation for
    // details.
    vector<double> jacobians, determinants, coords;
    model::mesh::getJacobians
      (t, localCoords, jacobians, determinants, coords);
    pp("Jacobian determinants at integration points", determinants, 1);
  }
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ extended tutorial 7
//
//  Additional mesh data: internal edges and faces
//
// -----------------------------------------------------------------------------
void x7(int argc, char **argv)
{
    model::add("x7");

  // Meshes are fully described in Gmsh by nodes and elements, both associated
  // to model entities. The API can be used to generate and handle other mesh
  // entities, i.e. mesh edges and faces, which are not stored by default.

  // Let's create a simple model and mesh it:
  model::occ::addBox(0, 0, 0, 1, 1, 1);
  model::occ::synchronize();
  option::setNumber("Mesh.MeshSizeMin", 2.);
  model::mesh::generate(3);

  // Like elements, mesh edges and faces are described by (an ordered list of)
  // their nodes. Let us retrieve the edges and the (triangular) faces of all
  // the first order tetrahedra in the mesh:
  int elementType = model::mesh::getElementType("tetrahedron", 1);
  vector<size_t> edgeNodes, faceNodes;
  model::mesh::getElementEdgeNodes(elementType, edgeNodes);
  model::mesh::getElementFaceNodes(elementType, 3, faceNodes);

  // Edges and faces are returned for each element as a list of nodes
  // corresponding to the canonical orientation of the edges and faces for a
  // given element type.

  // Gmsh can also identify unique edges and faces (a single edge or face
  // whatever the ordering of their nodes) and assign them a unique tag. This
  // identification can be done internally by Gmsh (e.g. when generating keys
  // for basis functions), or requested explicitly as follows:
  model::mesh::createEdges();
  model::mesh::createFaces();

  // Edge and face tags can then be retrieved by providing their nodes:
  vector<size_t> edgeTags, faceTags;
  vector<int> edgeOrientations, faceOrientations;
  model::mesh::getEdges(edgeNodes, edgeTags, edgeOrientations);
  model::mesh::getFaces(3, faceNodes, faceTags, faceOrientations);

  // Since element edge and face nodes are returned in the same order as the
  // elements, one can easily keep track of which element(s) each edge or face
  // is connected to:
  vector<size_t> elementTags, elementNodeTags;
  model::mesh::getElementsByType(elementType, elementTags, elementNodeTags);
  map<size_t, vector<size_t>> edges2Elements, faces2Elements;
  for(size_t i = 0; i < edgeTags.size(); i++) // 6 edges per tetrahedron
    edges2Elements[edgeTags[i]].push_back(elementTags[i / 6]);
  for(size_t i = 0; i < faceTags.size(); i++) // 4 faces per tetrahedron
    faces2Elements[faceTags[i]].push_back(elementTags[i / 4]);

  // New unique lower dimensional elements can also be easily created given the
  // edge or face nodes. This is especially useful for numerical methods that
  // require integrating or interpolating on internal edges or faces (like
  // e.g. Discontinuous Galerkin techniques), since creating elements for the
  // internal entities will make this additional mesh data readily available
  // (see `x6.cpp'). For example, we can create a new discrete surface...
  int s = model::addDiscreteEntity(2);

  // ... and fill it with unique triangles corresponding to the faces of the
  // tetrahedra:
  set<size_t> uniqueFaceTags;
  vector<size_t> tagsForTriangles, faceNodesForTriangles;
  size_t maxElementTag;
  model::mesh::getMaxElementTag(maxElementTag);
  for(size_t i = 0; i < faceTags.size(); i++) {
    if(uniqueFaceTags.find(faceTags[i]) == uniqueFaceTags.end()) {
      uniqueFaceTags.insert(faceTags[i]);
      tagsForTriangles.push_back(faceTags[i] + maxElementTag);
      faceNodesForTriangles.push_back(faceNodes[3 * i]);
      faceNodesForTriangles.push_back(faceNodes[3 * i + 1]);
      faceNodesForTriangles.push_back(faceNodes[3 * i + 2]);
    }
  }
  int elementType2D = model::mesh::getElementType("triangle", 1);
  model::mesh::addElementsByType(s, elementType2D, tagsForTriangles,
                                       faceNodesForTriangles);

  // Since the tags for the triangles have been created based on the face tags,
  // the information about neighboring elements can also be readily created,
  // useful e.g. in Finite Volume or Discontinuous Galerkin techniques:
  for(auto t : tagsForTriangles) {
    cout << "triangle " << t << " is connected to tetrahedra ";
    for(auto tt : faces2Elements[t - maxElementTag]) cout << tt << " ";
    cout << "\n";
  }

  // If all you need is the list of all edges or faces in terms of their nodes, you
  // can also directly call:
  model::mesh::getAllEdges(edgeTags, edgeNodes);
  model::mesh::getAllFaces(3, faceTags, faceNodes);

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 1
//
//  Geometry basics, elementary entities, physical groups
//
// -----------------------------------------------------------------------------
void t1(int argc, char **argv)
{
      // We now add a new model, named "t1". If model::add() is not called, a
  // new default (unnamed) model will be created on the fly, if necessary.
  model::add("t1");

  // The C++ API provides direct access to each supported geometry (CAD)
  // kernel. The built-in kernel is used in this first tutorial: the
  // corresponding API functions live in the "model::geo" namespace.

  // The first type of `elementary entity' in Gmsh is a `Point'. To create a
  // point with the built-in CAD kernel, the C++ API function is
  // model::geo::addPoint():
  // - the first 3 arguments are the point coordinates (x, y, z)
  // - the next (optional) argument is the target mesh size close to the point
  // - the last (optional) argument is the point tag (a stricly positive integer
  //   that uniquely identifies the point)
  double lc = 1e-2;
  model::geo::addPoint(0, 0, 0, lc, 1);

  // The distribution of the mesh element sizes will be obtained by
  // interpolation of these mesh sizes throughout the geometry. Another method
  // to specify mesh sizes is to use general mesh size Fields (see `t10.cpp'). A
  // particular case is the use of a background mesh (see `t7.cpp').
  //
  // If no target mesh size of provided, a default uniform coarse size will be
  // used for the model, based on the overall model size.
  //
  // We can then define some additional points. All points should have different
  // tags:
  model::geo::addPoint(.1, 0, 0, lc, 2);
  model::geo::addPoint(.1, .3, 0, lc, 3);

  // If the tag is not provided explicitly, a new tag is automatically created,
  // and returned by the function:
  int p4 = model::geo::addPoint(0, .3, 0, lc);

  // Curves are Gmsh's second type of elementery entities, and, amongst curves,
  // straight lines are the simplest. The API to create straight line segments
  // with the built-in kernel follows the same conventions: the first 2
  // arguments are point tags (the start and end points of the line), and the
  // last (optional one) is the line tag.
  //
  // In the commands below, for example, the line 1 starts at point 1 and ends
  // at point 2.
  //
  // Note that curve tags are separate from point tags - hence we can reuse tag
  // `1' for our first curve. And as a general rule, elementary entity tags in
  // Gmsh have to be unique per geometrical dimension.
  model::geo::addLine(1, 2, 1);
  model::geo::addLine(3, 2, 2);
  model::geo::addLine(3, p4, 3);
  model::geo::addLine(4, 1, p4);

  // The third elementary entity is the surface. In order to define a simple
  // rectangular surface from the four curves defined above, a curve loop has
  // first to be defined. A curve loop is defined by an ordered list of
  // connected curves, a sign being associated with each curve (depending on the
  // orientation of the curve to form a loop). The API function to create curve
  // loops takes a vector of integers as first argument, and the curve loop tag
  // (which must be unique amongst curve loops) as the second (optional)
  // argument:
  model::geo::addCurveLoop({4, 1, -2, 3}, 1);

  // We can then define the surface as a list of curve loops (only one here,
  // representing the external contour, since there are no holes--see `t4.cpp'
  // for an example of a surface with a hole):
  model::geo::addPlaneSurface({1}, 1);

  // Before they can be meshed (and, more generally, before they can be used by
  // API functions outside of the built-in CAD kernel functions), the CAD
  // entities must be synchronized with the Gmsh model, which will create the
  // relevant Gmsh data structures. This is achieved by the
  // gmsh.model.geo.synchronize() API call for the built-in CAD
  // kernel. Synchronizations can be called at any time, but they involve a non
  // trivial amount of processing; so while you could synchronize the internal
  // CAD data after every CAD command, it is usually better to minimize the
  // number of synchronization points.
  model::geo::synchronize();

  // At this level, Gmsh knows everything to display the rectangular surface 1
  // and to mesh it. An optional step is needed if we want to group elementary
  // geometrical entities into more meaningful groups, e.g. to define some
  // mathematical ("domain", "boundary"), functional ("left wing", "fuselage")
  // or material ("steel", "carbon") properties.
  //
  // Such groups are called "Physical Groups" in Gmsh. By default, if physical
  // groups are defined, Gmsh will export in output files only mesh elements
  // that belong to at least one physical group. (To force Gmsh to save all
  // elements, whether they belong to physical groups or not, set the
  // `Mesh.SaveAll' option to 1.) Physical groups are also identified by tags,
  // i.e. stricly positive integers, that should be unique per dimension (0D,
  // 1D, 2D or 3D). Physical groups can also be given names.
  //
  // Here we define a physical curve that groups the left, bottom and right
  // curves in a single group (with prescribed tag 5); and a physical surface
  // with name "My surface" (with an automatic tag) containing the geometrical
  // surface 1:
  model::addPhysicalGroup(1, {1, 2, 4}, 5);
  int ps = model::addPhysicalGroup(2, {1});
  model::setPhysicalName(2, ps, "My surface");

  // We can then generate a 2D mesh...
  model::mesh::generate(2);

  // ... and save it to disk
  write("t1.msh");

  // Remember that by default, if physical groups are defined, Gmsh will export
  // in the output mesh file only those elements that belong to at least one
  // physical group. To force Gmsh to save all elements, you can use
  //
  // option::setNumber("Mesh.SaveAll", 1);

  // By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
  // format (the `MSH' format). You can save meshes in other mesh formats by
  // specifying a filename with a different extension. For example
  //
  //   write("t1.unv");
  //
  // will save the mesh in the UNV format. You can also save the mesh in older
  // versions of the MSH format: simply set
  //
  //   option::setNumber("Mesh.MshFileVersion", x);
  //
  // for any version number `x'. As an alternative, you can also not specify
  // the format explicitly, and just choose a filename with the `.msh2' or
  // `.msh4' extension.

  // To visualize the model we can run the graphical user interface with
  // `fltk::run()'. Here we run it only if "-nopopup" is not provided in
  // the command line arguments:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();

  // Note that starting with Gmsh 3.0, models can be built using other geometry
  // kernels than the default "built-in" kernel. To use the OpenCASCADE CAD
  // kernel instead of the built-in kernel, you should use the functions in the
  // "model::occ" namespace.
  //
  // Different CAD kernels have different features. With OpenCASCADE, instead of
  // defining the surface by successively defining 4 points, 4 curves and 1
  // curve loop, one can define the rectangular surface directly with
  //
  // model::occ::addRectangle(.2, 0, 0, .1, .3);
  //
  // After synchronization with the Gmsh model with
  //
  // model::occ::synchronize();
  //
  // the underlying curves and points could be accessed with
  // model::getBoundary().
  //
  // See e.g `t16.cpp', `t18.cpp', `t19.cpp' or `t20.cpp' for complete examples
  // based on OpenCASCADE, and `examples/api' for more.

  // This should be called when you are done using the Gmsh C++ API:
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 2
//
//  Transformations, extruded geometries, volumes
//
// -----------------------------------------------------------------------------
void t2(int argc, char **argv)
{
    model::add("t2");

  // Copied from t1.cpp...
  double lc = 1e-2;
  model::geo::addPoint(0, 0, 0, lc, 1);
  model::geo::addPoint(.1, 0, 0, lc, 2);
  model::geo::addPoint(.1, .3, 0, lc, 3);
  model::geo::addPoint(0, .3, 0, lc, 4);
  model::geo::addLine(1, 2, 1);
  model::geo::addLine(3, 2, 2);
  model::geo::addLine(3, 4, 3);
  model::geo::addLine(4, 1, 4);
  model::geo::addCurveLoop({4, 1, -2, 3}, 1);
  model::geo::addPlaneSurface({1}, 1);
  model::geo::synchronize();
  model::addPhysicalGroup(1, {1, 2, 4}, 5);
  int ps = model::addPhysicalGroup(2, {1});
  model::setPhysicalName(2, ps, "My surface");

  // We can then add new points and curves in the same way as we did in
  // `t1.cpp':
  model::geo::addPoint(0, .4, 0, lc, 5);
  model::geo::addLine(4, 5, 5);

  // But Gmsh also provides tools to transform (translate, rotate, etc.)
  // elementary entities or copies of elementary entities.  Geometrical
  // transformations take a vector of pairs of integers as first argument, which
  // contains the list of entities, represented by (dimension, tag) pairs.  For
  // example, the point 5 (dimension=0, tag=5) can be moved by 0.02 to the left
  // (dx=-0.02, dy=0, dz=0) with
  model::geo::translate({{0, 5}}, -0.02, 0, 0);

  // And it can be further rotated by -Pi/4 around (0, 0.3, 0) (with the
  // rotation along the z axis) with:
  model::geo::rotate({{0, 5}}, 0, 0.3, 0, 0, 0, 1, -M_PI / 4);

  // Note that there are no units in Gmsh: coordinates are just numbers - it's
  // up to the user to associate a meaning to them.

  // Point 3 can be duplicated and translated by 0.05 along the y axis by using
  // the `copy()' function, which takes a vector of (dim, tag) pairs as input,
  // and returns another vector of (dim, tag) pairs:
  vector<pair<int, int> > ov;
  model::geo::copy({{0, 3}}, ov);
  model::geo::translate(ov, 0, 0.05, 0);

  // The new point tag is available in ov[0].second, and can be used to create
  // new lines:
  model::geo::addLine(3, ov[0].second, 7);
  model::geo::addLine(ov[0].second, 5, 8);
  model::geo::addCurveLoop({5, -8, -7, 3}, 10);
  model::geo::addPlaneSurface({10}, 11);

  // In the same way, we can translate copies of the two surfaces 1 and 11 to
  // the right with the following command:
  model::geo::copy({{2, 1}, {2, 11}}, ov);
  model::geo::translate(ov, 0.12, 0, 0);

  printf("New surfaces '%d' and '%d'\n", ov[0].second, ov[1].second);

  // Volumes are the fourth type of elementary entities in Gmsh. In the same way
  // one defines curve loops to build surfaces, one has to define surface loops
  // (i.e. `shells') to build volumes. The following volume does not have holes
  // and thus consists of a single surface loop:
  model::geo::addPoint(0., 0.3, 0.12, lc, 100);
  model::geo::addPoint(0.1, 0.3, 0.12, lc, 101);
  model::geo::addPoint(0.1, 0.35, 0.12, lc, 102);

  // We would like to retrieve the coordinates of point 5 to create point 103,
  // so we synchronize the model, and use `getValue()'
  model::geo::synchronize();
  vector<double> xyz;
  model::getValue(0, 5, {}, xyz);
  model::geo::addPoint(xyz[0], xyz[1], 0.12, lc, 103);

  model::geo::addLine(4, 100, 110);
  model::geo::addLine(3, 101, 111);
  model::geo::addLine(6, 102, 112);
  model::geo::addLine(5, 103, 113);
  model::geo::addLine(103, 100, 114);
  model::geo::addLine(100, 101, 115);
  model::geo::addLine(101, 102, 116);
  model::geo::addLine(102, 103, 117);

  model::geo::addCurveLoop({115, -111, 3, 110}, 118);
  model::geo::addPlaneSurface({118}, 119);
  model::geo::addCurveLoop({111, 116, -112, -7}, 120);
  model::geo::addPlaneSurface({120}, 121);
  model::geo::addCurveLoop({112, 117, -113, -8}, 122);
  model::geo::addPlaneSurface({122}, 123);
  model::geo::addCurveLoop({114, -110, 5, 113}, 124);
  model::geo::addPlaneSurface({124}, 125);
  model::geo::addCurveLoop({115, 116, 117, 114}, 126);
  model::geo::addPlaneSurface({126}, 127);

  model::geo::addSurfaceLoop({127, 119, 121, 123, 125, 11}, 128);
  model::geo::addVolume({128}, 129);

  // When a volume can be extruded from a surface, it is usually easier to use
  // the `extrude()' function directly instead of creating all the points,
  // curves and surfaces by hand. For example, the following command extrudes
  // the surface 11 along the z axis and automatically creates a new volume (as
  // well as all the needed points, curves and surfaces). As expected, the
  // function takes a vector of (dim, tag) pairs as input as well as the
  // translation vector, and returns a vector of (dim, tag) pairs as output:
  vector<pair<int, int> > ov2;
  model::geo::extrude({ov[1]}, 0, 0, 0.12, ov2);

  // Mesh sizes associated to geometrical points can be set by passing a vector
  // of (dim, tag) pairs for the corresponding points:
  model::geo::mesh::setSize(
    {{0, 103}, {0, 105}, {0, 109}, {0, 102}, {0, 28}, {0, 24}, {0, 6}, {0, 5}},
    lc * 3);

  // We finish by synchronizing the data from the built-in CAD kernel with the
  // Gmsh model:
  model::geo::synchronize();

  // We group volumes 129 and 130 in a single physical group with tag `1' and
  // name "The volume":
  model::addPhysicalGroup(3, {129, 130}, 1);
  model::setPhysicalName(3, 1, "The volume");

  // We finally generate and save the mesh:
  model::mesh::generate(3);
  write("t2.msh");

  // Note that, if the transformation tools are handy to create complex
  // geometries, it is also sometimes useful to generate the `flat' geometry,
  // with an explicit representation of all the elementary entities.
  //
  // With the built-in CAD kernel, this can be achieved by saving the model in
  // the `Gmsh Unrolled GEO' format:
  //
  // write("t2.geo_unrolled");
  //
  // With the OpenCASCADE CAD kernel, unrolling the geometry can be achieved by
  // exporting in the `OpenCASCADE BRep' format:
  //
  // write("t2.brep");
  //
  // (OpenCASCADE geometries can also be exported as STEP files.)

  // It is important to note that Gmsh never translates geometry data into a
  // common representation: all the operations on a geometrical entity are
  // performed natively with the associated CAD kernel. Consequently, one cannot
  // export a geometry constructed with the built-in kernel as an OpenCASCADE
  // BRep file; or export an OpenCASCADE model as an Unrolled GEO file.

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 3
//
//  Extruded meshes, ONELAB parameters, options
//
// -----------------------------------------------------------------------------
void t3(int argc, char **argv)
{
    auto createGeometryAndMesh = []()
  {
    // Clear all models and create a new one
    clear();
    model::add("t3");

    // Copied from t1.cpp...
    double lc = 1e-2;
    model::geo::addPoint(0, 0, 0, lc, 1);
    model::geo::addPoint(.1, 0, 0, lc, 2);
    model::geo::addPoint(.1, .3, 0, lc, 3);
    model::geo::addPoint(0, .3, 0, lc, 4);
    model::geo::addLine(1, 2, 1);
    model::geo::addLine(3, 2, 2);
    model::geo::addLine(3, 4, 3);
    model::geo::addLine(4, 1, 4);
    model::geo::addCurveLoop({4, 1, -2, 3}, 1);
    model::geo::addPlaneSurface({1}, 1);
    model::geo::synchronize();
    model::addPhysicalGroup(1, {1, 2, 4}, 5);
    int ps = model::addPhysicalGroup(2, {1});
    model::setPhysicalName(2, ps, "My surface");

    // As in `t2.cpp', we plan to perform an extrusion along the z axis.  But
    // here, instead of only extruding the geometry, we also want to extrude the
    // 2D mesh. This is done with the same `extrude()' function, but by
    // specifying element 'Layers' (2 layers in this case, the first one with 8
    // subdivisions and the second one with 2 subdivisions, both with a height
    // of h/2). The number of elements for each layer and the (end) height of
    // each layer are specified in two vectors:
    double h = 0.1;
    vector<pair<int, int> > ov;
    model::geo::extrude({{2, 1}}, 0, 0, h, ov, {8, 2}, {0.5, 1});

    // The extrusion can also be performed with a rotation instead of a
    // translation, and the resulting mesh can be recombined into prisms (we use
    // only one layer here, with 7 subdivisions). All rotations are specified by
    // an an axis point (-0.1, 0, 0.1), an axis direction (0, 1, 0), and a
    // rotation angle (-Pi/2):
    model::geo::revolve({{2, 28}}, -0.1, 0, 0.1, 0, 1, 0, -M_PI / 2, ov,
                              {7});

    // Using the built-in geometry kernel, only rotations with angles < Pi are
    // supported. To do a full turn, you will thus need to apply at least 3
    // rotations. The OpenCASCADE geometry kernel does not have this limitation.

    // A translation (-2 * h, 0, 0) and a rotation ((0, 0.15, 0.25), (1, 0, 0),
    // angle * Pi / 180) can also be combined to form a "twist". The last
    // (optional) argument for the extrude() and twist() functions specifies
    // whether the extruded mesh should be recombined or not. The `angle'
    // parameter is retrieved from the ONELAB database (it can be set
    // interactively in the GUI -- see below):
    vector<double> angle;
    onelab::getNumber("Parameters/Twisting angle", angle);
    model::geo::twist({{2, 50}}, 0, 0.15, 0.25, -2 * h, 0, 0, 1, 0, 0,
                            angle[0] * M_PI / 180., ov, {10}, {}, true);

    model::geo::synchronize();

    // All the extrusion functions return a vector of extruded entities: the
    // "top" of the extruded surface (in `ov[0]'), the newly created volume (in
    // `ov[1]') and the tags of the lateral surfaces (in `ov[2]', `ov[3]', ...).

    // We can then define a new physical volume (with tag 101) to group all the
    // elementary volumes:
    model::addPhysicalGroup(3, {1, 2, ov[1].second}, 101);

    model::mesh::generate(3);
    write("t3.msh");
  };

  // Let us now change some options... Since all interactive options are
  // accessible through the API, we can for example make point tags visible or
  // redefine some colors:
  option::setNumber("Geometry.PointNumbers", 1);
  option::setColor("Geometry.Color.Points", 255, 165, 0);
  option::setColor("General.Color.Text", 255, 255, 255);
  option::setColor("Mesh.Color.Points", 255, 0, 0);

  // Note that for conciseness "Color." can be ommitted in color options:
  int r, g, b, a;
  option::getColor("Geometry.Points", r, g, b, a);
  option::setColor("Geometry.Surfaces", r, g, b, a);

  // We create a ONELAB parameter to define the angle of the twist. ONELAB
  // parameters can be modified interactively in the GUI, and can be exchanged
  // with other codes connected to the same ONELAB database. The database can be
  // accessed through the Gmsh C++ API using JSON-formatted strings (see
  // https://gitlab.onelab.info/doc/tutorials/-/wikis/ONELAB-JSON-interface for
  // more information):
  onelab::set(R"( [
  {
    "type":"number",
    "name":"Parameters/Twisting angle",
    "values":[90],
    "min":0,
    "max":120,
    "step":1
  }
  ] )");

  // Create the geometry and mesh it:
  createGeometryAndMesh();

  // Launch the GUI and handle the "check" event (recorded in the
  // "ONELAB/Action" parameter) to recreate the geometry with a new twisting
  // angle if necessary:
  auto checkForEvent = [=]() -> bool {
    vector<string> action;
    onelab::getString("ONELAB/Action", action);
    if(action.size() && action[0] == "check") {
      onelab::setString("ONELAB/Action", {""});
      createGeometryAndMesh();
      graphics::draw();
    }
    return true;
  };

  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) {
    fltk::initialize();
    while(fltk::isAvailable() && checkForEvent())
      fltk::wait();
  }

  // When the GUI is launched, you can use the `Help->Current Options and
  // Workspace' menu to see the current values of all options. To save the
  // options in a file, use `File->Export->Gmsh Options', or through the api:

  // write("t3.opt");
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 4
//
//  Holes in surfaces, annotations, entity colors
//
// -----------------------------------------------------------------------------
void t4(int argc, char **argv)
{
    model::add("t4");

  double cm = 1e-02;
  double e1 = 4.5 * cm, e2 = 6 * cm / 2, e3 = 5 * cm / 2;
  double h1 = 5 * cm, h2 = 10 * cm, h3 = 5 * cm, h4 = 2 * cm, h5 = 4.5 * cm;
  double R1 = 1 * cm, R2 = 1.5 * cm, r = 1 * cm;
  double Lc1 = 0.01;
  double Lc2 = 0.003;

  double ccos =
    (-h5 * R1 + e2 * hypot(h5, hypot(e2, R1))) / (h5 * h5 + e2 * e2);
  double ssin = sqrt(1 - ccos * ccos);

  // We start by defining some points and some lines. To make the code shorter
  // we can redefine a namespace:
  namespace factory = model::geo;
  factory::addPoint(-e1 - e2, 0, 0, Lc1, 1);
  factory::addPoint(-e1 - e2, h1, 0, Lc1, 2);
  factory::addPoint(-e3 - r, h1, 0, Lc2, 3);
  factory::addPoint(-e3 - r, h1 + r, 0, Lc2, 4);
  factory::addPoint(-e3, h1 + r, 0, Lc2, 5);
  factory::addPoint(-e3, h1 + h2, 0, Lc1, 6);
  factory::addPoint(e3, h1 + h2, 0, Lc1, 7);
  factory::addPoint(e3, h1 + r, 0, Lc2, 8);
  factory::addPoint(e3 + r, h1 + r, 0, Lc2, 9);
  factory::addPoint(e3 + r, h1, 0, Lc2, 10);
  factory::addPoint(e1 + e2, h1, 0, Lc1, 11);
  factory::addPoint(e1 + e2, 0, 0, Lc1, 12);
  factory::addPoint(e2, 0, 0, Lc1, 13);

  factory::addPoint(R1 / ssin, h5 + R1 * ccos, 0, Lc2, 14);
  factory::addPoint(0, h5, 0, Lc2, 15);
  factory::addPoint(-R1 / ssin, h5 + R1 * ccos, 0, Lc2, 16);
  factory::addPoint(-e2, 0.0, 0, Lc1, 17);

  factory::addPoint(-R2, h1 + h3, 0, Lc2, 18);
  factory::addPoint(-R2, h1 + h3 + h4, 0, Lc2, 19);
  factory::addPoint(0, h1 + h3 + h4, 0, Lc2, 20);
  factory::addPoint(R2, h1 + h3 + h4, 0, Lc2, 21);
  factory::addPoint(R2, h1 + h3, 0, Lc2, 22);
  factory::addPoint(0, h1 + h3, 0, Lc2, 23);

  factory::addPoint(0, h1 + h3 + h4 + R2, 0, Lc2, 24);
  factory::addPoint(0, h1 + h3 - R2, 0, Lc2, 25);

  factory::addLine(1, 17, 1);
  factory::addLine(17, 16, 2);

  // Gmsh provides other curve primitives than straight lines: splines,
  // B-splines, circle arcs, ellipse arcs, etc. Here we define a new circle arc,
  // starting at point 14 and ending at point 16, with the circle's center being
  // the point 15:
  factory::addCircleArc(14, 15, 16, 3);

  // Note that, in Gmsh, circle arcs should always be smaller than Pi. The
  // OpenCASCADE geometry kernel does not have this limitation.

  // We can then define additional lines and circles, as well as a new surface:
  factory::addLine(14, 13, 4);
  factory::addLine(13, 12, 5);
  factory::addLine(12, 11, 6);
  factory::addLine(11, 10, 7);
  factory::addCircleArc(8, 9, 10, 8);
  factory::addLine(8, 7, 9);
  factory::addLine(7, 6, 10);
  factory::addLine(6, 5, 11);
  factory::addCircleArc(3, 4, 5, 12);
  factory::addLine(3, 2, 13);
  factory::addLine(2, 1, 14);
  factory::addLine(18, 19, 15);
  factory::addCircleArc(21, 20, 24, 16);
  factory::addCircleArc(24, 20, 19, 17);
  factory::addCircleArc(18, 23, 25, 18);
  factory::addCircleArc(25, 23, 22, 19);
  factory::addLine(21, 22, 20);

  factory::addCurveLoop({17, -15, 18, 19, -20, 16}, 21);
  factory::addPlaneSurface({21}, 22);

  // But we still need to define the exterior surface. Since this surface has a
  // hole, its definition now requires two curves loops:
  factory::addCurveLoop({11, -12, 13, 14, 1, 2, -3, 4, 5, 6, 7, -8, 9, 10}, 23);
  factory::addPlaneSurface({23, 21}, 24);

  // As a general rule, if a surface has N holes, it is defined by N+1 curve
  // loops: the first loop defines the exterior boundary; the other loops define
  // the boundaries of the holes.

  factory::synchronize();

  // Finally, we can add some comments by creating a post-processing view
  // containing some strings:
  int v = view::add("comments");

  // Add a text string in window coordinates, 10 pixels from the left and 10
  // pixels from the bottom:
  view::addListDataString(v, {10, -10}, {"Created with Gmsh"});

  // Add a text string in model coordinates centered at (X,Y,Z) = (0, 0.11, 0),
  // with some style attributes:
  view::addListDataString(v, {0, 0.11, 0}, {"Hole"},
                                {"Align", "Center", "Font", "Helvetica"});

  // If a string starts with `file://', the rest is interpreted as an image
  // file. For 3D annotations, the size in model coordinates can be specified
  // after a `@' symbol in the form `widthxheight' (if one of `width' or
  // `height' is zero, natural scaling is used; if both are zero, original image
  // dimensions in pixels are used):
  view::addListDataString(
    v, {0, 0.09, 0}, {"file://../t4_image.png@0.01x0"}, {"Align", "Center"});

  // The 3D orientation of the image can be specified by proving the direction
  // of the bottom and left edge of the image in model space:
  view::addListDataString(v, {-0.01, 0.09, 0},
                                {"file://../t4_image.png@0.01x0,0,0,1,0,1,0"});

  // The image can also be drawn in "billboard" mode, i.e. always parallel to
  // the camera, by using the `#' symbol:
  view::addListDataString(
    v, {0, 0.12, 0}, {"file://../t4_image.png@0.01x0#"}, {"Align", "Center"});

  // The size of 2D annotations is given directly in pixels:
  view::addListDataString(v, {150, -7}, {"file://../t4_image.png@20x0"});

  // These annotations are handled by a list-based post-processing view. For
  // large post-processing datasets, that contain actual field values defined on
  // a mesh, you should use model-based post-processing views instead, which
  // allow to efficiently store continuous or discontinuous scalar, vector and
  // tensor fields, or arbitrary polynomial order.

  // Views and geometrical entities can be made to respond to double-click
  // events, here to print some messages to the console:
  option::setString("View[0].DoubleClickedCommand",
                          "Printf('View[0] has been double-clicked!');");
  option::setString("Geometry.DoubleClickedLineCommand",
                          "Printf('Curve %g has been double-clicked!', "
                          "Geometry.DoubleClickedEntityTag);");

  // We can also change the color of some entities:
  model::setColor({{2, 22}}, 127, 127, 127); // Gray50
  model::setColor({{2, 24}}, 160, 32, 240); // Purple
  for(int i = 1; i <= 14; i++)
    model::setColor({{1, i}}, 255, 0, 0); // Red
  for(int i = 15; i <= 20; i++)
    model::setColor({{1, i}}, 255, 255, 0); // Yellow

  model::mesh::generate(2);

  write("t4.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 5
//
//  Mesh sizes, holes in volumes
//
// -----------------------------------------------------------------------------
void t5(int argc, char **argv)
{
  double lcar1 = .1;
  double lcar2 = .0005;
  double lcar3 = .055;

  // If we wanted to change these mesh sizes globally (without changing the
  // above definitions), we could give a global scaling factor for all mesh
  // sizes with e.g.
  //
  // option::setNumber("Mesh.MeshSizeFactor", 0.1);
  //
  // Since we pass `argc' and `argv' to `initialize()', we can also give
  // the option on the command line with the `-clscale' switch. For example,
  // with:
  //
  // > ./t5.exe -clscale 1
  //
  // this tutorial produces a mesh of approximately 3000 nodes and 14,000
  // tetrahedra. With
  //
  // > ./t5.exe -clscale 0.2
  //
  // the mesh counts approximately 231,000 nodes and 1,360,000 tetrahedra. You
  // can check mesh statistics in the graphical user interface
  // (fltk::run()) with the `Tools->Statistics' menu.
  //
  // See `t10.cpp' for more information about mesh sizes.

  // We proceed by defining some elementary entities describing a truncated
  // cube:

  model::geo::addPoint(0.5, 0.5, 0.5, lcar2, 1);
  model::geo::addPoint(0.5, 0.5, 0, lcar1, 2);
  model::geo::addPoint(0, 0.5, 0.5, lcar1, 3);
  model::geo::addPoint(0, 0, 0.5, lcar1, 4);
  model::geo::addPoint(0.5, 0, 0.5, lcar1, 5);
  model::geo::addPoint(0.5, 0, 0, lcar1, 6);
  model::geo::addPoint(0, 0.5, 0, lcar1, 7);
  model::geo::addPoint(0, 1, 0, lcar1, 8);
  model::geo::addPoint(1, 1, 0, lcar1, 9);
  model::geo::addPoint(0, 0, 1, lcar1, 10);
  model::geo::addPoint(0, 1, 1, lcar1, 11);
  model::geo::addPoint(1, 1, 1, lcar1, 12);
  model::geo::addPoint(1, 0, 1, lcar1, 13);
  model::geo::addPoint(1, 0, 0, lcar1, 14);

  model::geo::addLine(8, 9, 1);
  model::geo::addLine(9, 12, 2);
  model::geo::addLine(12, 11, 3);
  model::geo::addLine(11, 8, 4);
  model::geo::addLine(9, 14, 5);
  model::geo::addLine(14, 13, 6);
  model::geo::addLine(13, 12, 7);
  model::geo::addLine(11, 10, 8);
  model::geo::addLine(10, 13, 9);
  model::geo::addLine(10, 4, 10);
  model::geo::addLine(4, 5, 11);
  model::geo::addLine(5, 6, 12);
  model::geo::addLine(6, 2, 13);
  model::geo::addLine(2, 1, 14);
  model::geo::addLine(1, 3, 15);
  model::geo::addLine(3, 7, 16);
  model::geo::addLine(7, 2, 17);
  model::geo::addLine(3, 4, 18);
  model::geo::addLine(5, 1, 19);
  model::geo::addLine(7, 8, 20);
  model::geo::addLine(6, 14, 21);

  model::geo::addCurveLoop({-11, -19, -15, -18}, 22);
  model::geo::addPlaneSurface({22}, 23);
  model::geo::addCurveLoop({16, 17, 14, 15}, 24);
  model::geo::addPlaneSurface({24}, 25);
  model::geo::addCurveLoop({-17, 20, 1, 5, -21, 13}, 26);
  model::geo::addPlaneSurface({26}, 27);
  model::geo::addCurveLoop({-4, -1, -2, -3}, 28);
  model::geo::addPlaneSurface({28}, 29);
  model::geo::addCurveLoop({-7, 2, -5, -6}, 30);
  model::geo::addPlaneSurface({30}, 31);
  model::geo::addCurveLoop({6, -9, 10, 11, 12, 21}, 32);
  model::geo::addPlaneSurface({32}, 33);
  model::geo::addCurveLoop({7, 3, 8, 9}, 34);
  model::geo::addPlaneSurface({34}, 35);
  model::geo::addCurveLoop({-10, 18, -16, -20, 4, -8}, 36);
  model::geo::addPlaneSurface({36}, 37);
  model::geo::addCurveLoop({-14, -13, -12, 19}, 38);
  model::geo::addPlaneSurface({38}, 39);

  vector<int> shells, volumes;

  int sl =
    model::geo::addSurfaceLoop({35, 31, 29, 37, 33, 23, 39, 25, 27});
  shells.push_back(sl);

  // We create five holes in the cube:
  double x = 0, y = 0.75, z = 0, r = 0.09;
  for(int t = 1; t <= 5; t++) {
    x += 0.166;
    z += 0.166;
    cheeseHole(x, y, z, r, lcar3, shells, volumes);
    model::geo::addPhysicalGroup(3, {volumes.back()}, t);
    printf("Hole %d (center = {%g,%g,%g}, radius = %g) has number %d!\n",
                t, x, y, z, r, volumes.back());
  }

  // The volume of the cube, without the 5 holes, is defined by 6 surface loops:
  // the first surface loop defines the exterior surface; the surface loops
  // other than the first one define holes:
  int ve = model::geo::addVolume(shells);

  model::geo::synchronize();

  // Note that using solid modelling with the OpenCASCADE CAD kernel, the same
  // geometry could be built quite differently: see `t16.cpp'.

  // We finally define a physical volume for the elements discretizing the cube,
  // without the holes (for which physical groups were already defined in the
  // `cheeseHole()' function):
  model::addPhysicalGroup(3, {ve}, 10);

  // We could make only part of the model visible to only mesh this subset:
  // vector<pair<int, int> > ent;
  // model::getEntities(ent);
  // model::setVisibility(ent, false);
  // model::setVisibility({{3, 5}}, true, true);
  // option::setNumber("Mesh.MeshOnlyVisible", 1);

  // Meshing algorithms can changed globally using options:
  option::setNumber("Mesh.Algorithm",
                          6); // Frontal-Delaunay for 2D meshes

  // They can also be set for individual surfaces, e.g. for using `MeshAdapt' on
  // surface 1:
  model::mesh::setAlgorithm(2, 33, 1);

  // To generate a curvilinear mesh and optimize it to produce provably valid
  // curved elements (see A. Johnen, J.-F. Remacle and C. Geuzaine. Geometric
  // validity of curvilinear finite elements. Journal of Computational Physics
  // 233, pp. 359-372, 2013; and T. Toulorge, C. Geuzaine, J.-F. Remacle,
  // J. Lambrechts. Robust untangling of curvilinear meshes. Journal of
  // Computational Physics 254, pp. 8-26, 2013), you can uncomment the following
  // lines:
  //
  // option::setNumber("Mesh.ElementOrder", 2);
  // option::setNumber("Mesh.HighOrderOptimize", 2);

  model::mesh::generate(3);
  write("t5.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 6
//
//  Transfinite meshes
//
// -----------------------------------------------------------------------------
void t6(int argc, char **argv)
{
    model::add("t2");

  // Copied from t1.cpp...
  double lc = 1e-2;
  model::geo::addPoint(0, 0, 0, lc, 1);
  model::geo::addPoint(.1, 0, 0, lc, 2);
  model::geo::addPoint(.1, .3, 0, lc, 3);
  model::geo::addPoint(0, .3, 0, lc, 4);
  model::geo::addLine(1, 2, 1);
  model::geo::addLine(3, 2, 2);
  model::geo::addLine(3, 4, 3);
  model::geo::addLine(4, 1, 4);
  model::geo::addCurveLoop({4, 1, -2, 3}, 1);
  model::geo::addPlaneSurface({1}, 1);

  // Delete the surface and the left line, and replace the line with 3 new ones:
  model::geo::remove({{2, 1}, {1, 4}});

  int p1 = model::geo::addPoint(-0.05, 0.05, 0, lc);
  int p2 = model::geo::addPoint(-0.05, 0.1, 0, lc);
  int l1 = model::geo::addLine(1, p1);
  int l2 = model::geo::addLine(p1, p2);
  int l3 = model::geo::addLine(p2, 4);

  // Create surface:
  model::geo::addCurveLoop({2, -1, l1, l2, l3, -3}, 2);
  model::geo::addPlaneSurface({-2}, 1);

  // The `setTransfiniteCurve()' meshing constraints explicitly specifies the
  // location of the nodes on the curve. For example, the following command
  // forces 20 uniformly placed nodes on curve 2 (including the nodes on the two
  // end points):
  model::geo::mesh::setTransfiniteCurve(2, 20);

  // Let's put 20 points total on combination of curves `l1', `l2' and `l3'
  // (beware that the points `p1' and `p2' are shared by the curves, so we do
  // not create 6 + 6 + 10 = 22 nodes, but 20!)
  model::geo::mesh::setTransfiniteCurve(l1, 6);
  model::geo::mesh::setTransfiniteCurve(l2, 6);
  model::geo::mesh::setTransfiniteCurve(l3, 10);

  // Finally, we put 30 nodes following a geometric progression on curve 1
  // (reversed) and on curve 3: Put 30 points following a geometric progression
  model::geo::mesh::setTransfiniteCurve(1, 30, "Progression", -1.2);
  model::geo::mesh::setTransfiniteCurve(3, 30, "Progression", 1.2);

  // The `setTransfiniteSurface()' meshing constraint uses a transfinite
  // interpolation algorithm in the parametric plane of the surface to connect
  // the nodes on the boundary using a structured grid. If the surface has more
  // than 4 corner points, the corners of the transfinite interpolation have to
  // be specified by hand:
  model::geo::mesh::setTransfiniteSurface(1, "Left", {1, 2, 3, 4});

  // To create quadrangles instead of triangles, one can use the `setRecombine'
  // constraint:
  model::geo::mesh::setRecombine(2, 1);

  // When the surface has only 3 or 4 points on its boundary the list of corners
  // can be omitted in the `setTransfiniteSurface()' call:
  model::geo::addPoint(0.2, 0.2, 0, 1.0, 7);
  model::geo::addPoint(0.2, 0.1, 0, 1.0, 8);
  model::geo::addPoint(0, 0.3, 0, 1.0, 9);
  model::geo::addPoint(0.25, 0.2, 0, 1.0, 10);
  model::geo::addPoint(0.3, 0.1, 0, 1.0, 11);
  model::geo::addLine(8, 11, 10);
  model::geo::addLine(11, 10, 11);
  model::geo::addLine(10, 7, 12);
  model::geo::addLine(7, 8, 13);
  model::geo::addCurveLoop({13, 10, 11, 12}, 14);
  model::geo::addPlaneSurface({14}, 15);
  for(int i = 10; i <= 13; i++)
    model::geo::mesh::setTransfiniteCurve(i, 10);
  model::geo::mesh::setTransfiniteSurface(15);

  // The way triangles are generated can be controlled by specifying "Left",
  // "Right" or "Alternate" in `setTransfiniteSurface()' command. Try e.g.
  //
  // model::geo::mesh::setTransfiniteSurface(15, "Alternate");

  // Finally we apply an elliptic smoother to the grid to have a more regular
  // mesh:
  option::setNumber("Mesh.Smoothing", 100);

  model::geo::synchronize();
  model::mesh::generate(2);
  write("t6.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 7
//
//  Background meshes
//
// -----------------------------------------------------------------------------
// Mesh sizes can be specified very accurately by providing a background mesh,
// i.e., a post-processing view that contains the target mesh sizes.
void t7(int argc, char **argv)
{
    // Merge a list-based post-processing view containing the target mesh sizes:
  try {
    merge("../t7_bgmesh.pos");
  } catch(...) {
    logger::write("Could not load background mesh: bye!");
    finalize();
    return;
  }

  // If the post-processing view was model-based instead of list-based (i.e. if
  // it was based on an actual mesh), we would need to create a new model to
  // contain the geometry so that meshing it does not destroy the background
  // mesh. It's not necessary here since the view is list-based, but it does no
  // harm:
  model::add("t7");

  // Create a simple rectangular geometry:
  double lc = 1e-2;
  model::geo::addPoint(0, 0, 0, lc, 1);
  model::geo::addPoint(.1, 0, 0, lc, 2);
  model::geo::addPoint(.1, .3, 0, lc, 3);
  model::geo::addPoint(0, .3, 0, lc, 4);
  model::geo::addLine(1, 2, 1);
  model::geo::addLine(3, 2, 2);
  model::geo::addLine(3, 4, 3);
  model::geo::addLine(4, 1, 4);
  model::geo::addCurveLoop({4, 1, -2, 3}, 1);
  model::geo::addPlaneSurface({1}, 1);
  model::geo::synchronize();

  // Add the post-processing view as a new size field:
  int bg_field = model::mesh::field::add("PostView");
  model::mesh::field::setNumber(bg_field, "ViewIndex", 0);

  // Apply the view as the current background mesh size field:
  model::mesh::field::setAsBackgroundMesh(bg_field);

  // In order to compute the mesh sizes from the background mesh only, and
  // disregard any other size constraints, one can set:
  option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  option::setNumber("Mesh.MeshSizeFromPoints", 0);
  option::setNumber("Mesh.MeshSizeFromCurvature", 0);

  // See `t10.cpp' for additional information: background meshes are actually a
  // particular case of general "mesh size fields".

  model::mesh::generate(2);
  write("t7.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 8
//
//  Post-processing and animations
//
// -----------------------------------------------------------------------------
// In addition to creating geometries and meshes, the C++ API can also be used
// to manipulate post-processing datasets (called "views" in Gmsh).
void t8(int argc, char **argv)
{
    model::add("t8");

  // We first create a simple geometry
  double lc = 1e-2;
  model::geo::addPoint(0, 0, 0, lc, 1);
  model::geo::addPoint(.1, 0, 0, lc, 2);
  model::geo::addPoint(.1, .3, 0, lc, 3);
  model::geo::addPoint(0, .3, 0, lc, 4);
  model::geo::addLine(1, 2, 1);
  model::geo::addLine(3, 2, 2);
  model::geo::addLine(3, 4, 3);
  model::geo::addLine(4, 1, 4);
  model::geo::addCurveLoop({4, 1, -2, 3}, 1);
  model::geo::addPlaneSurface({1}, 1);
  model::geo::synchronize();

  // We merge some post-processing views to work on
  try {
    merge("../view1.pos");
    merge("../view1.pos");
    merge("../view4.pos"); // contains 2 views inside
  } catch(...) {
    logger::write("Could not load post-processing views: bye!");
    finalize();
    return;
  }

  // Gmsh can read post-processing views in various formats. Here the
  // `view1.pos' and `view4.pos' files are in the Gmsh "parsed" format, which is
  // interpreted by the GEO script parser. The parsed format should only be used
  // for relatively small datasets of course: for larger datasets using e.g. MSH
  // files is much more efficient. Post-processing views can also be created
  // directly from the C++ API.

  // We then set some general options:
  option::setNumber("General.Trackball", 0);
  option::setNumber("General.RotationX", 0);
  option::setNumber("General.RotationY", 0);
  option::setNumber("General.RotationZ", 0);

  int white[3] = {255, 255, 255};
  int black[3] = {0, 0, 0};
  option::setColor("General.Background", white[0], white[1], white[2]);

  // We can make our own shorter versions of repetitive methods
  auto set_color = [](string name, int c[3]) -> void {
    option::setColor(name, c[0], c[1], c[2]);
  };
  set_color("General.Foreground", black);
  set_color("General.Text", black);

  option::setNumber("General.Orthographic", 0);
  option::setNumber("General.Axes", 0);
  option::setNumber("General.SmallAxes", 0);

  // Show the GUI:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::initialize();

  // We also set some options for each post-processing view:
  option::setNumber("View[0].IntervalsType", 2);
  option::setNumber("View[0].OffsetZ", 0.05);
  option::setNumber("View[0].RaiseZ", 0);
  option::setNumber("View[0].Light", 1);
  option::setNumber("View[0].ShowScale", 0);
  option::setNumber("View[0].SmoothNormals", 1);

  option::setNumber("View[1].IntervalsType", 1);
  // We can't yet set the ColorTable through the API
  // option::setColorTable("View[1].ColorTable", "{ Green, Blue }");
  option::setNumber("View[1].NbIso", 10);
  option::setNumber("View[1].ShowScale", 0);

  option::setString("View[2].Name", "Test...");
  option::setNumber("View[2].Axes", 1);
  option::setNumber("View[2].IntervalsType", 2);
  option::setNumber("View[2].Type", 2);
  option::setNumber("View[2].IntervalsType", 2);
  option::setNumber("View[2].AutoPosition", 0);
  option::setNumber("View[2].PositionX", 85);
  option::setNumber("View[2].PositionY", 50);
  option::setNumber("View[2].Width", 200);
  option::setNumber("View[2].Height", 130);

  option::setNumber("View[3].Visible", 0);

  // You can save an MPEG movie directly by selecting `File->Export' in the
  // GUI. Several predefined animations are setup, for looping on all the time
  // steps in views, or for looping between views.

  // But the API can be used to build much more complex animations, by changing
  // options at run-time and re-rendering the graphics. Each frame can then be
  // saved to disk as an image, and multiple frames can be encoded to form a
  // movie. Below is an example of such a custom animation.

  int t = 0; // Initial step

  for(int num = 1; num <= 3; num++) {
    double nbt;
    option::getNumber("View[0].NbTimeStep", nbt);
    t = (t < nbt - 1) ? t + 1 : 0;

    // Set time step
    option::setNumber("View[0].TimeStep", t);
    option::setNumber("View[1].TimeStep", t);
    option::setNumber("View[2].TimeStep", t);
    option::setNumber("View[3].TimeStep", t);

    double max;
    option::getNumber("View[0].Max", max);
    option::setNumber("View[0].RaiseZ", 0.01 / max * t);

    if(num == 3) {
      double mw;
      option::getNumber("General.MenuWidth", mw);
      option::setNumber("General.GraphicsWidth", mw + 640);
      option::setNumber("General.GraphicsHeight", 480);
    }

    int frames = 50;
    for(int num2 = 1; num2 <= frames; num2++) {
      // Incrementally rotate the scene
      double rotx;
      option::getNumber("General.RotationX", rotx);
      option::setNumber("General.RotationX", rotx + 10);
      option::setNumber("General.RotationY", (rotx + 10) / 3.);
      double rotz;
      option::getNumber("General.RotationZ", rotz);
      option::setNumber("General.RotationZ", rotz + 0.1);

      // Draw the scene
      graphics::draw();

      if(num == 3) {
        // Uncomment the following lines to save each frame to an image file
        // write("t2-" + to_string(num2) + ".gif");
        // write("t2-" + to_string(num2) + ".ppm");
        // write("t2-" + to_string(num2) + ".jpg");
      }
    }

    if(num == 3) {
      // Here we could make a system call to generate a movie...
    }
  }

  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 9
//
//  Plugins
//
// -----------------------------------------------------------------------------

// Plugins can be added to Gmsh in order to extend its capabilities. For
// example, post-processing plugins can modify views, or create new views based
// on previously loaded views. Several default plugins are statically linked
// with Gmsh, e.g. Isosurface, CutPlane, CutSphere, Skin, Transform or Smooth.
//
// Plugins can be controlled through the API functions in the `plugin'
// namespace, or from the graphical interface (right click on the view button,
// then `Plugins').
void t9(int argc, char **argv)
{
    model::add("t9");

  // Let us for example include a three-dimensional scalar view:
  try {
    merge("../view3.pos");
  } catch(...) {
    logger::write("Could not load post-processing views: bye!");
    finalize();
    return;
  }

  // We then set some options for the `Isosurface' plugin (which extracts an
  // isosurface from a 3D scalar view), and run it:
  plugin::setNumber("Isosurface", "Value", 0.67); // Iso-value level
  plugin::setNumber("Isosurface", "View", 0); // Source view is View[0]
  plugin::run("Isosurface"); // Run the plugin!

  // We also set some options for the `CutPlane' plugin (which computes a
  // section of a 3D view using the plane A*x+B*y+C*z+D=0), and then run it:
  plugin::setNumber("CutPlane", "A", 0);
  plugin::setNumber("CutPlane", "B", 0.2);
  plugin::setNumber("CutPlane", "C", 1);
  plugin::setNumber("CutPlane", "D", 0);
  plugin::setNumber("CutPlane", "View", 0);
  plugin::run("CutPlane");

  // Add a title (By convention, for window coordinates a value greater than
  // 99999 represents the center. We could also use `General.GraphicsWidth / 2',
  // but that would only center the string for the current window size.):
  plugin::setString("Annotate", "Text", "A nice title");
  plugin::setNumber("Annotate", "X", 1.e5);
  plugin::setNumber("Annotate", "Y", 50);
  plugin::setString("Annotate", "Font", "Times-BoldItalic");
  plugin::setNumber("Annotate", "FontSize", 28);
  plugin::setString("Annotate", "Align", "Center");
  plugin::setNumber("Annotate", "View", 0);
  plugin::run("Annotate");

  plugin::setString("Annotate", "Text", "(and a small subtitle)");
  plugin::setNumber("Annotate", "Y", 70);
  plugin::setString("Annotate", "Font", "Times-Roman");
  plugin::setNumber("Annotate", "FontSize", 12);
  plugin::run("Annotate");

  // We finish by setting some options:
  option::setNumber("View[0].Light", 1);
  option::setNumber("View[0].IntervalsType", 1);
  option::setNumber("View[0].NbIso", 6);
  option::setNumber("View[0].SmoothNormals", 1);
  option::setNumber("View[1].IntervalsType", 2);
  option::setNumber("View[2].IntervalsType", 2);

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 10
//
//  Mesh size fields
//
// -----------------------------------------------------------------------------

// In addition to specifying target mesh sizes at the points of the geometry
// (see `t1.cpp') or using a background mesh (see `t7.cpp'), you can use general
// mesh size "Fields".
void t10(int argc, char **argv)
{
    model::add("t10");

  // Let's create a simple rectangular geometry:
  double lc = .15;
  model::geo::addPoint(0.0, 0.0, 0, lc, 1);
  model::geo::addPoint(1, 0.0, 0, lc, 2);
  model::geo::addPoint(1, 1, 0, lc, 3);
  model::geo::addPoint(0, 1, 0, lc, 4);
  model::geo::addPoint(0.2, .5, 0, lc, 5);

  model::geo::addLine(1, 2, 1);
  model::geo::addLine(2, 3, 2);
  model::geo::addLine(3, 4, 3);
  model::geo::addLine(4, 1, 4);

  model::geo::addCurveLoop({1, 2, 3, 4}, 5);
  model::geo::addPlaneSurface({5}, 6);

  model::geo::synchronize();

  // Say we would like to obtain mesh elements with size lc/30 near curve 2 and
  // point 5, and size lc elsewhere. To achieve this, we can use two fields:
  // "Distance", and "Threshold". We first define a Distance field (`Field[1]')
  // on points 5 and on curve 2. This field returns the distance to point 5 and
  // to (100 equidistant points on) curve 2.
  model::mesh::field::add("Distance", 1);
  model::mesh::field::setNumbers(1, "PointsList", {5});
  model::mesh::field::setNumbers(1, "CurvesList", {2});
  model::mesh::field::setNumber(1, "Sampling", 100);

  // We then define a `Threshold' field, which uses the return value of the
  // `Distance' field 1 in order to define a simple change in element size
  // depending on the computed distances
  //
  // SizeMax -                     /------------------
  //                              /
  //                             /
  //                            /
  // SizeMin -o----------------/
  //          |                |    |
  //        Point         DistMin  DistMax
  model::mesh::field::add("Threshold", 2);
  model::mesh::field::setNumber(2, "InField", 1);
  model::mesh::field::setNumber(2, "SizeMin", lc / 30);
  model::mesh::field::setNumber(2, "SizeMax", lc);
  model::mesh::field::setNumber(2, "DistMin", 0.15);
  model::mesh::field::setNumber(2, "DistMax", 0.5);

  // Say we want to modulate the mesh element sizes using a mathematical
  // function of the spatial coordinates. We can do this with the MathEval
  // field:
  model::mesh::field::add("MathEval", 3);
  model::mesh::field::setString(
    3, "F", "Cos(4*3.14*x) * Sin(4*3.14*y) / 10 + 0.101");

  // We could also combine MathEval with values coming from other fields. For
  // example, let's define a `Distance' field around point 1
  model::mesh::field::add("Distance", 4);
  model::mesh::field::setNumbers(4, "PointsList", {1});

  // We can then create a `MathEval' field with a function that depends on the
  // return value of the `Distance' field 4, i.e., depending on the distance to
  // point 1 (here using a cubic law, with minimum element size = lc / 100)
  model::mesh::field::add("MathEval", 5);
  stringstream stream;
  stream << "F4^3 + " << lc / 100;
  model::mesh::field::setString(5, "F", stream.str());

  // We could also use a `Box' field to impose a step change in element sizes
  // inside a box
  model::mesh::field::add("Box", 6);
  model::mesh::field::setNumber(6, "VIn", lc / 15);
  model::mesh::field::setNumber(6, "VOut", lc);
  model::mesh::field::setNumber(6, "XMin", 0.3);
  model::mesh::field::setNumber(6, "XMax", 0.6);
  model::mesh::field::setNumber(6, "YMin", 0.3);
  model::mesh::field::setNumber(6, "YMax", 0.6);
  model::mesh::field::setNumber(6, "Thickness", 0.3);

  // Many other types of fields are available: see the reference manual for a
  // complete list. You can also create fields directly in the graphical user
  // interface by selecting `Define->Size fields' in the `Mesh' module.

  // Let's use the minimum of all the fields as the background mesh field:
  model::mesh::field::add("Min", 7);
  model::mesh::field::setNumbers(7, "FieldsList", {2, 3, 5, 6});

  model::mesh::field::setAsBackgroundMesh(7);

  // The API also allows to set a global mesh size callback, which is called
  // each time the mesh size is queried
  auto meshSizeCallback = [](int dim, int tag, double x, double y, double z,
                             double lc) {
    return min(lc, 0.02 * x + 0.01);
  };
  model::mesh::setSizeCallback(meshSizeCallback);

  // To determine the size of mesh elements, Gmsh locally computes the minimum
  // of
  //
  // 1) the size of the model bounding box;
  // 2) if `Mesh.MeshSizeFromPoints' is set, the mesh size specified at
  //    geometrical points;
  // 3) if `Mesh.MeshSizeFromCurvature' is positive, the mesh size based on
  //    curvature (the value specifying the number of elements per 2 * pi rad);
  // 4) the background mesh field;
  // 5) any per-entity mesh size constraint;
  //
  // The value can then be further modified by the mesh size callback, if any,
  // before being constrained in the interval [`Mesh.MeshSizeMin',
  // `Mesh.MeshSizeMax'] and multiplied by `Mesh.MeshSizeFactor'. In addition,
  // boundary mesh sizes are interpolated inside surfaces and/or volumes
  // depending on the value of `Mesh.MeshSizeExtendFromBoundary' (which is set
  // by default).
  //
  // When the element size is fully specified by a background mesh (as it is in
  // this example), it is thus often desirable to set

  option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  option::setNumber("Mesh.MeshSizeFromPoints", 0);
  option::setNumber("Mesh.MeshSizeFromCurvature", 0);

  // This will prevent over-refinement due to small mesh sizes on the boundary.

  // Finally, while the default "Frontal-Delaunay" 2D meshing algorithm
  // (Mesh.Algorithm = 6) usually leads to the highest quality meshes, the
  // "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size
  // fields better - in particular size fields with large element size
  // gradients:

  option::setNumber("Mesh.Algorithm", 5);

  model::mesh::generate(2);
  write("t10.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 11
//
//  Unstructured quadrangular meshes
//
// -----------------------------------------------------------------------------
void t11(int argc, char **argv)
{
    model::add("t11");

  // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
  // transfinite meshes can be "recombined" into quads, prisms or
  // hexahedra. Unstructured meshes can be recombined in the same way. Let's
  // define a simple geometry with an analytical mesh size field:

  int p1 = model::geo::addPoint(-1.25, -.5, 0);
  int p2 = model::geo::addPoint(1.25, -.5, 0);
  int p3 = model::geo::addPoint(1.25, 1.25, 0);
  int p4 = model::geo::addPoint(-1.25, 1.25, 0);

  int l1 = model::geo::addLine(p1, p2);
  int l2 = model::geo::addLine(p2, p3);
  int l3 = model::geo::addLine(p3, p4);
  int l4 = model::geo::addLine(p4, p1);

  int cl = model::geo::addCurveLoop({l1, l2, l3, l4});
  int pl = model::geo::addPlaneSurface({cl});

  model::geo::synchronize();

  model::mesh::field::add("MathEval", 1);
  model::mesh::field::setString(
    1, "F", "0.01*(1.0+30.*(y-x*x)*(y-x*x) + (1-x)*(1-x))");
  model::mesh::field::setAsBackgroundMesh(1);

  // To generate quadrangles instead of triangles, we can simply add
  model::mesh::setRecombine(2, pl);

  // If we'd had several surfaces, we could have used the global option
  // "Mesh.RecombineAll":
  //
  // option::setNumber("Mesh.RecombineAll", 1);

  // The default recombination algorithm is called "Blossom": it uses a minimum
  // cost perfect matching algorithm to generate fully quadrilateral meshes from
  // triangulations. More details about the algorithm can be found in the
  // following paper: J.-F. Remacle, J. Lambrechts, B. Seny, E. Marchandise,
  // A. Johnen and C. Geuzaine, "Blossom-Quad: a non-uniform quadrilateral mesh
  // generator using a minimum cost perfect matching algorithm", International
  // Journal for Numerical Methods in Engineering 89, pp. 1102-1119, 2012.

  // For even better 2D (planar) quadrilateral meshes, you can try the
  // experimental "Frontal-Delaunay for quads" meshing algorithm, which is a
  // triangulation algorithm that enables to create right triangles almost
  // everywhere: J.-F. Remacle, F. Henrotte, T. Carrier-Baudouin, E. Bechet,
  // E. Marchandise, C. Geuzaine and T. Mouton. A frontal Delaunay quad mesh
  // generator using the L^inf norm. International Journal for Numerical Methods
  // in Engineering, 94, pp. 494-512, 2013. Uncomment the following line to try
  // the Frontal-Delaunay algorithms for quads:
  //
  // option::setNumber("Mesh.Algorithm", 8);

  // The default recombination algorithm might leave some triangles in the mesh,
  // if recombining all the triangles leads to badly shaped quads. In such
  // cases, to generate full-quad meshes, you can either subdivide the resulting
  // hybrid mesh (with `Mesh.SubdivisionAlgorithm' set to 1), or use the
  // full-quad recombination algorithm, which will automatically perform a
  // coarser mesh followed by recombination, smoothing and
  // subdivision. Uncomment the following line to try the full-quad algorithm:
  //
  // option::setNumber("Mesh.RecombinationAlgorithm", 2); // or 3

  // You can also set the subdivision step alone, with
  //
  // option::setNumber("Mesh.SubdivisionAlgorithm", 1);

  model::mesh::generate(2);

  // Note that you could also apply the recombination algorithm and/or the
  // subdivision step explicitly after meshing, as follows:
  //
  // model::mesh::generate(2);
  // model::mesh::recombine();
  // option::setNumber("Mesh.SubdivisionAlgorithm", 1);
  // model::mesh::refine();

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 12
//
//  Cross-patch meshing with compounds
//
// -----------------------------------------------------------------------------
// "Compound" meshing constraints allow to generate meshes across surface
// boundaries, which can be useful e.g. for imported CAD models (e.g. STEP) with
// undesired small features.

// When a `setCompound()' meshing constraint is given, at mesh generation time
// Gmsh
//  1. meshes the underlying elementary geometrical entities, individually
//  2. creates a discrete entity that combines all the individual meshes
//  3. computes a discrete parametrization (i.e. a piece-wise linear mapping)
//     on this discrete entity
//  4. meshes the discrete entity using this discrete parametrization instead
//     of the underlying geometrical description of the underlying elementary
//     entities making up the compound
//  5. optionally, reclassifies the mesh elements and nodes on the original
//     entities

// Step 3. above can only be performed if the mesh resulting from the
// combination of the individual meshes can be reparametrized, i.e. if the shape
// is "simple enough". If the shape is not amenable to reparametrization, you
// should create a full mesh of the geometry and first re-classify it to
// generate patches amenable to reparametrization (see `t13.cpp').

// The mesh of the individual entities performed in Step 1. should usually be
// finer than the desired final mesh; this can be controlled with the
// `Mesh.CompoundMeshSizeFactor' option.

// The optional reclassification on the underlying elementary entities in Step
// 5. is governed by the `Mesh.CompoundClassify' option.
void t12(int argc, char **argv)
{
    model::add("t12");

  double lc = 0.1;

  model::geo::addPoint(0, 0, 0, lc, 1);
  model::geo::addPoint(1, 0, 0, lc, 2);
  model::geo::addPoint(1, 1, 0.5, lc, 3);
  model::geo::addPoint(0, 1, 0.4, lc, 4);
  model::geo::addPoint(0.3, 0.2, 0, lc, 5);
  model::geo::addPoint(0, 0.01, 0.01, lc, 6);
  model::geo::addPoint(0, 0.02, 0.02, lc, 7);
  model::geo::addPoint(1, 0.05, 0.02, lc, 8);
  model::geo::addPoint(1, 0.32, 0.02, lc, 9);

  model::geo::addLine(1, 2, 1);
  model::geo::addLine(2, 8, 2);
  model::geo::addLine(8, 9, 3);
  model::geo::addLine(9, 3, 4);
  model::geo::addLine(3, 4, 5);
  model::geo::addLine(4, 7, 6);
  model::geo::addLine(7, 6, 7);
  model::geo::addLine(6, 1, 8);
  model::geo::addSpline({7, 5, 9}, 9);
  model::geo::addLine(6, 8, 10);

  model::geo::addCurveLoop({5, 6, 9, 4}, 11);
  model::geo::addSurfaceFilling({11}, 1);

  model::geo::addCurveLoop({-9, 3, 10, 7}, 13);
  model::geo::addSurfaceFilling({13}, 5);

  model::geo::addCurveLoop({-10, 2, 1, 8}, 15);
  model::geo::addSurfaceFilling({15}, 10);

  model::geo::synchronize();

  // Treat curves 2, 3 and 4 as a single curve when meshing (i.e. mesh across
  // points 6 and 7)
  model::mesh::setCompound(1, {2, 3, 4});

  // Idem with curves 6, 7 and 8
  model::mesh::setCompound(1, {6, 7, 8});

  // Treat surfaces 1, 5 and 10 as a single surface when meshing (i.e. mesh
  // across curves 9 and 10)
  model::mesh::setCompound(2, {1, 5, 10});

  model::mesh::generate(2);
  write("t12.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 13
//
//  Remeshing an STL file without an underlying CAD model
//
// -----------------------------------------------------------------------------
void t13(int argc, char **argv)
{
    auto createGeometryAndMesh = []()
  {
    model::add("t13");

    // Let's merge an STL mesh that we would like to remesh (from the parent
    // directory):
    try {
      merge("../t13_data.stl");
    } catch(...) {
      logger::write("Could not load STL mesh: bye!");
      return;
    }

    // We first classify ("color") the surfaces by splitting the original
    // surface along sharp geometrical features. This will create new discrete
    // surfaces, curves and points.

    // Angle between two triangles above which an edge is considered as sharp,
    // retrieved from the ONELAB database (see below):
    vector<double> n;
    onelab::getNumber("Parameters/Angle for surface detection", n);
    double angle = n[0];

    // For complex geometries, patches can be too complex, too elongated or too
    // large to be parametrized; setting the following option will force the
    // creation of patches that are amenable to reparametrization:
    onelab::getNumber
      ("Parameters/Create surfaces guaranteed to be parametrizable", n);
    bool forceParametrizablePatches = n[0] ? true : false;

    // For open surfaces include the boundary edges in the classification process:
    bool includeBoundary = true;

    // Force curves to be split on given angle:
    double curveAngle = 180;

    model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                        forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);

    // Create a geometry for all the discrete curves and surfaces in the mesh,
    // by computing a parametrization for each one
    model::mesh::createGeometry();

    // Note that if a CAD model (e.g. as a STEP file, see `t20.cpp') is
    // available instead of an STL mesh, it is usually better to use that CAD
    // model instead of the geometry created by reparametrizing the
    // mesh. Indeed, CAD geometries will in general be more accurate, with
    // smoother parametrizations, and will lead to more efficient and higher
    // quality meshing. Discrete surface remeshing in Gmsh is optimized to
    // handle dense STL meshes coming from e.g. imaging systems, where no CAD is
    // available; it is less well suited for the poor quality STL triangulations
    // (optimized for size, with e.g. very elongated triangles) that are usually
    // generated by CAD tools for e.g. 3D printing.

    // Create a volume from all the surfaces
    vector<pair<int, int> > s;
    model::getEntities(s, 2);
    vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = model::geo::addSurfaceLoop(sl);
    model::geo::addVolume({l});

    model::geo::synchronize();

    // We specify element sizes imposed by a size field, just because we can :-)
    int f = model::mesh::field::add("MathEval");
    onelab::getNumber("Parameters/Apply funny mesh size field?", n);
    if(n[0])
      model::mesh::field::setString(f, "F", "2*Sin((x+y)/5) + 3");
    else
      model::mesh::field::setString(f, "F", "4");
    model::mesh::field::setAsBackgroundMesh(f);

    model::mesh::generate(3);
  };

  // Create ONELAB parameters with remeshing options:
  onelab::set(R"( [
  {
    "type":"number",
    "name":"Parameters/Angle for surface detection",
    "values":[40],
    "min":20,
    "max":120,
    "step":1
  },
  {
    "type":"number",
    "name":"Parameters/Create surfaces guaranteed to be parametrizable",
    "values":[0],
    "choices":[0, 1]
  },
  {
    "type":"number",
    "name":"Parameters/Apply funny mesh size field?",
    "values":[0],
    "choices":[0, 1]
  }
  ] )");

  // Create the geometry and mesh it:
  createGeometryAndMesh();

  // Launch the GUI and handle the "check" event to recreate the geometry and mesh
  // with new parameters if necessary:
  auto checkForEvent = [=]() -> bool {
    vector<string> action;
    onelab::getString("ONELAB/Action", action);
    if(action.size() && action[0] == "check") {
      onelab::setString("ONELAB/Action", {""});
      createGeometryAndMesh();
      graphics::draw();
    }
    return true;
  };

  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) {
    fltk::initialize();
    while(fltk::isAvailable() && checkForEvent())
      fltk::wait();
  }
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 14
//
//  Homology and cohomology computation
//
// -----------------------------------------------------------------------------

// Homology computation in Gmsh finds representative chains of (relative)
// (co)homology space bases using a mesh of a model.  The representative basis
// chains are stored in the mesh as physical groups of Gmsh, one for each chain.
void t14(int argc, char **argv)
{
    model::add("t14");

  // Create an example geometry

  double m = 0.5; // mesh size
  double h = 2; // geometry height in the z-direction

  model::geo::addPoint(0, 0, 0, m, 1);
  model::geo::addPoint(10, 0, 0, m, 2);
  model::geo::addPoint(10, 10, 0, m, 3);
  model::geo::addPoint(0, 10, 0, m, 4);

  model::geo::addPoint(4, 4, 0, m, 5);
  model::geo::addPoint(6, 4, 0, m, 6);
  model::geo::addPoint(6, 6, 0, m, 7);
  model::geo::addPoint(4, 6, 0, m, 8);

  model::geo::addPoint(2, 0, 0, m, 9);
  model::geo::addPoint(8, 0, 0, m, 10);
  model::geo::addPoint(2, 10, 0, m, 11);
  model::geo::addPoint(8, 10, 0, m, 12);

  model::geo::addLine(1, 9, 1);
  model::geo::addLine(9, 10, 2);
  model::geo::addLine(10, 2, 3);

  model::geo::addLine(2, 3, 4);
  model::geo::addLine(3, 12, 5);
  model::geo::addLine(12, 11, 6);

  model::geo::addLine(11, 4, 7);
  model::geo::addLine(4, 1, 8);
  model::geo::addLine(5, 6, 9);

  model::geo::addLine(6, 7, 10);
  model::geo::addLine(7, 8, 11);
  model::geo::addLine(8, 5, 12);

  model::geo::addCurveLoop({6, 7, 8, 1, 2, 3, 4, 5}, 13);
  model::geo::addCurveLoop({11, 12, 9, 10}, 14);
  model::geo::addPlaneSurface({13, 14}, 15);

  vector<pair<int, int> > e;
  model::geo::extrude({{2, 15}}, 0, 0, h, e);

  model::geo::synchronize();

  // Create physical groups, which are used to define the domain of the
  // (co)homology computation and the subdomain of the relative (co)homology
  // computation.

  // Whole domain
  int domain_tag = e[1].second;
  int domain_physical_tag = 1001;
  model::addPhysicalGroup(3, {domain_tag}, domain_physical_tag);
  model::setPhysicalName(3, domain_physical_tag, "Whole domain");

  // Four "terminals" of the model
  vector<int> terminal_tags = {e[3].second, e[5].second, e[7].second,
                                    e[9].second};
  int terminals_physical_tag = 2001;
  model::addPhysicalGroup(2, terminal_tags, terminals_physical_tag);
  model::setPhysicalName(2, terminals_physical_tag, "Terminals");

  // Find domain boundary tags
  vector<pair<int, int> > boundary_dimtags;
  model::getBoundary({{3, domain_tag}}, boundary_dimtags, false, false);

  vector<int> boundary_tags, complement_tags;
  for(auto e : boundary_dimtags) {
    complement_tags.push_back(e.second);
    boundary_tags.push_back(e.second);
  }
  for(auto t : terminal_tags) {
    auto it = find(complement_tags.begin(), complement_tags.end(), t);
    if(it != complement_tags.end()) complement_tags.erase(it);
  }

  // Whole domain surface
  int boundary_physical_tag = 2002;
  model::addPhysicalGroup(2, boundary_tags, boundary_physical_tag);
  model::setPhysicalName(2, boundary_physical_tag, "Boundary");

  // Complement of the domain surface with respect to the four terminals
  int complement_physical_tag = 2003;
  model::addPhysicalGroup(2, complement_tags, complement_physical_tag);
  model::setPhysicalName(2, complement_physical_tag, "Complement");

  // Find bases for relative homology spaces of the domain modulo the four
  // terminals
  model::mesh::computeHomology({domain_physical_tag},
                                     {terminals_physical_tag}, {0, 1, 2, 3});

  // Find homology space bases isomorphic to the previous bases: homology spaces
  // modulo the non-terminal domain surface, a.k.a the thin cuts
  model::mesh::computeHomology({domain_physical_tag},
                                     {complement_physical_tag}, {0, 1, 2, 3});

  // Find cohomology space bases isomorphic to the previous bases: cohomology
  // spaces of the domain modulo the four terminals, a.k.a the thick cuts
  model::mesh::computeCohomology({domain_physical_tag},
                                       {terminals_physical_tag}, {0, 1, 2, 3});

  // More examples:
  // model::mesh::computeHomology();
  // model::mesh::computeHomology({domain_physical_tag});
  // model::mesh::computeHomology({domain_physical_tag},
  //                                    {boundary_physical_tag},
  //                                     {0,1,2,3});

  // Generate the mesh and perform the requested homology computations
  model::mesh::generate(3);

  // For more information, see M. Pellikka, S. Suuriniemi, L. Kettunen and
  // C. Geuzaine. Homology and cohomology computation in finite element
  // modeling. SIAM Journal on Scientific Computing 35(5), pp. 1195-1214, 2013.

  write("t14.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 15
//
//  Embedded points, lines and surfaces
//
// -----------------------------------------------------------------------------

// By default, across geometrical dimensions meshes generated by Gmsh are only
// conformal if lower dimensional entities are on the boundary of higher
// dimensional ones (i.e. if points, curves or surfaces are part of the boundary
// of volumes).

// Embedding constraints allow to force a mesh to be conformal to other lower
// dimensional entities.
void t15(int argc, char **argv)
{
    model::add("t15");

  // Copied from t1.cpp:
  double lc = 1e-2;
  model::geo::addPoint(0, 0, 0, lc, 1);
  model::geo::addPoint(.1, 0, 0, lc, 2);
  model::geo::addPoint(.1, .3, 0, lc, 3);
  model::geo::addPoint(0, .3, 0, lc, 4);
  model::geo::addLine(1, 2, 1);
  model::geo::addLine(3, 2, 2);
  model::geo::addLine(3, 4, 3);
  model::geo::addLine(4, 1, 4);
  model::geo::addCurveLoop({4, 1, -2, 3}, 1);
  model::geo::addPlaneSurface({1}, 1);

  // We change the mesh size to generate a coarser mesh
  lc *= 4;
  model::geo::mesh::setSize({{0, 1}, {0, 2}, {0, 3}, {0, 4}}, lc);

  // We define a new point
  model::geo::addPoint(0.02, 0.02, 0., lc, 5);

  // We have to synchronize before embedding entites:
  model::geo::synchronize();

  // One can force this point to be included ("embedded") in the 2D mesh, using
  // the `embed()' function:
  model::mesh::embed(0, {5}, 2, 1);

  // In the same way, one can use `embed()' to force a curve to be embedded in
  // the 2D mesh:
  model::geo::addPoint(0.02, 0.12, 0., lc, 6);
  model::geo::addPoint(0.04, 0.18, 0., lc, 7);
  model::geo::addLine(6, 7, 5);
  model::geo::synchronize();
  model::mesh::embed(1, {5}, 2, 1);

  // Points and curves can also be embedded in volumes
  vector<pair<int, int> > ext;
  model::geo::extrude({{2, 1}}, 0, 0, 0.1, ext);

  int p = model::geo::addPoint(0.07, 0.15, 0.025, lc);

  model::geo::synchronize();
  model::mesh::embed(0, {p}, 3, 1);

  model::geo::addPoint(0.025, 0.15, 0.025, lc, p + 1);
  int l = model::geo::addLine(7, p + 1);

  model::geo::synchronize();
  model::mesh::embed(1, {l}, 3, 1);

  // Finally, we can also embed a surface in a volume:
  model::geo::addPoint(0.02, 0.12, 0.05, lc, p + 2);
  model::geo::addPoint(0.04, 0.12, 0.05, lc, p + 3);
  model::geo::addPoint(0.04, 0.18, 0.05, lc, p + 4);
  model::geo::addPoint(0.02, 0.18, 0.05, lc, p + 5);

  model::geo::addLine(p + 2, p + 3, l + 1);
  model::geo::addLine(p + 3, p + 4, l + 2);
  model::geo::addLine(p + 4, p + 5, l + 3);
  model::geo::addLine(p + 5, p + 2, l + 4);

  int ll = model::geo::addCurveLoop({l + 1, l + 2, l + 3, l + 4});
  int s = model::geo::addPlaneSurface({ll});

  model::geo::synchronize();
  model::mesh::embed(2, {s}, 3, 1);

  // Note that with the OpenCASCADE kernel (see `t16.cpp'), when the
  // `fragment()' function is applied to entities of different dimensions, the
  // lower dimensional entities will be autmatically embedded in the higher
  // dimensional entities if necessary.

  model::mesh::generate(3);

  write("t15.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 16
//
//  Constructive Solid Geometry, OpenCASCADE geometry kernel
//
// -----------------------------------------------------------------------------

// Instead of constructing a model in a bottom-up fashion with Gmsh's built-in
// geometry kernel, starting with version 3 Gmsh allows you to directly use
// alternative geometry kernels. Here we will use the OpenCASCADE kernel.
void t16(int argc, char **argv)
{
      model::add("t16");

  // Let's build the same model as in `t5.cpp', but using constructive solid
  // geometry.

  // We can log all messages for further processing with:
  logger::start();

  // We first create two cubes:
  try {
    model::occ::addBox(0, 0, 0, 1, 1, 1, 1);
    model::occ::addBox(0, 0, 0, 0.5, 0.5, 0.5, 2);
  } catch(...) {
    logger::write("Could not create OpenCASCADE shapes: bye!");
    return;
  }

  // We apply a boolean difference to create the "cube minus one eigth" shape:
  vector<pair<int, int> > ov;
  vector<vector<pair<int, int> > > ovv;
  model::occ::cut({{3, 1}}, {{3, 2}}, ov, ovv, 3);

  // Boolean operations with OpenCASCADE always create new entities. By default
  // the extra arguments `removeObject' and `removeTool' in `cut()' are set to
  // `true', which will delete the original entities.

  // We then create the five spheres:
  double x = 0, y = 0.75, z = 0, r = 0.09;
  vector<pair<int, int> > holes;
  for(int t = 1; t <= 5; t++) {
    x += 0.166;
    z += 0.166;
    model::occ::addSphere(x, y, z, r, 3 + t);
    holes.push_back({3, 3 + t});
  }

  // If we had wanted five empty holes we would have used `cut()' again. Here we
  // want five spherical inclusions, whose mesh should be conformal with the
  // mesh of the cube: we thus use `fragment()', which intersects all volumes in
  // a conformal manner (without creating duplicate interfaces):
  model::occ::fragment({{3, 3}}, holes, ov, ovv);

  // ov contains all the generated entities of the same dimension as the input
  // entities:
  logger::write("fragment produced volumes:");
  for(auto e : ov)
    logger::write("(" + to_string(e.first) + "," +
                        to_string(e.second) + ")");

  // ovv contains the parent-child relationships for all the input entities:
  logger::write("before/after volume relations:");
  vector<pair<int, int> > in(1, pair<int, int>(3, 3));
  in.insert(in.end(), holes.begin(), holes.end());
  for(size_t i = 0; i < in.size(); i++) {
    string s = "parent (" + to_string(in[i].first) + "," +
                    to_string(in[i].second) + ") -> child";
    for(size_t j = 0; j < ovv[i].size(); j++) {
      s += " (" + to_string(ovv[i][j].first) + "," +
           to_string(ovv[i][j].second) + ")";
    }
    logger::write(s);
  }

  model::occ::synchronize();

  // When the boolean operation leads to simple modifications of entities, and
  // if one deletes the original entities, Gmsh tries to assign the same tag to
  // the new entities. (This behavior is governed by the
  // `Geometry.OCCBooleanPreserveNumbering' option.)

  // Here the `Physical Volume' definitions can thus be made for the 5 spheres
  // directly, as the five spheres (volumes 4, 5, 6, 7 and 8), which will be
  // deleted by the fragment operations, will be recreated identically (albeit
  // with new surfaces) with the same tags:
  for(int i = 1; i <= 5; i++) model::addPhysicalGroup(3, {3 + i}, i);

  // The tag of the cube will change though, so we need to access it
  // programmatically:
  model::addPhysicalGroup(3, {ov[ov.size() - 1].second}, 10);

  // Creating entities using constructive solid geometry is very powerful, but
  // can lead to practical issues for e.g. setting mesh sizes at points, or
  // identifying boundaries.

  // To identify points or other bounding entities you can take advantage of the
  // `getEntities()', `getBoundary()' and `getEntitiesInBoundingBox()'
  // functions:

  double lcar1 = .1;
  double lcar2 = .0005;
  double lcar3 = .055;

  // Assign a mesh size to all the points:
  model::getEntities(ov, 0);
  model::mesh::setSize(ov, lcar1);

  // Override this constraint on the points of the five spheres:
  model::getBoundary(holes, ov, false, false, true);
  model::mesh::setSize(ov, lcar3);

  // Select the corner point by searching for it geometrically:
  double eps = 1e-3;
  model::getEntitiesInBoundingBox(0.5 - eps, 0.5 - eps, 0.5 - eps,
                                        0.5 + eps, 0.5 + eps, 0.5 + eps, ov, 0);
  model::mesh::setSize(ov, lcar2);

  model::mesh::generate(3);

  write("t16.msh");

  // Additional examples created with the OpenCASCADE geometry kernel are
  // available in `t18.cpp', `t19.cpp' and `t20.cpp', as well as in the
  // `examples/api' directory.

  // Inspect the log:
  vector<string> log;
  logger::get(log);
  cout << "Logger has recorded " << log.size() << " lines" << endl;
  logger::stop();

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 17
//
//  Anisotropic background mesh
//
// -----------------------------------------------------------------------------

// As seen in `t7.cpp', mesh sizes can be specified very accurately by providing
// a background mesh, i.e., a post-processing view that contains the target mesh
// sizes.

// Here, the background mesh is represented as a metric tensor field defined on
// a square. One should use bamg as 2d mesh generator to enable anisotropic
// meshes in 2D.
void t17(int argc, char **argv)
{
    model::add("t17");

  // Create a square
  model::occ::addRectangle(-2, -2, 0, 4, 4);
  model::occ::synchronize();

  // Merge a post-processing view containing the target anisotropic mesh sizes
  try {
    merge("../t17_bgmesh.pos");
  } catch(...) {
    logger::write("Could not load background mesh: bye!");
    finalize();
    return;
  }

  // Apply the view as the current background mesh
  int bg_field = model::mesh::field::add("PostView");
  model::mesh::field::setNumber(bg_field, "ViewIndex", 0);
  model::mesh::field::setAsBackgroundMesh(bg_field);

  // Use bamg
  option::setNumber("Mesh.SmoothRatio", 3);
  option::setNumber("Mesh.AnisoMax", 1000);
  option::setNumber("Mesh.Algorithm", 7);

  model::mesh::generate(2);

  write("t17.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 18
//
//  Periodic meshes
//
// -----------------------------------------------------------------------------

// Periodic meshing constraints can be imposed on surfaces and curves.
void t18(int argc, char **argv)
{
    model::add("t18");

  // Let's use the OpenCASCADE geometry kernel to build two geometries.

  // The first geometry is very simple: a unit cube with a non-uniform mesh size
  // constraint (set on purpose to be able to verify visually that the
  // periodicity constraint works!):

  model::occ::addBox(0, 0, 0, 1, 1, 1, 1);
  model::occ::synchronize();

  vector<pair<int, int> > out;
  model::getEntities(out, 0);
  model::mesh::setSize(out, 0.1);
  model::mesh::setSize({{0, 1}}, 0.02);

  // To impose that the mesh on surface 2 (the right side of the cube) should
  // match the mesh from surface 1 (the left side), the following periodicity
  // constraint is set:
  vector<double> translation(
    {1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});
  model::mesh::setPeriodic(2, {2}, {1}, translation);

  // The periodicity transform is provided as a 4x4 affine transformation
  // matrix, given by row.

  // During mesh generation, the mesh on surface 2 will be created by copying
  // the mesh from surface 1.

  // Multiple periodicities can be imposed in the same way:
  model::mesh::setPeriodic(
    2, {6}, {5}, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1});
  model::mesh::setPeriodic(
    2, {4}, {3}, {1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1});

  // For more complicated cases, finding the corresponding surfaces by hand can
  // be tedious, especially when geometries are created through solid
  // modelling. Let's construct a slightly more complicated geometry.

  // We start with a cube and some spheres:
  model::occ::addBox(2, 0, 0, 1, 1, 1, 10);
  double x = 2 - 0.3, y = 0, z = 0;
  model::occ::addSphere(x, y, z, 0.35, 11);
  model::occ::addSphere(x + 1, y, z, 0.35, 12);
  model::occ::addSphere(x, y + 1, z, 0.35, 13);
  model::occ::addSphere(x, y, z + 1, 0.35, 14);
  model::occ::addSphere(x + 1, y + 1, z, 0.35, 15);
  model::occ::addSphere(x, y + 1, z + 1, 0.35, 16);
  model::occ::addSphere(x + 1, y, z + 1, 0.35, 17);
  model::occ::addSphere(x + 1, y + 1, z + 1, 0.35, 18);

  // We first fragment all the volumes, which will leave parts of spheres
  // protruding outside the cube:
  vector<vector<pair<int, int> > > out_map;
  vector<pair<int, int> > sph;
  for(int i = 11; i <= 18; i++) sph.push_back(pair<int, int>(3, i));
  model::occ::fragment({{3, 10}}, sph, out, out_map);
  model::occ::synchronize();

  // Ask OpenCASCADE to compute more accurate bounding boxes of entities using
  // the STL mesh:
  option::setNumber("Geometry.OCCBoundsUseStl", 1);

  // We then retrieve all the volumes in the bounding box of the original cube,
  // and delete all the parts outside it:
  double eps = 1e-3;
  vector<pair<int, int> > in;
  model::getEntitiesInBoundingBox(2 - eps, -eps, -eps, 2 + 1 + eps,
                                        1 + eps, 1 + eps, in, 3);
  for(auto i : in) {
    auto it = find(out.begin(), out.end(), i);
    if(it != out.end()) out.erase(it);
  }
  model::removeEntities(out, true); // Delete outside parts recursively

  // We now set a non-uniform mesh size constraint (again to check results
  // visually):
  vector<pair<int, int> > p;
  model::getBoundary(in, p, false, false, true); // Get all points
  model::mesh::setSize(p, 0.1);
  model::getEntitiesInBoundingBox(2 - eps, -eps, -eps, 2 + eps, eps, eps,
                                        p, 0);
  model::mesh::setSize(p, 0.001);

  // We now identify corresponding surfaces on the left and right sides of the
  // geometry automatically.

  // First we get all surfaces on the left:
  vector<pair<int, int> > sxmin;
  model::getEntitiesInBoundingBox(2 - eps, -eps, -eps, 2 + eps, 1 + eps,
                                        1 + eps, sxmin, 2);
  for(auto i : sxmin) {
    // Then we get the bounding box of each left surface
    double xmin, ymin, zmin, xmax, ymax, zmax;
    model::getBoundingBox(i.first, i.second, xmin, ymin, zmin, xmax, ymax,
                                zmax);
    // We translate the bounding box to the right and look for surfaces inside
    // it:
    vector<pair<int, int> > sxmax;
    model::getEntitiesInBoundingBox(xmin - eps + 1, ymin - eps,
                                          zmin - eps, xmax + eps + 1,
                                          ymax + eps, zmax + eps, sxmax, 2);
    // For all the matches, we compare the corresponding bounding boxes...
    for(auto j : sxmax) {
      double xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
      model::getBoundingBox(j.first, j.second, xmin2, ymin2, zmin2, xmax2,
                                  ymax2, zmax2);
      xmin2 -= 1;
      xmax2 -= 1;
      // ...and if they match, we apply the periodicity constraint
      if(abs(xmin2 - xmin) < eps && abs(xmax2 - xmax) < eps &&
         abs(ymin2 - ymin) < eps && abs(ymax2 - ymax) < eps &&
         abs(zmin2 - zmin) < eps && abs(zmax2 - zmax) < eps) {
        model::mesh::setPeriodic(2, {j.second}, {i.second}, translation);
      }
    }
  }

  model::mesh::generate(3);
  write("t18.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 19
//
//  Thrusections, fillets, pipes, mesh size from curvature
//
// -----------------------------------------------------------------------------

// The OpenCASCADE geometry kernel supports several useful features for solid
// modelling.
void t19(int argc, char **argv)
{
    model::add("t19");

  // Volumes can be constructed from (closed) curve loops thanks to the
  // `addThruSections()' function
  model::occ::addCircle(0, 0, 0, 0.5, 1);
  model::occ::addCurveLoop({1}, 1);
  model::occ::addCircle(0.1, 0.05, 1, 0.1, 2);
  model::occ::addCurveLoop({2}, 2);
  model::occ::addCircle(-0.1, -0.1, 2, 0.3, 3);
  model::occ::addCurveLoop({3}, 3);
  vector<pair<int, int> > out;
  model::occ::addThruSections({1, 2, 3}, out, 1);
  model::occ::synchronize();

  // We can also force the creation of ruled surfaces:
  model::occ::addCircle(2 + 0, 0, 0, 0.5, 11);
  model::occ::addCurveLoop({11}, 11);
  model::occ::addCircle(2 + 0.1, 0.05, 1, 0.1, 12);
  model::occ::addCurveLoop({12}, 12);
  model::occ::addCircle(2 - 0.1, -0.1, 2, 0.3, 13);
  model::occ::addCurveLoop({13}, 13);
  model::occ::addThruSections({11, 12, 13}, out, 11, true, true);
  model::occ::synchronize();

  // We copy the first volume, and fillet all its edges:
  model::occ::copy({{3, 1}}, out);
  model::occ::translate(out, 4, 0, 0);
  model::occ::synchronize();
  vector<pair<int, int> > f;
  model::getBoundary(out, f);
  vector<pair<int, int> > e;
  model::getBoundary(f, e, false);
  vector<int> c;
  for(auto i : e) c.push_back(abs(i.second));
  model::occ::fillet({out[0].second}, c, {0.1}, out);
  model::occ::synchronize();

  // OpenCASCADE also allows general extrusions along a smooth path. Let's first
  // define a spline curve:
  double nturns = 1;
  int npts = 20;
  double r = 1;
  double h = 1 * nturns;
  vector<int> p;
  for(int i = 0; i < npts; i++) {
    double theta = i * 2 * M_PI * nturns / npts;
    model::occ::addPoint(r * cos(theta), r * sin(theta), i * h / npts, 1,
                               1000 + i);
    p.push_back(1000 + i);
  }
  model::occ::addSpline(p, 1000);

  // A wire is like a curve loop, but open:
  model::occ::addWire({1000}, 1000);

  // We define the shape we would like to extrude along the spline (a disk):
  model::occ::addDisk(1, 0, 0, 0.2, 0.2, 1000);
  model::occ::rotate({{2, 1000}}, 0, 0, 0, 1, 0, 0, M_PI / 2);

  // We extrude the disk along the spline to create a pipe (other sweeping types
  // can be specified; try e.g. "Frenet" instead of "DiscreteTrihedron"):
  model::occ::addPipe({{2, 1000}}, 1000, out, "DiscreteTrihedron");

  // We delete the source surface, and increase the number of sub-edges for a
  // nicer display of the geometry:
  model::occ::remove({{2, 1000}});
  option::setNumber("Geometry.NumSubEdges", 1000);

  model::occ::synchronize();

  // We can activate the calculation of mesh element sizes based on curvature
  // (here with a target of 20 elements per 2*Pi radians):
  option::setNumber("Mesh.MeshSizeFromCurvature", 20);

  // We can constraint the min and max element sizes to stay within reasonnable
  // values (see `t10.cpp' for more details):
  option::setNumber("Mesh.MeshSizeMin", 0.001);
  option::setNumber("Mesh.MeshSizeMax", 0.3);

  model::mesh::generate(3);
  write("t19.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 20
//
//  STEP import and manipulation, geometry partitioning
//
// -----------------------------------------------------------------------------

// The OpenCASCADE CAD kernel allows to import STEP files and to modify them. In
// this tutorial we will load a STEP geometry and partition it into slices.
void t20(int argc, char **argv)
{
    model::add("t20");

  // Load a STEP file (using `importShapes' instead of `merge' allows to
  // directly retrieve the tags of the highest dimensional imported entities):
  vector<pair<int, int> > v;
  try {
    model::occ::importShapes("../t20_data.step", v);
  } catch(...) {
    logger::write("Could not load STEP file: bye!");
    finalize();
    return;
  }

  // If we had specified
  //
  // option::setString("Geometry.OCCTargetUnit", "M");
  //
  // before merging the STEP file, OpenCASCADE would have converted the units to
  // meters (instead of the default, which is millimeters).

  // Get the bounding box of the volume:
  double xmin, ymin, zmin, xmax, ymax, zmax;
  model::occ::getBoundingBox(v[0].first, v[0].second, xmin, ymin, zmin,
                                   xmax, ymax, zmax);

  // We want to slice the model into N slices, and either keep the volume slices
  // or just the surfaces obtained by the cutting:

  int N = 5; // Number of slices
  string dir = "X"; // Direction: "X", "Y" or "Z"
  bool surf = false; // Keep only surfaces?

  double dx = (xmax - xmin);
  double dy = (ymax - ymin);
  double dz = (zmax - zmin);
  double L = (dir == "X") ? dz : dx;
  double H = (dir == "Y") ? dz : dy;

  // Create the first cutting plane
  vector<pair<int, int> > s;
  s.push_back({2, model::occ::addRectangle(xmin, ymin, zmin, L, H)});
  if(dir == "X") {
    model::occ::rotate({s[0]}, xmin, ymin, zmin, 0, 1, 0, -M_PI / 2);
  }
  else if(dir == "Y") {
    model::occ::rotate({s[0]}, xmin, ymin, zmin, 1, 0, 0, M_PI / 2);
  }
  double tx = (dir == "X") ? dx / N : 0;
  double ty = (dir == "Y") ? dy / N : 0;
  double tz = (dir == "Z") ? dz / N : 0;
  model::occ::translate({s[0]}, tx, ty, tz);

  // Create the other cutting planes:
  vector<pair<int, int> > tmp;
  for(int i = 1; i < N - 1; i++) {
    model::occ::copy({s[0]}, tmp);
    s.push_back(tmp[0]);
    model::occ::translate({s.back()}, i * tx, i * ty, i * tz);
  }

  // Fragment (i.e. intersect) the volume with all the cutting planes:
  vector<pair<int, int> > ov;
  vector<vector<pair<int, int> > > ovv;
  model::occ::fragment(v, s, ov, ovv);

  // Now remove all the surfaces (and their bounding entities) that are not on
  // the boundary of a volume, i.e. the parts of the cutting planes that "stick
  // out" of the volume:
  model::occ::getEntities(tmp, 2);
  model::occ::remove(tmp, true);

  model::occ::synchronize();

  if(surf) {
    // If we want to only keep the surfaces, retrieve the surfaces in bounding
    // boxes around the cutting planes...
    double eps = 1e-4;
    vector<pair<int, int> > s;
    for(int i = 1; i < N; i++) {
      double xx = (dir == "X") ? xmin : xmax;
      double yy = (dir == "Y") ? ymin : ymax;
      double zz = (dir == "Z") ? zmin : zmax;
      vector<pair<int, int> > e;
      model::getEntitiesInBoundingBox(
        xmin - eps + i * tx, ymin - eps + i * ty, zmin - eps + i * tz,
        xx + eps + i * tx, yy + eps + i * ty, zz + eps + i * tz, e, 2);
      s.insert(s.end(), e.begin(), e.end());
    }
    // ...and remove all the other entities (here directly in the model, as we
    // won't modify any OpenCASCADE entities later on):
    vector<pair<int, int> > dels;
    model::getEntities(dels, 2);
    for(auto it = s.begin(); it != s.end(); ++it) {
      auto it2 = find(dels.begin(), dels.end(), *it);
      if(it2 != dels.end()) dels.erase(it2);
    }
    model::getEntities(tmp, 3);
    model::removeEntities(tmp);
    model::removeEntities(dels);
    model::getEntities(tmp, 1);
    model::removeEntities(tmp);
    model::getEntities(tmp, 0);
    model::removeEntities(tmp);
  }

  // Finally, let's specify a global mesh size and mesh the partitioned model:
  option::setNumber("Mesh.MeshSizeMin", 3);
  option::setNumber("Mesh.MeshSizeMax", 3);
  model::mesh::generate(3);
  write("t20.msh");

  // Launch the GUI to see the results:
  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) fltk::run();
}

// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 21
//
//  Mesh partitioning
//
// -----------------------------------------------------------------------------
// Gmsh can partition meshes using different algorithms, e.g. the graph
// partitioner Metis or the `SimplePartition' plugin. For all the partitioning
// algorithms, the relationship between mesh elements and mesh partitions is
// encoded through the creation of new (discrete) elementary entities, called
// "partition entities".
//
// Partition entities behave exactly like other discrete elementary entities;
// the only difference is that they keep track of both a mesh partition index
// and their parent elementary entity.
//
// The major advantage of this approach is that it allows to maintain a full
// boundary representation of the partition entities, which Gmsh creates
// automatically if `Mesh.PartitionCreateTopology' is set.
void t21(int argc, char **argv)
{
    model::add("t21");

  // Let us start by creating a simple geometry with two adjacent squares
  // sharing an edge:
  model::add("t21");
  model::occ::addRectangle(0, 0, 0, 1, 1, 1);
  model::occ::addRectangle(1, 0, 0, 1, 1, 2);
  vector<pair<int, int> > ov;
  vector<vector<pair<int, int> > > ovv;
  model::occ::fragment({{2, 1}}, {{2, 2}}, ov, ovv);
  model::occ::synchronize();
  vector<pair<int, int> > tmp;
  model::getEntities(tmp, 0);
  model::mesh::setSize(tmp, 0.05);

  // We create one physical group for each square, and we mesh the resulting
  // geometry:
  model::addPhysicalGroup(2, {1}, 100);
  model::setPhysicalName(2, 100, "Left");
  model::addPhysicalGroup(2, {2}, 200);
  model::setPhysicalName(2, 200, "Right");
  model::mesh::generate(2);

  // We now define several ONELAB parameters to fine-tune how the mesh will be
  // partitioned:
  onelab::set(R"( [
  {
    "type":"number",
    "name":"Parameters/0Mesh partitioner",
    "values":[0],
    "choices":[0, 1],
    "valueLabels":{"Metis":0, "SimplePartition":1}
  },
  {
    "type":"number",
    "name":"Parameters/1Number of partitions",
    "values":[3],
    "min":1,
    "max":256,
    "step":1
  },
  {
    "type":"number",
    "name":"Parameters/2Create partition topology (BRep)?",
    "values":[1],
    "choices":[0, 1]
  },
  {
    "type":"number",
    "name":"Parameters/3Create ghost cells?",
    "values":[0],
    "choices":[0, 1]
  },
  {
    "type":"number",
    "name":"Parameters/3Create new physical groups?",
    "values":[0],
    "choices":[0, 1]
  },
  {
    "type":"number",
    "name":"Parameters/3Write file to disk?",
    "values":[1],
    "choices":[0, 1]
  },
  {
    "type":"number",
    "name":"Parameters/4Write one file per partition?",
    "values":[0],
    "choices":[0, 1]
  }
  ] )");

  auto partitionMesh = []()
  {
    // Number of partitions
    vector<double> n;
    onelab::getNumber("Parameters/1Number of partitions", n);
    int N = static_cast<int>(n[0]);

    // Should we create the boundary representation of the partition entities?
    onelab::getNumber("Parameters/2Create partition topology (BRep)?", n);
    option::setNumber("Mesh.PartitionCreateTopology", n[0]);

    // Should we create ghost cells?
    onelab::getNumber("Parameters/3Create ghost cells?", n);
    option::setNumber("Mesh.PartitionCreateGhostCells", n[0]);

    // Should we automatically create new physical groups on the partition
    // entities?
    onelab::getNumber("Parameters/3Create new physical groups?", n);
    option::setNumber("Mesh.PartitionCreatePhysicals", n[0]);

    // Should we keep backward compatibility with pre-Gmsh 4, e.g. to save the
    // mesh in MSH2 format?
    option::setNumber("Mesh.PartitionOldStyleMsh2", 0);

    // Should we save one mesh file per partition?
    onelab::getNumber("Parameters/4Write one file per partition?", n);
    option::setNumber("Mesh.PartitionSplitMeshFiles", n[0]);

    onelab::getNumber("Parameters/0Mesh partitioner", n);
    if(n[0] == 0) {
      // Use Metis to create N partitions
      model::mesh::partition(N);
      // Several options can be set to control Metis: `Mesh.MetisAlgorithm' (1:
      // Recursive, 2: K-way), `Mesh.MetisObjective' (1: min. edge-cut, 2:
      // min. communication volume), `Mesh.PartitionTriWeight' (weight of
      // triangles), `Mesh.PartitionQuadWeight' (weight of quads), ...
    }
    else {
      // Use the `SimplePartition' plugin to create chessboard-like partitions
      plugin::setNumber("SimplePartition", "NumSlicesX", N);
      plugin::setNumber("SimplePartition", "NumSlicesY", 1);
      plugin::setNumber("SimplePartition", "NumSlicesZ", 1);
      plugin::run("SimplePartition");
    }

    // Save mesh file (or files, if `Mesh.PartitionSplitMeshFiles' is set):
    onelab::getNumber("Parameters/3Write file to disk?", n);
    if(n[0]) write("t21.msh");

    // Iterate over partitioned entities and print some info (see the first
    // extended tutorial `x1.cpp' for additional information):
    vector<pair<int, int> > entities;
    model::getEntities(entities);

    for(auto e : entities) {
      vector<int> partitions;
      model::getPartitions(e.first, e.second, partitions);
      if(partitions.size()) {
        string type;
        model::getType(e.first, e.second, type);
        cout << "Entity (" << e.first << "," << e.second << ") "
                  << "of type " << type << "\n";
        cout << " - Partition(s):";
        for(auto p : partitions) cout << " " << p;
        cout << "\n";
        int pdim, ptag;
        model::getParent(e.first, e.second, pdim, ptag);
        cout << " - Parent: (" << pdim << "," << ptag << ")\n";
        vector<pair<int, int> > bnd;
        model::getBoundary({e}, bnd);
        cout << " - Boundary:";
        for(auto b : bnd) cout << " (" << b.first << "," << b.second << ")";
        cout << "\n";
      }
    }
  };

  partitionMesh();

  // Launch the GUI and handle the "check" event to re-partition the mesh
  // according to the choices made in the GUI
  auto checkForEvent = [=]() -> bool {
    vector<string> action;
    onelab::getString("ONELAB/Action", action);
    if(action.size() && action[0] == "check") {
      onelab::setString("ONELAB/Action", {""});
      partitionMesh();
      graphics::draw();
    }
    return true;
  };

  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) {
    fltk::initialize();
    while(fltk::isAvailable() && checkForEvent())
      fltk::wait();
  }
}

void test_group(int argc, char **argv)
{
    open(argv[1]);
    model::mesh::generate(2);
    map<string, pair<int, int>> groups;
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    string groupname = "left_edge";

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return;
    }

    // get entities of the chosen physical group
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            cout << "\telement type " << elementTypes[i] << ":\n";
            cout << "\t\telements: ";
            for (int j = 0; j < elementTags[i].size(); ++j)
                cout << elementTags[i][j] << " ";
            cout << '\n';
            cout << "\t\tnodes: ";
            for (int j = 0; j < nodeTags[i].size(); ++j)
                cout << nodeTags[i][j] << " ";
            cout << '\n';
        }
    }
}

void functions_test(int argc, char **argv)
{
    open(argv[1]);
    model::mesh::generate(2);
    //code numero 4
    map<string, pair<int, int>> groups;
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    string groupname = "left_edge";

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return;
    }

    // get entities of the chosen physical group
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            cout << "\telement type " << elementTypes[i] << ":\n";
            cout << "\t\telements: ";
            for (int j = 0; j < elementTags[i].size(); ++j)
                cout << elementTags[i][j] << " ";
            cout << '\n';
            cout << "\t\tnodes: ";
            for (int j = 0; j < nodeTags[i].size(); ++j)
                cout << nodeTags[i][j] << " ";
            cout << '\n';
        }
    }

//code numero 3

    model::getPhysicalGroups(dimTags);

    for(int i=0; i<dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;

        string name;
        model::getPhysicalName(dim, tag, name);

        cout << "Physical group (" << dim << "D) "
                  << "named '"<< name << "' tag=" << tag << '\n';

        // -- code added to previous code (begin) ------------------------------
        // set value (dim, tag) for key 'name'
        groups[name] = {dim, tag};
        // -- code added to previous code (end) --------------------------------
    }

    // -- code added to previous code (begin) ----------------------------------
    // print map
    cout << "Physical Groups:\n";
    for (auto &x : groups)
        cout << "\t '" << x.first << "' => (" << x.second.first << ',' 
                                                   << x.second.second << ")\n";
    // -- code added to previous code (end) ------------------------------------

// code numero 2
    model::getPhysicalGroups(dimTags); // vector is filled

    for(int i=0; i<dimTags.size(); ++i) // size gives size of the vector
    {
        int dim = dimTags[i].first;   // dimTags[i] is a pair<int,int>
        int tag = dimTags[i].second;

        string name; // output parameter of "getPhysicalName"
        model::getPhysicalName(dim, tag, name);

        // prints the data to the console
        cout << "Physical group (" << dim << "D) "
                  << "named '"<< name << "' tag=" << tag << '\n';
    }

//code numero 5
    // input parameter
    int elementType = 1;

    // outputs
    string elementName;
    int order, numNodes, numPrimaryNodes;
    vector<double> localNodeCoord;

    // gmsh call
    model::mesh::getElementProperties(elementType,
                                            elementName, dim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);

    // display data returned by gmsh
    cout << "Element type " << elementType << ":\n";
    cout << "\tname = '" << elementName << "'\n";
    cout << "\tdim = " << dim << '\n';
    cout << "\torder = " << order << '\n';
    cout << "\tnumNodes = " << numNodes << '\n';
    cout << "\tlocalNodeCoord = ";
    for(int i=0; i<localNodeCoord.size(); ++i)
        cout << localNodeCoord[i] << " ";
    cout << '\n';
    cout << "\tnumPrimaryNodes = " << numPrimaryNodes << '\n';

//code numero 6
// Renumber the node tags in a continuous sequence.
    model::mesh::renumberNodes(); 

    // output parameters
    vector<size_t> nodeTags;
    vector<double> coord;
    vector<double> parametricCoord;

    // gmsh call
    model::mesh::getNodes(nodeTags, coord, parametricCoord);

    // display the vectors
    cout << "\tnodeTags = ";
    for(int i=0; i<nodeTags.size(); ++i)
        cout << nodeTags[i] << " ";
    cout << '\n';

    cout << "\tcoord = ";
    for(int i=0; i<coord.size(); ++i)
        cout << coord[i] << " ";
    cout << '\n';   

    // code numero 7
    // get selected keys from the onelab database
    vector<string> keys;
    onelab::getNames(keys, "(Boundary Conditions|Materials).+"); //see rectangle.geo file, retrieves boundary conditions and mat
    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);

        // print received data
        cout << key << ":";
        for (auto &v : value)
            cout << " " << v;
        cout << '\n';

        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            cout << "\t.type of data = '" << words[0] << "'\n";
            cout << "\t.physical group name = '" << words[1] << "'\n";
            cout << "\t.field = '" << words[2] << "'\n";
            cout << "\t.value = " << value[0] << '\n';
        }
    }

//code numero 8
for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    groupname = "domain";  // <= we choose a physical group here!

    dim = groups[groupname].first;
    tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return;
    }

    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            // -- CODE ADDED TO PREVIOUS EXAMPLE (begin) -----------------------

            // get element name as a string
            string elementName;
            int edim, eorder, enumNodes, enumPrimaryNodes;
            vector<double> localNodeCoord;
            model::mesh::getElementProperties(elementTypes[i],
                                                    elementName, dim, eorder,
                                                    enumNodes, localNodeCoord, 
                                                    enumPrimaryNodes);
            cout << '\t' << elementName << " order=" << eorder 
                      << " nodes=" << enumNodes << '\n';

            // get Gauss points coordinates and weights
            vector<double> localCoords, weights;
            model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    "CompositeGauss2", // "Gauss2"
                                                    localCoords, weights);
            cout << "\t" << weights.size() << " integration points\n";

            cout << "\tweights: ";
            for (int j = 0; j < weights.size(); ++j)
                cout << weights[j] << " ";
            cout << '\n';

            cout << "\tlocalCoords (gmsh): ";
            for (int j = 0; j < localCoords.size(); ++j)
                cout << localCoords[j] << " ";
            cout << '\n';

            // convert to Eigen format
            Map<Matrix<double, Dynamic, Dynamic>>
                localCoords_E(&localCoords[0], localCoords.size()/3, 3);
            cout << "\tlocalCoords (Eigen):\n";
            IOFormat fmt(4, 0, ", ", "\n", "\t\t[", "]");
            cout << localCoords_E.format(fmt) << '\n';         

            // get Jacobians
            vector<double> jacobians, determinants, coords;
            model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords);
            cout << "\tGot " << determinants.size() 
                    << " Jacobians for type " << elementTypes[i] << "\n";

            cout << "\tjacobians (gmsh): ";
            for (int j = 0; j < jacobians.size(); ++j)
                cout << jacobians[j] << " ";
            cout << '\n';

            // convert first jacobian to Eigen format
            Map<Matrix<double, Dynamic, Dynamic> >
                jacob(&jacobians[0], 3, 3);
            cout << "\tfirst jacobian (Eigen):\n" 
                      << jacob.format(fmt) << '\n';

            // compute the inverse with Eigen
            Matrix<double, 3, 3> jacobinv = jacob.inverse();
            cout << "\tinverse of the first jacobian (Eigen):\n" 
                      << jacobinv.format(fmt) << "\n";

            // Compute J*Jinv with Eigen "*" operator
            cout << "\tVerification (J*Jinv):\n" 
                      << (jacob*jacobinv).format(fmt) << '\n';

            // shape functions and derivatives at the Gauss points
            int numComponents, numOrientations;
            vector<double> basisFunctions;
            model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);
            cout << "\tGot " << basisFunctions.size() 
                      << " basis functions for type "
                      << elementTypes[i] << '\n';

            model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);
            cout << "\tGot " << basisFunctions.size() 
                      << " basis function gradients for type "
                      << elementTypes[i] << '\n';

            // -- CODE ADDED TO PREVIOUS EXAMPLE (end) -------------------------
        }
    }

    //code numero 9
    vector<string> names;
    model::list(names);
    cout << "the file contains " << names.size() << " model names\n";

    // get all the nodes
    vector<double> coords;
    model::mesh::getNodes(nodeTags, coords, parametricCoord);

    cout << "the mesh has " << nodeTags.size() << " nodes\n";

    // get all the elements for dim=2
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elnodeTags;
    model::mesh::getElements(elementTypes, elementTags, elnodeTags, 2);
    
    // compute total number of 2D elements
    int nbelems=0;
    for(size_t i=0; i< elementTags.size(); ++i)
        nbelems+=elementTags[i].size();
    cout << "the mesh has " << nbelems << " elements\n";

    // Create a new post-processing view
    int viewtag1 = view::add("NodeData");
    int viewtag2 = view::add("ElementNodeData");

    size_t nstep = 20;
    double Lx = 20.;
    double v = 10.0;
    double tend = Lx / v;
    for (size_t step = 0; step < nstep; step++)
    {
        double time = tend * ((double)step / nstep);

        // ** NodeData
        // loop over the nodes and compute some values from their coordinates
        vector<vector<double>> data1(nodeTags.size()); // check this size!!
        for (size_t i = 0; i < nodeTags.size(); ++i)
        {
            size_t tag = nodeTags[i];
            double x = coords[i * 3 + 0];
            double y = coords[i * 3 + 1];
            double z = coords[i * 3 + 2];

            double val = sin(2 * M_PI * (x - v * time) / Lx * 2) + y;
            data1[i].resize(1);
            data1[i][0] = val;
        }
        view::addModelData(viewtag1, step, names[0], "NodeData",
                                 nodeTags, data1, time);

        // ** ElementNodeData
        // loop over the elements and compute some values at the coords of its nodes

        vector<vector<double>> data2(nbelems);
        vector<size_t> elems2D(nbelems);

        // for each element type
        //cout << "there are " << elementTypes.size() << " element type(s)\n";
        int k=0;
        for (size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
        {
            // get info about this type of element
            string name;
            int dim, order, numNodes, numPrimaryNodes;
            vector<double> paramCoord;
            model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                                    numNodes, paramCoord, numPrimaryNodes);
            //cout << "element type has " << numNodes << " nodes\n";

            // select element/node tags for this type of element
            auto &etags = elementTags[ityp];
            auto &enods = elnodeTags[ityp];

            // loop over elements of type "ityp"
            //cout << "loop over the " << etags.size() << " elements\n";
            for (size_t i = 0; i < etags.size(); ++i)
            {
                elems2D[k] = etags[i];

                double offset = (elems2D[k]%2) *0.5; // so that fields are discontinious across elements

                data2[k].resize(numNodes);
                // loop over nodes of the element
                for (size_t j = 0; j < numNodes; ++j)
                {
                    int nodeTag = enods[numNodes * i + j];

                    vector<double> nodecoord;
                    vector<double> nodepcoord;
                    int dim, tag;
                    model::mesh::getNode(nodeTag, nodecoord, nodepcoord, dim, tag);
                    double x = nodecoord[0];
                    double y = nodecoord[1];
                    double z = nodecoord[2];
                    double val = sin(2 * M_PI * (x - v * time) / Lx * 2) + y;
                    data2[k][j] = val + offset;
                }
                k++;
            }
        }
          view::addModelData(viewtag2, step, names[0], "ElementNodeData",
                                   elems2D, data2, time, 1); // the last ,1 is important!
      
    }

    // save the results to disk
    view::write(viewtag1, "results.msh");
    view::write(viewtag2, "results.msh", true);

    // set display options & view results
    option::setNumber("View[0].IntervalsType", 3);
    option::setNumber("View[0].AdaptVisualizationGrid", 1);
    option::setNumber("View[0].MaxRecursionLevel", 3);
    option::setNumber("View[0].TargetError", -0.0001);
    
    option::setNumber("View[1].IntervalsType", 3);
    option::setNumber("View[1].AdaptVisualizationGrid", 1);
    option::setNumber("View[1].MaxRecursionLevel", 3);
    option::setNumber("View[1].TargetError", -0.0001);

    fltk::run();

//code numero 1

    fltk::run(); // opens the gmsh window
    
    // -- code added to previous code (end) ------------------------------------
}

void sphere_avec_trous(int argc, char **argv)
{
    model::add (" boolean ");
    double R = 1.4, Rs = R *.7, Rt = R *1.25;
    vector<pair<int, int> > ov; 
    vector<vector<pair<int, int> > > ovv;
    model::occ::addBox( -R , -R , -R , 2* R ,2* R ,2* R , 1);
    model::occ::addSphere(0 ,0 ,0 , Rt , 2);
    model::occ::intersect({{3 , 1}} , {{3 , 2}}, ov, ovv, 3);
    model::occ::addCylinder( -2* R ,0 ,0 , 4* R ,0 ,0 , Rs , 4);
    model::occ::addCylinder (0 , -2* R ,0 , 0 ,4* R ,0 , Rs , 5);
    model::occ::addCylinder (0 ,0 , -2* R , 0 ,0 ,4* R , Rs , 6);
    model::occ::fuse({{3 , 4} , {3 , 5}} , {{3 , 6}}, ov, ovv, 7);
    model::occ::cut({{3 , 3}} , {{3 , 7}}, ov, ovv, 8);
    model::occ::synchronize();
    model::mesh::generate(3);
    fltk::run();
}

void structure_avec_torsion(int argc, char **argv)
{
    auto createGeometryAndMesh = []()
  {
    // Clear all models and create a new one
    clear();
    model::add("t3");

    // Copied from t1.cpp...
    double lc = 1e-2;
    model::geo::addPoint(0, 0, 0, lc, 1);
    model::geo::addPoint(.1, 0, 0, lc, 2);
    model::geo::addPoint(.1, .3, 0, lc, 3);
    model::geo::addPoint(0, .3, 0, lc, 4);
    model::geo::addLine(1, 2, 1);
    model::geo::addLine(3, 2, 2);
    model::geo::addLine(3, 4, 3);
    model::geo::addLine(4, 1, 4);
    model::geo::addCurveLoop({4, 1, -2, 3}, 1);
    model::geo::addPlaneSurface({1}, 1);
    model::geo::synchronize();
    model::addPhysicalGroup(1, {1, 2, 4}, 5);
    int ps = model::addPhysicalGroup(2, {1});
    model::setPhysicalName(2, ps, "My surface");

    // As in `t2.cpp', we plan to perform an extrusion along the z axis.  But
    // here, instead of only extruding the geometry, we also want to extrude the
    // 2D mesh. This is done with the same `extrude()' function, but by
    // specifying element 'Layers' (2 layers in this case, the first one with 8
    // subdivisions and the second one with 2 subdivisions, both with a height
    // of h/2). The number of elements for each layer and the (end) height of
    // each layer are specified in two vectors:
    double h = 0.1;
    vector<pair<int, int> > ov;
    model::geo::extrude({{2, 1}}, 0, 0, h, ov, {8, 2}, {0.5, 1});

    // The extrusion can also be performed with a rotation instead of a
    // translation, and the resulting mesh can be recombined into prisms (we use
    // only one layer here, with 7 subdivisions). All rotations are specified by
    // an an axis point (-0.1, 0, 0.1), an axis direction (0, 1, 0), and a
    // rotation angle (-Pi/2):
    model::geo::revolve({{2, 28}}, -0.1, 0, 0.1, 0, 1, 0, -M_PI / 2, ov,
                              {7});

    // Using the built-in geometry kernel, only rotations with angles < Pi are
    // supported. To do a full turn, you will thus need to apply at least 3
    // rotations. The OpenCASCADE geometry kernel does not have this limitation.

    // A translation (-2 * h, 0, 0) and a rotation ((0, 0.15, 0.25), (1, 0, 0),
    // angle * Pi / 180) can also be combined to form a "twist". The last
    // (optional) argument for the extrude() and twist() functions specifies
    // whether the extruded mesh should be recombined or not. The `angle'
    // parameter is retrieved from the ONELAB database (it can be set
    // interactively in the GUI -- see below):
    vector<double> angle;
    onelab::getNumber("Parameters/Twisting angle", angle);
    model::geo::twist({{2, 50}}, 0, 0.15, 0.25, -2 * h, 0, 0, 1, 0, 0,
                            angle[0] * M_PI / 180., ov, {10}, {}, true);

    model::geo::synchronize();

    // All the extrusion functions return a vector of extruded entities: the
    // "top" of the extruded surface (in `ov[0]'), the newly created volume (in
    // `ov[1]') and the tags of the lateral surfaces (in `ov[2]', `ov[3]', ...).

    // We can then define a new physical volume (with tag 101) to group all the
    // elementary volumes:
    model::addPhysicalGroup(3, {1, 2, ov[1].second}, 101);

    model::mesh::generate(3);
    write("t3.msh");
  };

  // Let us now change some options... Since all interactive options are
  // accessible through the API, we can for example make point tags visible or
  // redefine some colors:
  option::setNumber("Geometry.PointNumbers", 1);
  option::setColor("Geometry.Points", 255, 165, 0);
  option::setColor("General.Text", 255, 255, 255);
  option::setColor("Mesh.Points", 255, 0, 0);

  // Note that color options are special: setting a color option of "X.Y"
  // actually sets the option "X.Color.Y".
  int r, g, b, a;
  option::getColor("Geometry.Points", r, g, b, a);
  option::setColor("Geometry.Surfaces", r, g, b, a);

  // We create a ONELAB parameter to define the angle of the twist. ONELAB
  // parameters can be modified interactively in the GUI, and can be exchanged
  // with other codes connected to the same ONELAB database. The database can be
  // accessed through the Gmsh C++ API using JSON-formatted strings (see
  // https://gitlab.onelab.info/doc/tutorials/-/wikis/ONELAB-JSON-interface for
  // more information):
  onelab::set(R"( [
  {
    "type":"number",
    "name":"Parameters/Twisting angle",
    "values":[90],
    "min":0,
    "max":120,
    "step":1
  }
  ] )");

  // Create the geometry and mesh it:
  createGeometryAndMesh();

  // Launch the GUI and handle the "check" event (recorded in the
  // "ONELAB/Action" parameter) to recreate the geometry with a new twisting
  // angle if necessary:
  auto checkForEvent = [=]() -> bool {
    vector<string> action;
    onelab::getString("ONELAB/Action", action);
    if(action.size() and action[0] == "check") {
      onelab::setString("ONELAB/Action", {""});
      createGeometryAndMesh();
      graphics::draw();
    }
    return true;
  };

  set<string> args(argv, argv + argc);
  if(!args.count("-nopopup")) {
    fltk::initialize();
    while(fltk::isAvailable() && checkForEvent())
      fltk::wait();
  }

  // When the GUI is launched, you can use the `Help->Current Options and
  // Workspace' menu to see the current values of all options. To save the
  // options in a file, use `File->Export->Gmsh Options', or through the api:

  // write("t3.opt");

}

// This program takes a .geo file as argument, generates the mesh of the 
// surfaces and displays the result in the Gmsh window.
// 
// This example introduces:
// - the gmsh.h header
// - initialize()/finalize() functions
// - C++ namespaces / standard library (cout)
//
// How to build/run? (windows) 
//   [in examples/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code1_loadgeo.exe ..\rectangle.geo
void code1_loadgeo(int argc, char **argv)
{
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]); // similar to "File / Open" in the GUI
    model::mesh::generate(2);  // similar to "Mesh / 2D" in the GUI

    fltk::run(); // opens the gmsh window
}

// This program loops over the physical groups defined in a .geo file
// and prints their names, dimensions and tags.
//
// This example introduces:
// - (dim,tags)
// - pair<T1,T2>
// - vector<T>
// - vectorpair == vector<pair<int,int>>
// - references / const references
// - how to call a gmsh function (input vs output arguments)
// - string
//
// How to build/run? (windows) 
//   [in examples/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code2_getgroups.exe ..\rectangle.geo
void code2_getgroups(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]);
    model::mesh::generate(2);

    // -- code added to code1.cpp (begin) --------------------------------------

    vectorpair dimTags; // output parameter of "getPhysicalGroups", defining an empty vector of pairs
    model::getPhysicalGroups(dimTags); // vector is filled

    for(int i=0; i<dimTags.size(); ++i) // size gives size of the vector
    {
        int dim = dimTags[i].first;   // dimTags[i] is a pair<int,int>
        int tag = dimTags[i].second;

        string name; // output parameter of "getPhysicalName"
        model::getPhysicalName(dim, tag, name);

        // prints the data to the console
        cout << "Physical group (" << dim << "D) "
                  << "named '"<< name << "' tag=" << tag << '\n';
    }

    // -- code added to code1.cpp (end) ----------------------------------------
}

// This program loops over the physical groups defined in a .geo file
// and prints their names, dimensions and tags.
//
// This example introduces:
// - map
// - Physical Group 
//      => Entities 
//          => Element types 
//              => Elements/nodes hierarchy
//
// How to build/run? (windows) 
//   [in example/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code3_buildphymap.exe ..\rectangle.geo
void code3_buildphysmap(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }
    open(argv[1]);
    model::mesh::generate(2);

    // -- code added to previous code (begin) ----------------------------------
    // definition of the map: "group_name" => (dim,tag)
    map< string, pair<int,int> > groups;
    // -- code added to previous code (end) ------------------------------------

    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);

    for(int i=0; i<dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;

        string name;
        model::getPhysicalName(dim, tag, name);

        cout << "Physical group (" << dim << "D) "
                  << "named '"<< name << "' tag=" << tag << '\n';

        // -- code added to previous code (begin) ------------------------------
        // set value (dim, tag) for key 'name'
        groups[name] = {dim, tag};
        // -- code added to previous code (end) --------------------------------
    }

    // -- code added to previous code (begin) ----------------------------------
    // print map
    cout << "Physical Groups:\n";
    for (auto &x : groups)
        cout << "\t '" << x.first << "' => (" << x.second.first << ',' 
                                                   << x.second.second << ")\n";
    // -- code added to previous code (end) ------------------------------------
}

// This program loops over the physical groups defined in a .geo file
// and prints their names, dimensions and tags.
//
// This example loops over the gmsh data structure hierarchy:
// - Physical Group
//      => contains 1 or more "Entities"
//          => contains 1 or more "Element types"
//              => access to Elements
//
// How to build/run? (windows) 
//   [in example/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code4_hierarchy.exe ..\rectangle.geo
void code4_hierarchy(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]);
    model::mesh::generate(2);

    map<string, pair<int, int>> groups;
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    // -- code added to previous code (begin) ----------------------------------

    string groupname = "left_edge";

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return;
    }

    // get entities of the chosen physical group
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            cout << "\telement type " << elementTypes[i] << ":\n";
            cout << "\t\telements: ";
            for (int j = 0; j < elementTags[i].size(); ++j)
                cout << elementTags[i][j] << " ";
            cout << '\n';
            cout << "\t\tnodes: ";
            for (int j = 0; j < nodeTags[i].size(); ++j)
                cout << nodeTags[i][j] << " ";
            cout << '\n';
        }
    }
    // -- code added to previous code (end) ------------------------------------
}

// This program displays info about "element type 1".
//
// How to build/run? (windows) 
//   [in example/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code5_elprops.exe
void code5_elprops(int argc, char **argv)
{
    // input parameter
    int elementType = 1;

    // outputs
    string elementName;
    int dim, order, numNodes, numPrimaryNodes;
    vector<double> localNodeCoord;

    // gmsh call
    model::mesh::getElementProperties(elementType,
                                            elementName, dim, order,
                                            numNodes, localNodeCoord, 
                                            numPrimaryNodes);

    // display data returned by gmsh
    cout << "Element type " << elementType << ":\n";
    cout << "\tname = '" << elementName << "'\n";
    cout << "\tdim = " << dim << '\n';
    cout << "\torder = " << order << '\n';
    cout << "\tnumNodes = " << numNodes << '\n';
    cout << "\tlocalNodeCoord = ";
    for(int i=0; i<localNodeCoord.size(); ++i)
        cout << localNodeCoord[i] << " ";
    cout << '\n';
    cout << "\tnumPrimaryNodes = " << numPrimaryNodes << '\n';
}

// This program gets all the node tags and the coordinates of the model.
//
// How to build/run? (windows) 
//   [in example/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code6_allnodes.exe ..\rectangle.geo
void code6_allnodes(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]);
    model::mesh::generate(2);

    // Renumber the node tags in a continuous sequence.
    model::mesh::renumberNodes(); 

    // output parameters
    vector<size_t> nodeTags;
    vector<double> coord;
    vector<double> parametricCoord;

    // gmsh call
    model::mesh::getNodes(nodeTags, coord, parametricCoord);


    // display the vectors
    cout << "\tnodeTags = ";
    for(int i=0; i<nodeTags.size(); ++i)
        cout << nodeTags[i] << " ";
    cout << '\n';

    cout << "\tcoord = ";
    for(int i=0; i<coord.size(); ++i)
        cout << coord[i] << " ";
    cout << '\n';   
}

// This program takes a .geo file as argument, and reads the onelab parameters
// defined by the "SetNumber" command.
// This technique can be used to define additional parameters for the solver
// such as the values for the boundary conditions or material parameters.
//
// How to build/run? (windows)
//   [in example/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code7_parameters.exe ..\rectangle.geo
void code7_parameters(int argc, char **argv)
{
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]);

    // get selected keys from the onelab database
    vector<string> keys;
    onelab::getNames(keys, "(Boundary Conditions|Materials).+"); //see rectangle.geo file, retrieves boundary conditions and mat
    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);

        // print received data
        cout << key << ":";
        for (auto &v : value)
            cout << " " << v;
        cout << '\n';

        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            cout << "\t.type of data = '" << words[0] << "'\n";
            cout << "\t.physical group name = '" << words[1] << "'\n";
            cout << "\t.field = '" << words[2] << "'\n";
            cout << "\t.value = " << value[0] << '\n';
        }
    }
}

// This program demonstrates how to convert gmsh vectors to Eigen matrices.
// It also shows:
// - how to select a quadrature type and get the positions and weights of 
//   the Gauss points.
// - how to compute the inverse of the Jacobians.
// - how to get the values of the shape functions and their derivatives at
//   the Gauss points.
//
// How to build/run? (windows)
//   [in example/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code8_jacobians.exe ..\mono.geo
void code8_jacobians(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]);
    model::mesh::generate(2);

    map<string, pair<int, int>> groups;
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    string groupname = "domain";  // <= we choose a physical group here!

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return;
    }

    // get entities of the chosen physical group
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            // -- CODE ADDED TO PREVIOUS EXAMPLE (begin) -----------------------

            // get element name as a string
            string elementName;
            int edim, eorder, enumNodes, enumPrimaryNodes;
            vector<double> localNodeCoord;
            model::mesh::getElementProperties(elementTypes[i],
                                                    elementName, dim, eorder,
                                                    enumNodes, localNodeCoord, 
                                                    enumPrimaryNodes);
            cout << '\t' << elementName << " order=" << eorder 
                      << " nodes=" << enumNodes << '\n';

            // get Gauss points coordinates and weights
            vector<double> localCoords, weights;
            model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    "CompositeGauss2", // "Gauss2"
                                                    localCoords, weights);
            cout << "\t" << weights.size() << " integration points\n";

            cout << "\tweights: ";
            for (int j = 0; j < weights.size(); ++j)
                cout << weights[j] << " ";
            cout << '\n';

            cout << "\tlocalCoords (gmsh): ";
            for (int j = 0; j < localCoords.size(); ++j)
                cout << localCoords[j] << " ";
            cout << '\n';

            // convert to Eigen format
            Map<Matrix<double, Dynamic, Dynamic>>
                localCoords_E(&localCoords[0], localCoords.size()/3, 3);
            cout << "\tlocalCoords (Eigen):\n";
            IOFormat fmt(4, 0, ", ", "\n", "\t\t[", "]");
            cout << localCoords_E.format(fmt) << '\n';         

            // get Jacobians
            vector<double> jacobians, determinants, coords;
            model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords);
            cout << "\tGot " << determinants.size() 
                    << " Jacobians for type " << elementTypes[i] << "\n";

            cout << "\tjacobians (gmsh): ";
            for (int j = 0; j < jacobians.size(); ++j)
                cout << jacobians[j] << " ";
            cout << '\n';

            // convert first jacobian to Eigen format
            Map<Matrix<double, Dynamic, Dynamic> >
                jacob(&jacobians[0], 3, 3);
            cout << "\tfirst jacobian (Eigen):\n" 
                      << jacob.format(fmt) << '\n';

            // compute the inverse with Eigen
            Matrix<double, 3, 3> jacobinv = jacob.inverse();
            cout << "\tinverse of the first jacobian (Eigen):\n" 
                      << jacobinv.format(fmt) << "\n";

            // Compute J*Jinv with Eigen "*" operator
            cout << "\tVerification (J*Jinv):\n" 
                      << (jacob*jacobinv).format(fmt) << '\n';

            // shape functions and derivatives at the Gauss points
            int numComponents, numOrientations;
            vector<double> basisFunctions;
            model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "Lagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);
            cout << "\tGot " << basisFunctions.size() 
                      << " basis functions for type "
                      << elementTypes[i] << '\n';

            model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);
            cout << "\tGot " << basisFunctions.size() 
                      << " basis function gradients for type "
                      << elementTypes[i] << '\n';

            // -- CODE ADDED TO PREVIOUS EXAMPLE (end) -------------------------
        }
    }
}

// This example generate some results over a mesh
//
// How to use this example?
//
// 1.build the code:
//     cd build
//     cmake ..
//     make
// 2.generate a mesh from a geo file:
//     gmsh -2 -order 3 ..\src\mymesh.geo
// 3.run the program with the msh as argument
//     bin\myview.exe ..\src\mymesh.msh
// 4.display the generated msh data
//     gmsh results.msh
void code9_views(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }

    open(argv[1]);
    model::mesh::generate(2); 

    vector<string> names;
    model::list(names);
    cout << "the file contains " << names.size() << " model names\n";

    // get all the nodes
    vector<size_t> nodeTags;
    vector<double> coords;
    vector<double> parametricCoord;
    model::mesh::getNodes(nodeTags, coords, parametricCoord);

    cout << "the mesh has " << nodeTags.size() << " nodes\n";

    // get all the elements for dim=2
    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elnodeTags;
    model::mesh::getElements(elementTypes, elementTags, elnodeTags, 2);
    
    // compute total number of 2D elements
    int nbelems=0;
    for(size_t i=0; i< elementTags.size(); ++i)
        nbelems+=elementTags[i].size();
    cout << "the mesh has " << nbelems << " elements\n";

    // Create a new post-processing view
    int viewtag1 = view::add("NodeData");
    int viewtag2 = view::add("ElementNodeData");

    size_t nstep = 20;
    double Lx = 20.;
    double v = 10.0;
    double tend = Lx / v;
    for (size_t step = 0; step < nstep; step++)
    {
        double time = tend * ((double)step / nstep);

        // ** NodeData
        // loop over the nodes and compute some values from their coordinates
        vector<vector<double>> data1(nodeTags.size()); // check this size!!
        for (size_t i = 0; i < nodeTags.size(); ++i)
        {
            size_t tag = nodeTags[i];
            double x = coords[i * 3 + 0];
            double y = coords[i * 3 + 1];
            double z = coords[i * 3 + 2];

            double val = sin(2 * M_PI * (x - v * time) / Lx * 2) + y;
            data1[i].resize(1);
            data1[i][0] = val;
        }
        view::addModelData(viewtag1, step, names[0], "NodeData",
                                 nodeTags, data1, time);

        // ** ElementNodeData
        // loop over the elements and compute some values at the coords of its nodes

        vector<vector<double>> data2(nbelems);
        vector<size_t> elems2D(nbelems);

        // for each element type
        //cout << "there are " << elementTypes.size() << " element type(s)\n";
        int k=0;
        for (size_t ityp = 0; ityp < elementTypes.size(); ++ityp)
        {
            // get info about this type of element
            string name;
            int dim, order, numNodes, numPrimaryNodes;
            vector<double> paramCoord;
            model::mesh::getElementProperties(elementTypes[ityp], name, dim, order,
                                                    numNodes, paramCoord, numPrimaryNodes);
            //cout << "element type has " << numNodes << " nodes\n";

            // select element/node tags for this type of element
            auto &etags = elementTags[ityp];
            auto &enods = elnodeTags[ityp];

            // loop over elements of type "ityp"
            //cout << "loop over the " << etags.size() << " elements\n";
            for (size_t i = 0; i < etags.size(); ++i)
            {
                elems2D[k] = etags[i];

                double offset = (elems2D[k]%2) *0.5; // so that fields are discontinious across elements

                data2[k].resize(numNodes);
                // loop over nodes of the element
                for (size_t j = 0; j < numNodes; ++j)
                {
                    int nodeTag = enods[numNodes * i + j];

                    vector<double> nodecoord;
                    vector<double> nodepcoord;
                    int dim, tag;
                    model::mesh::getNode(nodeTag, nodecoord, nodepcoord, dim, tag);
                    double x = nodecoord[0];
                    double y = nodecoord[1];
                    double z = nodecoord[2];
/*
                    double x = coords[nod * 3 + 0];
                    double y = coords[nod * 3 + 1];
                    double z = coords[nod * 3 + 2];
*/
                    double val = sin(2 * M_PI * (x - v * time) / Lx * 2) + y;
                    data2[k][j] = val + offset;
                }
                k++;
            }
        }
          view::addModelData(viewtag2, step, names[0], "ElementNodeData",
                                   elems2D, data2, time, 1); // the last ,1 is important!
      
    }

    // save the results to disk
    view::write(viewtag1, "results.msh");
    view::write(viewtag2, "results.msh", true);

    // set display options & view results
    option::setNumber("View[0].IntervalsType", 3);
    option::setNumber("View[0].AdaptVisualizationGrid", 1);
    option::setNumber("View[0].MaxRecursionLevel", 3);
    option::setNumber("View[0].TargetError", -0.0001);
    
    option::setNumber("View[1].IntervalsType", 3);
    option::setNumber("View[1].AdaptVisualizationGrid", 1);
    option::setNumber("View[1].MaxRecursionLevel", 3);
    option::setNumber("View[1].TargetError", -0.0001);

    fltk::run();
}

// This program takes a .geo file as argument, generates the mesh of the 
// surfaces and displays the result in the Gmsh window.
// 
// This example introduces:
// - the gmsh.h header
// - initialize()/finalize() functions
// - C++ namespaces / standard library (cout)
//
// How to build/run? (windows) 
//   [in examples/gmsh_api]
//   mkdir build
//   cd build
//   cmake .. && make && code1_loadgeo.exe ..\rectangle.geo
void mycode(int argc, char **argv)
{
    /*--------------INIT----------------*/
    // the program requires 1 argument: a .geo file.
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <geo_file>\n";
        return;
    }
    open(argv[1]); // similar to "File / Open" in the GUI
    model::mesh::generate(2);  // similar to "Mesh / 2D" in the GUI

    /*--------get physical groups (more particularly: the domain)--------------*/
    map<string, pair<int, int>> groups;
    vectorpair dimTags;
    model::getPhysicalGroups(dimTags);
    for (int i = 0; i < dimTags.size(); ++i)
    {
        int dim = dimTags[i].first;
        int tag = dimTags[i].second;
        string name;
        model::getPhysicalName(dim, tag, name);
        groups[name] = {dim, tag};
    }

    /*-----------get physical properties of domain----------------------*/
    vector<string> keys;
    onelab::getNames(keys, "(Materials).+");
    double E;
    double nu;
    for (auto &key : keys)
    {
        // get corresponding value
        vector<double> value;
        onelab::getNumber(key, value);
        // expected key structure is "type/group_name/field"
        // => split the key string into components
        stringstream ss(key);
        vector<string> words;
        string word;
        while(getline(ss, word, '/')) // read string until '/'
            words.push_back(word);
        if(words.size()==3) 
        {
            if(words[2].compare("Young") == 0){
                E = value[0];
            }
            if(words[2].compare("Poisson") == 0){
                nu = value[0];
            }
        }
    }
    // cout << "E = " << E << "\n";
    // cout << "nu = " << nu << "\n";

    /*-----------------H matrix------------------*/
    Matrix<double, 3, 3> H_matrix;
    H_matrix(0,0) = 1.;
    H_matrix(1,0) = nu;
    H_matrix(2,0) = 0.;
    H_matrix(0,1) = nu;
    H_matrix(1,1) = 1.;
    H_matrix(2,1) = 0.;
    H_matrix(0,2) = 0.;
    H_matrix(1,2) = 0.;
    H_matrix(2,2) = (1-nu)/2;
    H_matrix = E/(1-nu*nu)*H_matrix;
    // cout << H_matrix << '\n'; 

    /*------------loop over elements of domain------------------------*/
    string groupname = "domain";

    int dim = groups[groupname].first;
    int tag = groups[groupname].second;
    if (tag == 0)
    {
        cerr << "Group '" << groupname << "' does not exist!\n";
        return;
    }

    // get entities of the domain
    vector<int> tags;
    model::getEntitiesForPhysicalGroup(dim, tag, tags);

    cout << "Entities in group named '" << groupname << "': ";
    for (int i = 0; i < tags.size(); ++i)
        cout << tags[i] << " ";
    cout << '\n';

    // loop over entities
    for (int k = 0; k < tags.size(); ++k)
    {
        cout << "Entity (" << dim << "," << tags[k] << "):\n";
        vector<int> elementTypes;
        vector<vector<size_t>> elementTags;
        vector<vector<size_t>> nodeTags;
        model::mesh::getElements(elementTypes, elementTags,
                                       nodeTags, dim, tags[k]);
        
        // loop over element types
        for (int i = 0; i < elementTypes.size(); ++i)
        {
            cout << "\telement type " << elementTypes[i] << ":\n";
            cout << "\t\telements: ";
            for (int j = 0; j < elementTags[i].size(); ++j)
                cout << elementTags[i][j] << " ";
            cout << '\n';
            
            // get Gauss points coordinates and weights
            vector<double> localCoords, weights;
            model::mesh::getIntegrationPoints(elementTypes[i], 
                                                    "CompositeGauss2", // "Gauss2"
                                                    localCoords, weights);
            cout << "\t" << weights.size() << " integration points\n";

            cout << "\tweights: ";
            for (int j = 0; j < weights.size(); ++j)
                cout << weights[j] << " ";
            cout << '\n';

            cout << "\tlocalCoords (gmsh): ";
            for (int j = 0; j < localCoords.size(); ++j)
                cout << localCoords[j] << " ";
            cout << '\n';

            // convert to Eigen format
            Map<Matrix<double, Dynamic, Dynamic>>
                localCoords_E(&localCoords[0], localCoords.size()/3, 3);
            cout << "\tlocalCoords (Eigen):\n";
            IOFormat fmt(4, 0, ", ", "\n", "\t\t[", "]");
            cout << localCoords_E.format(fmt) << '\n';

            // get determinants and Jacobians
            vector<double> jacobians, determinants, coords;
            model::mesh::getJacobians(elementTypes[i], 
                                            localCoords, jacobians, 
                                            determinants, coords);
            cout << "\tGot " << determinants.size() 
                    << " Jacobians for type " << elementTypes[i] << "\n";

            cout << "\tdeterminants (gmsh): ";
            for (int j = 0; j < determinants.size(); ++j)
                cout << determinants[j] << " ";
            cout << '\n';

            cout << "\tjacobians (gmsh): ";
            for (int j = 0; j < jacobians.size(); ++j)
                cout << jacobians[j] << " ";
            cout << '\n';

            // derivatives of the shape functions at the Gauss points of the reference element
            // we are here working with derivatives in the reference space, same for all elements
            int numComponents, numOrientations;
            vector<double> basisFunctions;
            model::mesh::getBasisFunctions(elementTypes[i], localCoords, 
                                                 "GradLagrange", numComponents, 
                                                 basisFunctions, 
                                                 numOrientations);
            cout << "\tbasis functions (gmsh): ";
            for (int j = 0; j < basisFunctions.size(); ++j)
                cout << basisFunctions[j] << " ";
            cout << '\n';

            /*---------computing the B matrix of the first element's (i.e elementTags[i][0]) first GP:------*/
            // ne serait-ça pas plutôt utile de directement convertir les Jacobiens en 2D ? travailler en full 2D
            // et faire pareil pour l'évaluation des gradients des basis functions ci-dessus ?
            vector<double> b_matrix;

            // convert first jacobian to Eigen format
            Map<Matrix<double, Dynamic, Dynamic> >
                jacob(&jacobians[0], 3, 3);
            // passage en 2D xD
            Matrix<double, 2, 2> jacob2D;
            jacob2D(0,0) = jacob(0,0);
            jacob2D(1,0) = jacob(1,0);
            jacob2D(0,1) = jacob(0,1);
            jacob2D(1,1) = jacob(1,1);
            cout << "\tfirst jacobian (Eigen):\n" 
                      << jacob2D.format(fmt) << '\n';
                    
            // compute the inverse with Eigen
            Matrix<double, 2, 2> jacobinv = jacob2D.inverse();
            cout << "\tinverse of the first jacobian (Eigen):\n" 
                      << jacobinv.format(fmt) << "\n";

            // compute the transpose of the inverse with Eigen
            Matrix<double, 2, 2> jacobinvtrans = jacobinv.transpose();
            cout << "\ttranspose of the inverse of the first jacobian (Eigen):\n" 
                      << jacobinvtrans.format(fmt) << "\n";
            Matrix<double, 3, Dynamic> B_matrix;
            B_matrix(0,0) = 0.;
            cout << "\ttranspose of the inverse of the first jacobian (Eigen):\n" 
                      << B_matrix.format(fmt) << "\n";
            for (int j = 0; j < weights.size(); ++j){ // looping over the number of gauss points ( = number of shape functions (A VERIF))
                
            }
        }
    }
}

void Tutorial_sparse_example_details(int argc, char **argv)
{
    typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
    typedef Triplet<double> T;
}

void buildProblem(vector<T>& coefficients, VectorXd& b, int n)
{
  b.setZero();
  ArrayXd boundary = ArrayXd::LinSpaced(n, 0,M_PI).sin().pow(2);
  for(int j=0; j<n; ++j)
  {
    for(int i=0; i<n; ++i)
    {
      int id = i+j*n;
      insertCoefficient(id, i-1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i+1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j-1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j+1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j,    4, coefficients, b, boundary);
    }
  }
}

void insertCoefficient(int id, int i, int j, double w, vector<T>& coeffs, VectorXd& b, const VectorXd& boundary)
{
    int n = int(boundary.size());
    int id1 = i+j*n;

    if(i==-1 || i==n) b(id) -= w * boundary(j); // constrained coefficient
    else  if(j==-1 || j==n) b(id) -= w * boundary(i); // constrained coefficient
    else  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}

// void saveAsBitmap(const VectorXd& x, int n, const char* filename)
// {
//   Array<unsigned char,Dynamic,Dynamic> bits = (x*255).cast<unsigned char>();
//   QImage img(bits.data(), n,n,QImage::Format_Indexed8);
//   img.setColorCount(256);
//   for(int i=0;i<256;i++) img.setColor(i,qRgb(i,i,i));
//   img.save(filename);
// }

void Tutorial_sparse_example(int argc, char **argv)
{
    typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
    typedef Triplet<double> T;

    void buildProblem(vector<T>& coefficients, VectorXd& b, int n);
    void saveAsBitmap(const VectorXd& x, int n, const char* filename);

    // if(argc!=2) {
    //   cerr << "Error: expected one and only one argument.\n";
    //   return -1;
    // }
  
    int n = 300;  // size of the image
    int m = n*n;  // number of unknows (=number of pixels)

    // Assembly:
    vector<T> coefficients;            // list of non-zeros coefficients
    VectorXd b(m);                   // the right hand side-vector resulting from the constraints
    buildProblem(coefficients, b, n);

    SpMat A(m,m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Solving:
    SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
    VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side

    // Export the result to a file:
    // saveAsBitmap(x, n, argv[1]);
}