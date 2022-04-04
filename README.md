# Multiphysics 0471 

Multiphysics Integrated Computational Project

FOR WINDOWS: 

To build the FEM solver:
- go to the srcs/FEM folder (cd srcs/FEM)
- create a build folder (mkdir build)
- go to this build folder (cd build)
- cmake ..
- make

To execute the FEM solver:
- stay in the build folder and solver.exe ..\my_geo.geo
- for simple tension configuration: solver.exe ..\simple_tension.geo
- for second validation configuration: solver.exe ..\uniform_charge.geo