# ⚡ 2D Multi-Domain Boundary Element Method (BEM) Solver for Electrostatics

An open-source C++ framework designed to simulate 2D electrostatic fields ($\phi$, $\mathbf{E}$) and electric displacement discontinuities across multi-material interfaces using the **Boundary Element Method (BEM)**.

The solver couples with **Gmsh** (via `ONELAB` and `.geo` geometry scripts) and **Eigen**, utilizing **OpenMP** multi-threading for high-performance matrix assembly and post-processing field computations.

## 🌟 Key Features

* **Multi-Domain Electrostatics**: Computes potential $\phi$ and normal electric flux derivatives across heterogeneous dielectric media, interfaces, inclusions, and holes.

* **Interfacial Jump Conditions**: Accurately handles continuity of electric potential ($\phi_A = \phi_B$) and displacement jumps ($\varepsilon_A \frac{\partial \phi_A}{\partial n} - \varepsilon_B \frac{\partial \phi_B}{\partial n} = 0$) across continuous/dispersed phase boundaries.

* **Analytical BEM Kernels**: Uses fast and accurate analytical integrations for 2D Green's function potential ($G$) and normal derivative ($H$) kernels, as well as spatial gradients ($\nabla G$, $\nabla H$).

* **OpenMP Parallelization**: Multi-threaded assembly of linear systems and evaluation of internal post-processing grids.

* **Gmsh & ONELAB Integration**: Geometry and boundary conditions are driven directly through Gmsh `.geo` files, with dynamic physical group identification.

## 🏗️ Project Structure

```
.
├── boundaryElements.cpp   # Mesh parsing, element creation, and interface grouping
├── computations.cpp       # Analytical BEM kernel integrations (G, H, gradG, gradH)
├── functionsBEM.hpp       # Core data structures (elementStruct, domainData, etc.)
├── postProcessing.cpp     # Internal field evaluation (potential & electric field E)
├── main.cpp               # Solver workflow and ONELAB interface execution
└── README.md              # Documentation
```

## 🏷️ Gmsh Naming Conventions

The solver automatically parses physical groups from `.geo` files based on strict naming prefixes:

| Physical Group | Naming Convention | Description |
| :--- | :--- | :--- |
| **Domain Surface** | `BEM_domain_<Name>` | Defines 2D dielectric domains |
| **Outer Boundary** | `<BoundaryName>` | External boundaries (Dirichlet/Neumann configured in ONELAB) |
| **Inclusion Interface** | `BEM_interface_<DomainA>_<DomainB>_inc` | Closed boundary of an inclusion/droplet inside a domain |
| **Juxtaposed Interface**| `BEM_interface_<DomainA>_<DomainB>_jux` | Interface separating two adjoining domains |

> ⚠️ **Note**: Correct physical group naming in your `.geo` file is required for automatic interface normal orientation and jump condition setup.

## 🚀 Getting Started

### Prerequisites

* **C++17** compatible compiler (`g++`, `clang++`, or Xcode/Clang)
* [Gmsh SDK](https://gmsh.info/) (C++ API headers and linked library)
* [Eigen 3](https://eigen.tuxfamily.org/) header-only linear algebra library
* [OpenMP](https://www.openmp.org/) library for parallel processing

### Build & Compilation

If compiling directly via command line:

```bash
g++ -O3 -std=c++17 main.cpp boundaryElements.cpp computations.cpp postProcessing.cpp \
    -I/path/to/eigen3 \
    -I/path/to/gmsh/include \
    -L/path/to/gmsh/lib -lgmsh \
    -fopenmp \
    -o bem_solver
```

*(If using **Xcode**, ensure that the Eigen and Gmsh header search paths are included in Build Settings, and `libgmsh` is linked to your target).*

### Execution

Launch the solver with an input Gmsh geometry file:

```bash
# Set thread count for OpenMP
export OMP_NUM_THREADS=8

# Run the solver
./bem_solver path/to/geometry.geo
```

## 📜 References

* **Boundary Element Method in Electrostatics**: Brebbia, C. A., & Dominguez, J. *Boundary Elements: An Introductory Course*.
* **Microfluidic BEM Electrostatics**: Schütz et al., *Lab on a Chip*, 2019.
