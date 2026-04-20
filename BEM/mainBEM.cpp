#include "functionsBEM.hpp"
#include <gmsh.h>
#include <iostream>

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <geo_file>\n";
        return 0;
    }

    int nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads);

    gmsh::initialize();
    omp_set_num_threads(nthreads);

    gmsh::open(argv[1]);
    gmsh::model::mesh::generate(2);

    std::map<int, double> electrostaticPressure;
    int nbViews = 0;
    std::map<int, std::pair<double,double>> nodalDisplacements;
    const bool postProcessing = true;
    const bool meshUntangler  = false;

    // solverBEM automatically selects single-domain or multi-domain path.
    solverBEM(electrostaticPressure, nbViews, nodalDisplacements,
              postProcessing, 1, meshUntangler);

    if(postProcessing)
        gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
