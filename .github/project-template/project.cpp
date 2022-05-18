#include <ibtk/IBTKInit.h>
#include <ibtk/SAMRAIDataCache.h>

#include <VariableDatabase.h>

#include <Eigen/Core>

#ifdef IBTK_HAVE_LIBMESH
#include <libmesh/point.h>
#endif

#include <iostream>
#include <memory>

// Test IBAMR project. This verifies that we can correctly compile and link against IBAMR

int main(int argc, char **argv)
{
    // Make sure we set up the right C++ version by explicitly using C++11 features:
    auto ibtk_init = std::make_shared<IBTK::IBTKInit>(argc, argv);
#ifdef IBTK_HAVE_LIBMESH
    const libMesh::Point point(1.0, 2.0, 3.0);
#endif

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    std::cout << "Number of contexts: " << var_db->getNumberOfRegisteredPatchDataIndices() << '\n';

    std::cout << "hello, world\n";

    std::cout << "eigen version: " << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.'
              << EIGEN_MINOR_VERSION << '\n';
}
