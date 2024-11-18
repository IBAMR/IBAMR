#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IndexUtilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <HierarchyDataOpsManager.h>

#include <cmath>
#include <fstream>

#include "../tests.h"

// Test some functions in IndexUtilities by computing cells and cell centers.

int
main(int argc, char** argv)
{
    using namespace SAMRAI;
    using namespace IBTK;

    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    tbox::Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);

    auto tuple = setup_hierarchy<NDIM>(app_initializer);
    auto patch_hierarchy = std::get<0>(tuple);

    auto print_index = [&](const IBTK::Point& p)
    {
        tbox::pout << "Point = " << p << '\n';
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::pout << "  Level = " << ln << '\n';
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (hier::PatchLevel<NDIM>::Iterator it(level); it; it++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(it());
                const auto index = IndexUtilities::getCellIndex(p, patch->getPatchGeometry(), patch->getBox());
                IBTK::Vector c0, c1;
                c0 = IndexUtilities::getCellCenter(*patch, index);
                c1 = IndexUtilities::getCellCenter<IBTK::Vector>(
                    patch_hierarchy->getGridGeometry(), level->getRatio(), index);
                TBOX_ASSERT(c0 == c1);
                tbox::pout << "    Box         = " << patch->getBox() << '\n'
                           << "    Index       = " << index << '\n'
                           << "    cell center = " << c0.transpose() << '\n'
                           << "    cell center = " << c1.transpose() << '\n'
                           << "    contains    = " << patch->getBox().contains(index) << '\n';
            }
        }
    };

    {
        IBTK::Point p0;
        p0.setZero();
        p0[0] = -0.35355339059327373086;
        p0[1] = -0.35355339059327373086;

        print_index(p0);
    }
    {
        IBTK::Point p0;
        p0.setZero();
        p0[0] = -5.5511151231257827e-17;
        p0[1] = -0.49999999999999994;

        print_index(p0);
    }
    {
        IBTK::Point p0;
        p0.setZero();
        p0[0] = -0.25;
        p0[1] = -0.25;

        print_index(p0);
    }
    {
        IBTK::Point p0;
        p0.setZero();
        p0[0] = -2.7755575615628914e-17;
        p0[1] = -0.32322330470336308;

        print_index(p0);
    }
    {
        IBTK::Point p0;
        p0.setZero();
        p0[0] = -4.163336342344337e-17;
        p0[1] = -0.41161165235168151;

        print_index(p0);
    }

    app_initializer->getVisItDataWriter()->writePlotData(patch_hierarchy, 0, 0.0);
}
