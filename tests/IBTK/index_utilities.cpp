#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscao.h>

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

    if (app_initializer->getInputDatabase()->getBoolWithDefault("test_mapping", false))
    {
        hier::Index<NDIM> lower(-3), extent(4);
        lower(1) = 2;
        extent(1) = 3;
#if (NDIM == 3)
        lower(2) = -1;
        extent(2) = 2;
#endif
        const hier::Index<NDIM> upper = lower + extent - 1;
        const int volume = extent.getProduct();
        const int offset = 17;
        const int depth = 2;
        const int first = offset + depth * volume;
        int mismatches = 0;
        mismatches += IndexUtilities::mapIndexToInteger(lower, lower, extent, depth, offset) != first;
        mismatches += IndexUtilities::mapIndexToInteger(upper, lower, extent, depth, offset) != first + volume - 1;
        mismatches += IndexUtilities::mapIndexToInteger(lower - 1, lower, extent, depth, offset) != -1;
        mismatches += IndexUtilities::mapIndexToInteger(upper + 1, lower, extent, depth, offset) != -1;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            hier::Index<NDIM> below = lower, above = upper;
            --below(axis);
            ++above(axis);
            mismatches += IndexUtilities::mapIndexToInteger(below, lower, extent, depth, offset) != -1;
            mismatches += IndexUtilities::mapIndexToInteger(above, lower, extent, depth, offset) != -1;
            hier::IntVector<NDIM> periodic(0);
            periodic(axis) = extent(axis);
            hier::Index<NDIM> partner = below;
            partner(axis) = upper(axis);
            mismatches += IndexUtilities::mapIndexToInteger(below, lower, extent, depth, offset, periodic) !=
                          IndexUtilities::mapIndexToInteger(partner, lower, extent, depth, offset);
            partner = above;
            partner(axis) = lower(axis);
            mismatches += IndexUtilities::mapIndexToInteger(above, lower, extent, depth, offset, periodic) !=
                          IndexUtilities::mapIndexToInteger(partner, lower, extent, depth, offset);
            --below((axis + 1) % NDIM);
            mismatches += IndexUtilities::mapIndexToInteger(below, lower, extent, depth, offset, periodic) != -1;
        }
        hier::IntVector<NDIM> periodic(extent);
        mismatches +=
            IndexUtilities::mapIndexToInteger(lower - 1, lower, extent, depth, offset, periodic) != first + volume - 1;
        mismatches += IndexUtilities::mapIndexToInteger(upper + 1, lower, extent, depth, offset, periodic) != first;
        PetscInt application[2] = { first, first + volume - 1 };
        PetscInt petsc[2] = { 1, 0 };
        AO ordering = nullptr;
        int ierr = AOCreateMapping(PETSC_COMM_SELF, 2, application, petsc, &ordering);
        IBTK_CHKERRQ(ierr);
        PetscInt indices[3] = { -1, application[0], application[1] };
        ierr = AOApplicationToPetsc(ordering, 3, indices);
        IBTK_CHKERRQ(ierr);
        mismatches += indices[0] != -1 || indices[1] != 1 || indices[2] != 0;
        ierr = AODestroy(&ordering);
        IBTK_CHKERRQ(ierr);
        tbox::plog << "mapping mismatches = " << mismatches << '\n';
        return mismatches ? 1 : 0;
    }

    auto tuple = setup_hierarchy<NDIM>(app_initializer);
    auto patch_hierarchy = std::get<0>(tuple);

    auto print_index = [&](const IBTK::Point& p)
    {
        tbox::pout << "Point = " << p << '\n';
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::pout << "  Level = " << ln << '\n';
            tbox::Pointer<hier::PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (hier::PatchLevel<NDIM>::Iterator it(level); it; it++)
            {
                tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(it());
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
