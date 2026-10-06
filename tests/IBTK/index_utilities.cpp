#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IndexUtilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <HierarchyDataOpsManager.h>

#include <array>
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
        const hier::Index<NDIM> upper = lower + extent - 1;
        const int offset = 17;
        const int depth = 2;
        auto map_index = [&](const hier::Index<NDIM>& i, const hier::IntVector<NDIM>& periodic_shift)
        { return IndexUtilities::mapIndexToInteger(i, lower, extent, depth, offset, periodic_shift); };
        const hier::IntVector<NDIM> no_shift(0);
        tbox::pout << "lower = " << map_index(lower, no_shift) << '\n'
                   << "upper = " << map_index(upper, no_shift) << '\n';
        for (int axis = 0; axis < NDIM; ++axis)
        {
            hier::Index<NDIM> below = lower, above = upper;
            --below(axis);
            ++above(axis);
            hier::IntVector<NDIM> periodic_shift(0);
            periodic_shift(axis) = extent(axis);
            tbox::pout << "axis " << axis << ":\n"
                       << "  below = " << map_index(below, no_shift) << '\n'
                       << "  above = " << map_index(above, no_shift) << '\n'
                       << "  below, periodic = " << map_index(below, periodic_shift) << '\n'
                       << "  above, periodic = " << map_index(above, periodic_shift) << '\n';
            // An index that is also outside the array along a nonperiodic axis.
            --below((axis + 1) % NDIM);
            tbox::pout << "  below along two axes, periodic = " << map_index(below, periodic_shift) << '\n';
        }
        const hier::IntVector<NDIM> periodic_shift(extent);
        tbox::pout << "lower - 1, periodic = " << map_index(lower - 1, periodic_shift) << '\n'
                   << "upper + 1, periodic = " << map_index(upper + 1, periodic_shift) << '\n';
        return 0;
    }

    if (app_initializer->getInputDatabase()->getBoolWithDefault("test_clamp", false))
    {
        // clampToDomain() returns a point in [x_lower, x_upper), and that point
        // is in an interior cell. The last three domains are so far from the
        // origin, or so short, that the offset of clampToDomain() is smaller
        // than the spacing of the floating-point values near x_upper.
        const double domains[][2] = { { 0.0, 1.0 },
                                      { -4.0, 4.0 },
                                      { 1.0e9, 1.0e9 + 1.0 },
                                      { -1.0e9 - 1.0, -1.0e9 },
                                      { 1.0, std::nextafter(1.0, 2.0) } };
        const int n_cells = 16;
        const hier::Index<NDIM> ilower(0), iupper(n_cells - 1);
        auto in_domain = [](const double x, const double x_lower, const double x_upper)
        { return x_lower <= x && x < x_upper; };
        int domain_n = 0;
        for (const auto& domain : domains)
        {
            const double x_lower = domain[0], x_upper = domain[1];
            const double x_below = IndexUtilities::clampToDomain(x_lower - 1.0, x_lower, x_upper);
            const double x_on_lower = IndexUtilities::clampToDomain(x_lower, x_lower, x_upper);
            const double x_on_upper = IndexUtilities::clampToDomain(x_upper, x_lower, x_upper);
            const double x_above = IndexUtilities::clampToDomain(x_upper + 1.0, x_lower, x_upper);
            std::array<double, NDIM> X, X_lower, X_upper, dx;
            X.fill(x_on_upper);
            X_lower.fill(x_lower);
            X_upper.fill(x_upper);
            dx.fill((x_upper - x_lower) / n_cells);
            const hier::Index<NDIM> index =
                IndexUtilities::getCellIndex(X, X_lower.data(), X_upper.data(), dx.data(), ilower, iupper);
            tbox::pout << "domain " << domain_n++ << ":\n"
                       << "  below the domain is moved to the lower boundary = " << (x_below == x_lower) << '\n'
                       << "  the lower boundary is not moved = " << (x_on_lower == x_lower) << '\n'
                       << "  the upper boundary is moved into the domain = " << in_domain(x_on_upper, x_lower, x_upper)
                       << '\n'
                       << "  above the domain is moved into the domain = " << in_domain(x_above, x_lower, x_upper)
                       << '\n'
                       << "  cell of the upper boundary is interior = "
                       << (ilower(0) <= index(0) && index(0) <= iupper(0)) << '\n';
        }
        return 0;
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
