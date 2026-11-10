// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>

#include <petscsys.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <BoundaryBox.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <StandardTagAndInitialize.h>

#include <functional>
#include <memory>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// Dirichlet data equal to a function at the physical boundary, evaluated at the cell centers along the boundary. For
// side-centered data, the boundary operator shifts the patch geometry.
class FunctionDirichletBc : public RobinBcCoefStrategy<NDIM>
{
public:
    explicit FunctionDirichletBc(std::function<double(const VectorNd&)> fcn) : d_fcn(std::move(fcn))
    {
    }

    void setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                    Pointer<ArrayData<NDIM, double>>& bcoef_data,
                    Pointer<ArrayData<NDIM, double>>& gcoef_data,
                    const Pointer<SAMRAI::hier::Variable<NDIM>>& /*variable*/,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double /*fill_time*/) const override
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const x_upper = pgeom->getXUpper();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& patch_lower = patch.getBox().lower();
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int normal_axis = location_index / 2;
        const Box<NDIM>& box = (acoef_data ? acoef_data->getBox() : gcoef_data->getBox());
        for (Box<NDIM>::Iterator b(box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            VectorNd x;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == normal_axis)
                    x[d] = location_index % 2 == 0 ? x_lower[d] : x_upper[d];
                else
                    x[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
            }
            if (acoef_data) (*acoef_data)(i, 0) = 1.0;
            if (bcoef_data) (*bcoef_data)(i, 0) = 0.0;
            if (gcoef_data) (*gcoef_data)(i, 0) = d_fcn(x);
        }
    }

    IntVector<NDIM> numberOfExtensionsFillable() const override
    {
        return IntVector<NDIM>(128);
    }

private:
    std::function<double(const VectorNd&)> d_fcn;
};
} // namespace

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress the SAMRAI warning for patches smaller than the ghost cell width, because the warning contains a path.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    // Record the message of an expected error in the output.
    Logger::getInstance()->setAbortAppender(Pointer<Logger::Appender>(new TestAppender()));

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        // Create cell-centered and node-centered quantities, and initialize them with a function read from the input
        // file
        Pointer<SAMRAI::hier::Variable<NDIM>> Q_var;
        if (input_db->getString("CENTERING").compare("CELL") == 0)
            Q_var = new CellVariable<NDIM, double>("Q");
        else if (input_db->getString("CENTERING").compare("SIDE") == 0)
            Q_var = new SideVariable<NDIM, double>("Q");

        // With EXTRAP_TYPE, the test fills the ghost cells with linear refinement and the extrapolation at the physical
        // boundary instead of the coarse-fine interpolation, for a linear function that both reproduce.
        const std::string extrap_type = input_db->getStringWithDefault("EXTRAP_TYPE", "NONE");
        const bool test_extrap = extrap_type != "NONE";
        // With UNRELATED_ROBIN_DATA, fill a linear function with Robin data unrelated to it, like the pressure: only
        // the linear extension of the coarse data reproduces it.
        const bool unrelated_robin_data = input_db->getBoolWithDefault("UNRELATED_ROBIN_DATA", false);
        const bool extend_linearly = input_db->getBoolWithDefault("EXTEND_COARSE_DATA_LINEARLY", false);
        const bool Q_is_linear = test_extrap || unrelated_robin_data;
        auto Q_fcn = [Q_is_linear](const VectorNd& x) -> double
        {
            double val = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                val += x[d] + (Q_is_linear ? 0.0 : x[d] * x[d]);
                for (int e = d + 1; e < NDIM; ++e)
                {
                    val += Q_is_linear ? 0.0 : x[d] * x[e];
                }
            }
            return val;
        };

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        // The ghost cell width is an input. The extrapolation and split inputs specify 2.
        const int ghost_width = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        const int Q_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));
        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // Allocate patch data
        int coarsest_ln = 0, finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_idx);
        }

        // Fill in initial data
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(Q_idx);
                Pointer<SideData<NDIM, double>> Q_sc_data = patch->getPatchData(Q_idx);
                if (Q_cc_data)
                {
                    for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        (*Q_cc_data)(idx) = Q_fcn(x);
                    }
                }
                else if (Q_sc_data)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                                x[d] = xlow[d] +
                                       dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                            (*Q_sc_data)(idx) = Q_fcn(x);
                        }
                    }
                }
            }
        }

        // Test the coarse-fine interface operator alone, which fills one layer of ghost cells by quadratic
        // interpolation, also where the interface meets the physical boundary.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        std::vector<std::unique_ptr<FunctionDirichletBc>> bcs;
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs;
        if (!test_extrap && input_db->getStringWithDefault("BC_TYPE", "EXTRAPOLATION") == "ROBIN")
        {
            const int n_components = Q_var.isNull() || !Pointer<SideVariable<NDIM, double>>(Q_var) ? 1 : NDIM;
            for (int c = 0; c < n_components; ++c)
            {
                bcs.push_back(std::make_unique<FunctionDirichletBc>(
                    unrelated_robin_data ?
                        std::function<double(const VectorNd&)>([](const VectorNd&) { return 1000.0; }) :
                        std::function<double(const VectorNd&)>(Q_fcn)));
                bc_coefs.push_back(bcs.back().get());
            }
        }
        ghost_cell_comps[0] = test_extrap ?
                                  InterpolationTransactionComponent(Q_idx,
                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                    false,
                                                                    "CONSERVATIVE_COARSEN",
                                                                    extrap_type,
                                                                    false,
                                                                    nullptr /*bdry_conds*/) :
                                  InterpolationTransactionComponent(
                                      Q_idx, "NONE", true, "NONE", "LINEAR", false, bc_coefs, nullptr /*fill_pattern*/);
        ghost_cell_comps[0].d_extend_coarse_data_linearly = extend_linearly;
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);

        if (!test_extrap)
        {
            int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(Q_idx);
                Pointer<SideData<NDIM, double>> Q_sc_data = patch->getPatchData(Q_idx);
                // Only check the first layer of ghost cells on the left side, the layer that the interpolation fills.
                Box<NDIM> ghost_box = patch->getBox();
                ghost_box.upper(0) = ghost_box.lower(0);
                ghost_box.shift(0, -1);

                if (Q_cc_data)
                {
                    for (CellIterator<NDIM> ci(ghost_box); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        if (!IBTK::abs_equal_eps((*Q_cc_data)(idx)-Q_fcn(x), 0.0))
                        {
                            pout << "Incorrect ghost value!\n";
                            pout << "On ghost cell " << idx << "\n";
                            pout << "Computed value " << (*Q_cc_data)(idx) << "\n";
                            pout << "Exact value:   " << Q_fcn(x) << "\n";
                        }
                    }
                }
                else if (Q_sc_data)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(ghost_box, axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                                x[d] = xlow[d] +
                                       dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                            if (!IBTK::abs_equal_eps((*Q_sc_data)(idx)-Q_fcn(x), 0.0))
                            {
                                pout << "Incorrect ghost value!\n";
                                pout << "On ghost cell " << idx << " and axis " << axis << "\n";
                                pout << "Computed value " << (*Q_sc_data)(idx) << "\n";
                                pout << "Exact value:   " << Q_fcn(x) << "\n";
                            }
                        }
                    }
                }
            }
        }

        // Check that the ghost values one cell beyond each end of a patch along the left interface equal the exact
        // quadratic field. Such cells lie inside the domain only if several patches touch the interface.
        if (!test_extrap)
        {
            int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            const Box<NDIM> domain = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], level->getRatio());
            int num_checked = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(Q_idx);
                Pointer<SideData<NDIM, double>> Q_sc_data = patch->getPatchData(Q_idx);
                for (int t = 1; t < NDIM; ++t)
                {
                    for (const int end : { -1, +1 })
                    {
                        // The ghost cells beyond the end of the patch, away from its other tangential ends.
                        Box<NDIM> extension = patch->getBox();
                        extension.upper(0) = extension.lower(0) - 1;
                        extension.lower(0) = extension.upper(0);
                        for (int d = 1; d < NDIM; ++d)
                        {
                            if (d == t)
                            {
                                extension.lower(d) = extension.upper(d) =
                                    end < 0 ? patch->getBox().lower(d) - 1 : patch->getBox().upper(d) + 1;
                            }
                            else
                            {
                                extension.lower(d) += 2;
                                extension.upper(d) -= 2;
                            }
                        }
                        extension = extension * domain;
                        if (extension.empty()) continue;
                        if (Q_cc_data)
                        {
                            for (CellIterator<NDIM> ci(extension); ci; ci++)
                            {
                                const CellIndex<NDIM>& idx = ci();
                                VectorNd x;
                                for (int d = 0; d < NDIM; ++d)
                                    x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                                ++num_checked;
                                if (!IBTK::abs_equal_eps((*Q_cc_data)(idx)-Q_fcn(x), 0.0))
                                {
                                    pout << "Incorrect ghost value!\n";
                                    pout << "On extension cell " << idx << "\n";
                                    pout << "Computed value " << (*Q_cc_data)(idx) << "\n";
                                    pout << "Exact value:   " << Q_fcn(x) << "\n";
                                }
                            }
                        }
                        else if (Q_sc_data)
                        {
                            // The faces that are not normal to the direction of the extension.
                            for (int axis = 0; axis < NDIM; ++axis)
                            {
                                if (axis == t) continue;
                                for (SideIterator<NDIM> si(extension, axis); si; si++)
                                {
                                    const SideIndex<NDIM>& idx = si();
                                    VectorNd x;
                                    for (int d = 0; d < NDIM; ++d)
                                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) +
                                                                  (d == axis ? 0.0 : 0.5));
                                    ++num_checked;
                                    if (!IBTK::abs_equal_eps((*Q_sc_data)(idx)-Q_fcn(x), 0.0))
                                    {
                                        pout << "Incorrect ghost value!\n";
                                        pout << "On extension face " << idx << " and axis " << axis << "\n";
                                        pout << "Computed value " << (*Q_sc_data)(idx) << "\n";
                                        pout << "Exact value:   " << Q_fcn(x) << "\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (num_checked > 0) pout << "Checked " << num_checked << " values beyond the ends of patches.\n";
        }

        // Check the ghost cells between the lower boundary and a patch within the ghost cell width of it. The cells are
        // filled from a temporary coarse patch whose ghost cells are set by extrapolation.
        if (test_extrap)
        {
            const int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            const Box<NDIM> domain = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], level->getRatio());
            int num_checked = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(Q_idx);
                TBOX_ASSERT(Q_cc_data);
                // The ghost cells below the patch, away from its corners and from the other boundaries.
                Box<NDIM> below = patch->getBox();
                below.upper(1) = below.lower(1) - 1;
                below.lower(1) -= ghost_width;
                for (int d = 0; d < NDIM; ++d)
                {
                    if (d != 1)
                    {
                        below.lower(d) += 2;
                        below.upper(d) -= 2;
                    }
                }
                below = below * domain;
                if (below.empty()) continue;
                for (CellIterator<NDIM> ci(below); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    ++num_checked;
                    if (!IBTK::abs_equal_eps((*Q_cc_data)(idx)-Q_fcn(x), 0.0))
                    {
                        pout << "Incorrect ghost value!\n";
                        pout << "On ghost cell " << idx << "\n";
                        pout << "Computed value " << (*Q_cc_data)(idx) << "\n";
                        pout << "Exact value:   " << Q_fcn(x) << "\n";
                    }
                }
            }
            pout << "Checked " << num_checked << " values below patches.\n";
        }

        // Deallocate patch data
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
