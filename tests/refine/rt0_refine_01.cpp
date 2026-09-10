// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SAMRAIVectorReal.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <algorithm>
#include <cmath>
#include <iomanip>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// Check specialized-linear refinement of constant and globally linear fields.
// Also verify that we can correctly refine a piecewise linear solution with the
// RT0 refinement class. Since the RT0 element is a vector-valued element that
// is, on Cartesian grids,
//
//     RT0_K = (P^1(x) * P^0(y), P^0(x) * P^1(y))
//
// we expect refining a vector field which is in that space to have zero error.

namespace
{
void
test_specialized_linear_refine()
{
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u", 2);
    const int u_idx = var_db->registerVariableAndContext(u_var, var_db->getContext("specialized"), IntVector<NDIM>(1));
    CartSideDoubleSpecializedLinearRefine refine_op;
    constexpr double SENTINEL = 1.e30;

    for (int ratio_number = 0; ratio_number < 3; ++ratio_number)
    {
        IntVector<NDIM> ratio(ratio_number == 0 ? 2 : 4);
        Box<NDIM> fine_box;
        for (int d = 0; d < NDIM; ++d)
        {
            if (ratio_number == 2) ratio(d) = d + 2;
            fine_box.lower(d) = (-3 + d) * ratio(d);
            fine_box.upper(d) = (2 + d) * ratio(d) - 1;
        }
        const Box<NDIM> fine_ghost_box = Box<NDIM>::grow(fine_box, IntVector<NDIM>(1));
        // Cover the requested fine ghosts before adding the coarse slope stencil.
        const Box<NDIM> coarse_box = Box<NDIM>::coarsen(fine_ghost_box, ratio);
        Patch<NDIM> coarse(coarse_box, var_db->getPatchDescriptor());
        Patch<NDIM> fine(fine_box, var_db->getPatchDescriptor());
        coarse.allocatePatchData(u_idx, 0.0);
        fine.allocatePatchData(u_idx, 0.0);
        Pointer<SideData<NDIM, double>> coarse_data = coarse.getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> fine_data = fine.getPatchData(u_idx);

        for (const bool linear : { false, true })
        {
            // Side coordinates are integral in the normal direction and half-integral otherwise.
            auto exact = [&](const hier::Index<NDIM>& i, const int axis, const int depth, const IntVector<NDIM>& scale)
            {
                double value = 10.0 * (axis + 1) + depth;
                if (linear)
                    for (int d = 0; d < NDIM; ++d)
                        value += (axis + 1) * (d + 2) * (depth + 1) * (i(d) + (d == axis ? 0.0 : 0.5)) / scale(d);
                return value;
            };
            for (int axis = 0; axis < NDIM; ++axis)
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(coarse_data->getGhostBox(), axis)); it; it++)
                    for (int depth = 0; depth < coarse_data->getDepth(); ++depth)
                        (*coarse_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower), depth) =
                            exact(it(), axis, depth, IntVector<NDIM>(1));

            double max_error = 0.0, upper_plane_max_error = 0.0;
            int outside_errors = 0;
            // Test the interior, the full ghost box, and each lower/upper ghost slab.
            for (int region = 0; region < 2 + 2 * NDIM; ++region)
            {
                Box<NDIM> destination = region == 0 ? fine_box : fine_ghost_box;
                if (region >= 2)
                {
                    const int normal = (region - 2) / 2;
                    const int boundary = region % 2 == 0 ? fine_ghost_box.lower(normal) : fine_ghost_box.upper(normal);
                    destination.lower(normal) = destination.upper(normal) = boundary;
                }
                fine_data->fillAll(SENTINEL);
                refine_op.refine(fine, coarse, u_idx, u_idx, destination, ratio);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(destination, axis);
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(fine_data->getGhostBox(), axis)); it;
                         it++)
                        for (int depth = 0; depth < fine_data->getDepth(); ++depth)
                        {
                            const double value =
                                (*fine_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower), depth);
                            if (side_box.contains(it()))
                            {
                                if (!std::isfinite(value)) TBOX_ERROR("Refinement produced a nonfinite value.\n");
                                const double error = std::abs(value - exact(it(), axis, depth, ratio));
                                max_error = std::max(max_error, error);
                                if (it()(axis) == side_box.upper(axis))
                                    upper_plane_max_error = std::max(upper_plane_max_error, error);
                            }
                            else
                            {
                                outside_errors += value != SENTINEL;
                            }
                        }
                }
            }
            plog << (linear ? "linear" : "constant") << " ratio " << ratio << ": max error = " << max_error
                 << ", upper plane max error = " << upper_plane_max_error << ", outside writes = " << outside_errors
                 << '\n';
        }
    }
}

// The fine-cell divergence equals the parent-cell divergence for RT0 refinement.
int
test_rt0_schedule(Pointer<PatchHierarchy<NDIM>> hierarchy, int u_idx, int refined_idx)
{
    Pointer<PatchLevel<NDIM>> coarse = hierarchy->getPatchLevel(0);
    Pointer<PatchLevel<NDIM>> fine = hierarchy->getPatchLevel(1);
    if (coarse->getNumberOfPatches() != 1 || fine->getNumberOfPatches() != 1)
    {
        TBOX_ERROR("RT0 invariant fixture requires one patch per level.\n");
    }
    Pointer<Patch<NDIM>> coarse_patch = coarse->getPatch(0);
    Pointer<Patch<NDIM>> fine_patch = fine->getPatch(0);
    Pointer<SideData<NDIM, double>> coarse_data = coarse_patch->getPatchData(u_idx);
    Pointer<SideData<NDIM, double>> refined_data = fine_patch->getPatchData(refined_idx);
    Pointer<CartesianPatchGeometry<NDIM>> coarse_geom = coarse_patch->getPatchGeometry();
    Pointer<CartesianPatchGeometry<NDIM>> fine_geom = fine_patch->getPatchGeometry();
    const IntVector<NDIM> ratio = fine->getRatioToCoarserLevel();
    Pointer<RefineAlgorithm<NDIM>> algorithm = new RefineAlgorithm<NDIM>();
    algorithm->registerRefine(refined_idx, u_idx, refined_idx, new CartSideDoubleRT0Refine());
    Pointer<RefineSchedule<NDIM>> schedule = algorithm->createSchedule(fine, Pointer<PatchLevel<NDIM>>(), 0, hierarchy);
    double affine_error = 0.0, divergence_error = 0.0, divergence_norm = 0.0;
    for (int profile = 0; profile < 3; ++profile)
    {
        coarse_data->fillAll(0.0);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(coarse_patch->getBox(), axis)); b; b++)
            {
                double value = axis + 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double x = coarse_geom->getXLower()[d] +
                                     coarse_geom->getDx()[d] *
                                         (b()(d) - coarse_patch->getBox().lower()(d) + (d == axis ? 0.0 : 0.5));
                    if (profile == 0 && d == axis)
                    {
                        value += 0.3 * (axis + 1) * x;
                    }
                    else if (profile == 1)
                    {
                        value += 0.2 * (d + 1) * std::sin(2.0 * std::acos(-1.0) * x);
                    }
                    else if (profile == 2)
                    {
                        value += d == axis ? 0.3 * std::abs(2.0 * x - 1.0) : 0.2 * std::floor(4.0 * x);
                    }
                }
                (*coarse_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) = value;
            }
        }
        refined_data->fillAll(0.0);
        schedule->fillData(0.0);
        if (profile == 0)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_patch->getBox(), axis)); b; b++)
                {
                    const double x = fine_geom->getXLower()[axis] +
                                     fine_geom->getDx()[axis] * (b()(axis) - fine_patch->getBox().lower()(axis));
                    const double exact = axis + 1.0 + 0.3 * (axis + 1) * x;
                    const double error =
                        std::abs((*refined_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) - exact);
                    if (!std::isfinite(error))
                    {
                        TBOX_ERROR("Nonfinite RT0 affine error.\n");
                    }
                    affine_error = std::max(affine_error, error);
                }
            }
        }
        for (Box<NDIM>::Iterator b(fine_patch->getBox()); b; b++)
        {
            const hier::Index<NDIM> coarse_index = IndexUtilities::coarsen(b(), ratio);
            double coarse_div = 0.0, fine_div = 0.0;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                coarse_div += ((*coarse_data)(SideIndex<NDIM>(coarse_index, axis, SideIndex<NDIM>::Upper)) -
                               (*coarse_data)(SideIndex<NDIM>(coarse_index, axis, SideIndex<NDIM>::Lower))) /
                              coarse_geom->getDx()[axis];
                fine_div += ((*refined_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Upper)) -
                             (*refined_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower))) /
                            fine_geom->getDx()[axis];
            }
            const double error = std::abs(fine_div - coarse_div);
            if (!std::isfinite(error))
            {
                TBOX_ERROR("Nonfinite RT0 divergence error.\n");
            }
            divergence_error = std::max(divergence_error, error);
            divergence_norm = std::max(divergence_norm, std::abs(coarse_div));
        }
    }
    plog << std::setprecision(12) << "affine schedule error = " << affine_error << '\n'
         << "divergence error = " << divergence_error << '\n'
         << "coarse divergence norm = " << divergence_norm << '\n';
    return affine_error < 1.0e-12 && divergence_error < 1.0e-12 && divergence_norm > 0.0 ? 0 : 1;
}

} // namespace

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // this test only works in serial
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "rt0.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        if (input_db->getBoolWithDefault("test_specialized_linear_refine", false))
        {
            test_specialized_linear_refine();
            return 0;
        }

        // Create major algorithm and data objects that comprise the
        // application.
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

        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");
        Pointer<SideVariable<NDIM, double>> u_sc_var = new SideVariable<NDIM, double>("u_sc");
        const IntVector<NDIM> ghosts(input_db->getBoolWithDefault("test_rt0_schedule", false) ? 1 : 0);
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, ghosts);
        Pointer<SideVariable<NDIM, double>> exact_sc_var = new SideVariable<NDIM, double>("exact_sc");
        const int exact_sc_idx = var_db->registerVariableAndContext(exact_sc_var, ctx, ghosts);
        // u_cc_var is only for plotting (and testing): uncomment if output is desired
// #define DO_PLOT
#ifdef DO_PLOT
        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx);
        Pointer<CellVariable<NDIM, double>> exact_cc_var = new CellVariable<NDIM, double>("exact_cc", NDIM);
        const int exact_cc_idx = var_db->registerVariableAndContext(exact_cc_var, ctx);
#endif

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(exact_sc_idx, 0.0);
#ifdef DO_PLOT
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(exact_cc_idx, 0.0);
#endif
        }

        Pointer<VisItDataWriter<NDIM>> visit_writer = app_initializer->getVisItDataWriter();

        if (input_db->getBoolWithDefault("test_rt0_schedule", false))
        {
            return test_rt0_schedule(patch_hierarchy, u_sc_idx, exact_sc_idx);
        }

        // The rest is just book-keeping, this is the actual test:
        auto do_test = [&](const std::string& db_u_fcn_name, const int coarse_level_n)
        {
            muParserCartGridFunction u_fcn(
                db_u_fcn_name, app_initializer->getComponentDatabase(db_u_fcn_name), grid_geometry);
            u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
            u_fcn.setDataOnPatchHierarchy(exact_sc_idx, exact_sc_var, patch_hierarchy, 0.0);

            solv::SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, 1);
            u_vec.addComponent(u_sc_var, u_sc_idx);

            const int fine_level_n = coarse_level_n + 1;
            Pointer<PatchLevel<NDIM>> level_0 = patch_hierarchy->getPatchLevel(coarse_level_n);
            Pointer<PatchLevel<NDIM>> level_1 = patch_hierarchy->getPatchLevel(fine_level_n);

            // there should only be one patch on each patch level
            Pointer<SideData<NDIM, double>> u_sc_0_data = level_0->getPatch(0)->getPatchData(u_sc_idx);
            Pointer<SideData<NDIM, double>> u_sc_1_data = level_1->getPatch(0)->getPatchData(u_sc_idx);
            const Box<NDIM> patch_box_0 = level_0->getPatch(0)->getBox();
            const Box<NDIM> patch_box_1 = level_1->getPatch(0)->getBox();

            Pointer<SideData<NDIM, double>> exact_sc_0_data = level_0->getPatch(0)->getPatchData(exact_sc_idx);
            Pointer<SideData<NDIM, double>> exact_sc_1_data = level_1->getPatch(0)->getPatchData(exact_sc_idx);

            const IntVector<NDIM> ratio = level_1->getRatioToCoarserLevel();
            IBTK::CartSideDoubleRT0Refine refine_op;
            constexpr double SENTINEL = 1.e30;
            u_sc_1_data->fillAll(SENTINEL);
            refine_op.refine(*level_1->getPatch(0), *level_0->getPatch(0), u_sc_idx, u_sc_idx, patch_box_1, ratio);

            for (int axis = 0; axis < NDIM; ++axis)
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box_1, axis)); it; it++)
                {
                    const SideIndex<NDIM> side(it(), axis, SideIndex<NDIM>::Lower);
                    const double value = (*u_sc_1_data)(side);
                    if (!std::isfinite(value) || value == SENTINEL)
                        TBOX_ERROR("RT0 refinement left a destination side unfilled.\n");
                }

            solv::SAMRAIVectorReal<NDIM, double> exact_vec("e", patch_hierarchy, coarse_level_n, fine_level_n);
            exact_vec.addComponent(exact_sc_var, exact_sc_idx);
            exact_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false),
                               Pointer<SAMRAIVectorReal<NDIM, double>>(&exact_vec, false));

            pout << "test results for " << db_u_fcn_name << '\n';
            pout << "max norm of u_sc: " << u_vec.maxNorm() << '\n';
            pout << "max norm of exact - refined: " << std::abs(exact_vec.maxNorm()) << '\n';
        };

        do_test("constant_function", 0);
        do_test("linear_x", 0);
        do_test("linear_y", 0);
        if (NDIM == 3)
        {
            do_test("linear_z", 0);
        }
        do_test("both_linear", 0);

        for (int coarse_level_n = 0; coarse_level_n < gridding_algorithm->getMaxLevels() - 1; ++coarse_level_n)
        {
            do_test("non_rt", coarse_level_n);
        }

#ifdef DO_PLOT
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        visit_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, nullptr, 0.0, false);
        visit_writer->registerPlotQuantity(exact_cc_var->getName(), "VECTOR", exact_cc_idx);
        hier_math_ops.interp(exact_cc_idx, exact_cc_var, exact_sc_idx, exact_sc_var, nullptr, 0.0, false);
#endif

#ifdef DO_PLOT
        visit_writer->writePlotData(patch_hierarchy, 0, 0.0);
#endif
    }
}
