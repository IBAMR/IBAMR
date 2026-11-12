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
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAIVectorReal.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleLinearRefine.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// Check specialized-linear and linear refinement of constant and globally linear fields.
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
test_exact_refine(const RefineOperator<NDIM>& refine_op)
{
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u", 2);
    const int u_idx = var_db->registerVariableAndContext(u_var, var_db->getContext("exact"), IntVector<NDIM>(1));
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
            // Report only errors above round-off.
            constexpr double ROUNDOFF = 1.0e-12;
            plog << (linear ? "linear" : "constant") << " ratio " << ratio
                 << ": max error = " << (max_error > ROUNDOFF ? max_error : 0.0)
                 << ", upper plane max error = " << (upper_plane_max_error > ROUNDOFF ? upper_plane_max_error : 0.0)
                 << ", outside writes = " << outside_errors << '\n';
        }
    }
}

// Refining a smooth field converges at second order.
void
test_smooth_refine(const RefineOperator<NDIM>& refine_op)
{
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u_smooth");
    const int u_idx = var_db->registerVariableAndContext(u_var, var_db->getContext("smooth"), IntVector<NDIM>(1));

    auto exact = [](const hier::Index<NDIM>& i, const int axis, const double dx)
    {
        double x[NDIM];
        for (int d = 0; d < NDIM; ++d) x[d] = (i(d) + (d == axis ? 0.0 : 0.5)) * dx;
        double value = std::sin(2.0 * M_PI * x[axis] + 0.3) * std::cos(2.0 * M_PI * x[(axis + 1) % NDIM] + 0.5);
        if (NDIM == 3) value *= std::cos(2.0 * M_PI * x[(axis + 2) % NDIM] + 0.1);
        return value;
    };

    for (const int r : { 2, 4 })
    {
        const IntVector<NDIM> ratio(r);
        double errors[2];
        for (int k = 0; k < 2; ++k)
        {
            const int N = 16 << k;
            Box<NDIM> fine_box;
            for (int d = 0; d < NDIM; ++d)
            {
                fine_box.lower(d) = (N / 4) * r;
                fine_box.upper(d) = (3 * N / 4) * r - 1;
            }
            const Box<NDIM> fine_ghost_box = Box<NDIM>::grow(fine_box, IntVector<NDIM>(1));
            Patch<NDIM> coarse(Box<NDIM>::coarsen(fine_ghost_box, ratio), var_db->getPatchDescriptor());
            Patch<NDIM> fine(fine_box, var_db->getPatchDescriptor());
            coarse.allocatePatchData(u_idx, 0.0);
            fine.allocatePatchData(u_idx, 0.0);
            Pointer<SideData<NDIM, double>> coarse_data = coarse.getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> fine_data = fine.getPatchData(u_idx);
            for (int axis = 0; axis < NDIM; ++axis)
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(coarse_data->getGhostBox(), axis)); it; it++)
                    (*coarse_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = exact(it(), axis, 1.0 / N);

            refine_op.refine(fine, coarse, u_idx, u_idx, fine_ghost_box, ratio);
            errors[k] = 0.0;
            for (int axis = 0; axis < NDIM; ++axis)
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(fine_ghost_box, axis)); it; it++)
                {
                    const double value = (*fine_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower));
                    errors[k] = std::max(errors[k], std::abs(value - exact(it(), axis, 1.0 / (N * r))));
                }
        }
        plog << "smooth ratio " << ratio << ": max errors = " << errors[0] << ", " << errors[1]
             << ", observed order = " << std::log2(errors[0] / errors[1]) << '\n';
    }
}

// Refinement is linear in the data and keeps values within the range of the coarse data.
void
test_random_refine(const RefineOperator<NDIM>& refine_op)
{
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    constexpr int N_FIELDS = 3;
    std::vector<int> u_idxs;
    for (int k = 0; k < N_FIELDS; ++k)
    {
        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u_random_" + std::to_string(k));
        u_idxs.push_back(var_db->registerVariableAndContext(u_var, var_db->getContext("random"), IntVector<NDIM>(1)));
    }
    constexpr double A = 0.7, B = -1.3;
    constexpr double TOLERANCE = 1.0e-13;
    std::mt19937 generator(0);

    for (const int r : { 2, 4 })
    {
        const IntVector<NDIM> ratio(r);
        Box<NDIM> fine_box;
        for (int d = 0; d < NDIM; ++d)
        {
            fine_box.lower(d) = -2 * r;
            fine_box.upper(d) = 3 * r - 1;
        }
        const Box<NDIM> fine_ghost_box = Box<NDIM>::grow(fine_box, IntVector<NDIM>(1));
        Patch<NDIM> coarse(Box<NDIM>::coarsen(fine_ghost_box, ratio), var_db->getPatchDescriptor());
        Patch<NDIM> fine(fine_box, var_db->getPatchDescriptor());
        std::vector<Pointer<SideData<NDIM, double>>> coarse_data, fine_data;
        for (const int u_idx : u_idxs)
        {
            coarse.allocatePatchData(u_idx, 0.0);
            fine.allocatePatchData(u_idx, 0.0);
            coarse_data.push_back(coarse.getPatchData(u_idx));
            fine_data.push_back(fine.getPatchData(u_idx));
        }

        // Fields 0 and 1 are random; field 2 is their linear combination.
        double lower = 1.0e30, upper = -1.0e30;
        for (int axis = 0; axis < NDIM; ++axis)
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(coarse_data[0]->getGhostBox(), axis)); it; it++)
            {
                const SideIndex<NDIM> i(it(), axis, SideIndex<NDIM>::Lower);
                (*coarse_data[0])(i) = 2.0 * static_cast<double>(generator() % 1000) / 999.0 - 1.0;
                (*coarse_data[1])(i) = 2.0 * static_cast<double>(generator() % 1000) / 999.0 - 1.0;
                (*coarse_data[2])(i) = A * (*coarse_data[0])(i) + B * (*coarse_data[1])(i);
                lower = std::min(lower, (*coarse_data[0])(i));
                upper = std::max(upper, (*coarse_data[0])(i));
            }
        for (int k = 0; k < N_FIELDS; ++k) refine_op.refine(fine, coarse, u_idxs[k], u_idxs[k], fine_ghost_box, ratio);

        double defect = 0.0, scale = 0.0;
        int n_sides = 0;
        for (int axis = 0; axis < NDIM; ++axis)
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(fine_ghost_box, axis)); it; it++)
            {
                const SideIndex<NDIM> i(it(), axis, SideIndex<NDIM>::Lower);
                const double u0 = (*fine_data[0])(i);
                if (u0 < lower - TOLERANCE || u0 > upper + TOLERANCE)
                    TBOX_ERROR("Refinement produced a value outside the range of the coarse data.\n");
                defect = std::max(defect, std::abs((*fine_data[2])(i)-A * u0 - B * (*fine_data[1])(i)));
                scale = std::max(scale, std::abs((*fine_data[2])(i)));
                ++n_sides;
            }
        if (!(defect <= TOLERANCE * scale)) TBOX_ERROR("Refinement is not linear in the data.\n");
        plog << "random ratio " << ratio << ": refined " << n_sides << " sides\n";
    }
}

// Fill the ghost cells of side-centered data on a hierarchy that has coarse-fine interfaces and either physical
// boundaries (the field is linear, so the fill is exact) or periodic boundaries (the field is smooth).
void
test_linear_fill(Pointer<PatchHierarchy<NDIM>> patch_hierarchy, const int u_idx, const bool smooth)
{
    constexpr double SENTINEL = 1.e30;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln) patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_idx, 0.0);

    auto exact = [&](const double* const x, const int axis)
    {
        if (!smooth)
        {
            double value = axis + 1.0;
            for (int d = 0; d < NDIM; ++d) value += (axis + 2.0) * (d + 1.0) * x[d];
            return value;
        }
        double value = std::sin(2.0 * M_PI * x[axis] + 0.3) * std::cos(2.0 * M_PI * x[(axis + 1) % NDIM] + 0.5);
        if (NDIM == 3) value *= std::cos(2.0 * M_PI * x[(axis + 2) % NDIM] + 0.1);
        return value;
    };
    auto side_coordinates = [](const Patch<NDIM>& patch, const hier::Index<NDIM>& i, const int axis, double* const x)
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
        for (int d = 0; d < NDIM; ++d)
            x[d] =
                pgeom->getXLower()[d] + pgeom->getDx()[d] * (i(d) - patch.getBox().lower(d) + (d == axis ? 0.0 : 0.5));
    };

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            u_data->fillAll(SENTINEL);
            for (int axis = 0; axis < NDIM; ++axis)
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                {
                    double x[NDIM];
                    side_coordinates(*patch, it(), axis, x);
                    (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = exact(x, axis);
                }
        }
    }

    HierarchyGhostCellInterpolation fill_op;
    fill_op.initializeOperatorState(HierarchyGhostCellInterpolation::InterpolationTransactionComponent(
                                        u_idx, "LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr),
                                    patch_hierarchy);
    fill_op.fillData(0.0);

    for (int ln = 1; ln <= finest_ln; ++ln)
    {
        double max_error = 0.0;
        int n_ghosts = 0;
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> interior = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); it; it++)
                {
                    if (interior.contains(it())) continue;
                    double x[NDIM];
                    side_coordinates(*patch, it(), axis, x);
                    const double value = (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower));
                    if (!std::isfinite(value) || value == SENTINEL) TBOX_ERROR("A ghost value was not filled.\n");
                    max_error = std::max(max_error, std::abs(value - exact(x, axis)));
                    ++n_ghosts;
                }
            }
        }
        if (!smooth && max_error > 1.0e-12) TBOX_ERROR("The fill of a linear field is not exact.\n");
        plog << "level " << ln << ": ghost sides = " << n_ghosts << ", max error = " << max_error << '\n';
    }

    for (int ln = 0; ln <= finest_ln; ++ln) patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(u_idx);
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
            IBTK_DISABLE_EXTRA_WARNINGS
            CartSideDoubleSpecializedLinearRefine refine_op;
            IBTK_ENABLE_EXTRA_WARNINGS
            test_exact_refine(refine_op);
            return 0;
        }
        if (input_db->getBoolWithDefault("test_linear_refine", false))
        {
            CartSideDoubleLinearRefine refine_op;
            test_exact_refine(refine_op);
            test_smooth_refine(refine_op);
            test_random_refine(refine_op);
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
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx);
        Pointer<SideVariable<NDIM, double>> exact_sc_var = new SideVariable<NDIM, double>("exact_sc");
        const int exact_sc_idx = var_db->registerVariableAndContext(exact_sc_var, ctx);
        const bool linear_fill = input_db->getBoolWithDefault("test_linear_fill", false);
        int u_fill_idx = IBTK::invalid_index;
        if (linear_fill)
        {
            grid_geometry->addSpatialRefineOperator(new CartSideDoubleLinearRefine());
            Pointer<SideVariable<NDIM, double>> u_fill_var = new SideVariable<NDIM, double>("u_fill");
            u_fill_idx = var_db->registerVariableAndContext(u_fill_var, ctx, IntVector<NDIM>(2));
        }
        // u_cc_var is only for plotting (and testing): uncomment if output is desired
// #define DO_PLOT
#ifdef DO_PLOT
        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx);
        Pointer<CellVariable<NDIM, double>> exact_cc_var = new CellVariable<NDIM, double>("exact_cc", NDIM);
        const int exact_cc_idx = var_db->registerVariableAndContext(exact_cc_var, ctx);
#endif

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = linear_fill ? 0 : std::numeric_limits<int>::max();
        int level_number = 0;
        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        if (linear_fill)
        {
            test_linear_fill(patch_hierarchy, u_fill_idx, input_db->getBool("smooth_field"));
            return 0;
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
