// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/*
 * This test retains the legacy affine/non-RT RT0 checks and adds dedicated
 * input-driven coverage for nontrivial RT0 refine exactness and divergence
 * preservation.
 */

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
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <array>
#include <cmath>
#include <string>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

namespace
{
enum class ProfileType
{
    AFFINE,
    PIECEWISE_RT0,
    NONLINEAR,
    UNKNOWN
};

ProfileType
string_to_profile_type(const std::string& profile)
{
    ProfileType profile_type = ProfileType::UNKNOWN;
    if (profile == "affine") profile_type = ProfileType::AFFINE;
    if (profile == "piecewise_rt0") profile_type = ProfileType::PIECEWISE_RT0;
    if (profile == "nonlinear") profile_type = ProfileType::NONLINEAR;
    if (profile_type == ProfileType::UNKNOWN) TBOX_ERROR("Unknown profile_type = " << profile << "\n");
    return profile_type;
}

int
compute_piecewise_region(const double x)
{
    if (x < 1.0 / 3.0) return 0;
    if (x < 2.0 / 3.0) return 1;
    return 2;
}

double
compute_piecewise_rt0_value(const VectorNd& X, const int axis)
{
    int transverse_sum = 0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis) continue;
        transverse_sum += compute_piecewise_region(X[d]);
    }
    const int mod3 = transverse_sum % 3;
    const int mod5 = (2 * transverse_sum + axis) % 5;
    const double c0 = 0.1 * static_cast<double>(axis + 1) + 0.05 * static_cast<double>(mod3);
    const double c1 = 0.2 + 0.03 * static_cast<double>(mod5);
    return c0 + c1 * X[axis];
}

void
set_affine_side_field(Pointer<SideData<NDIM, double>> u_data,
                      Pointer<Patch<NDIM>> patch,
                      const std::array<double, NDIM>& coeffs)
{
    const Box<NDIM>& patch_box = patch->getBox();

    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const VectorNd X = IBTK::IndexUtilities::getSideCenter(*patch, i_s);
            (*u_data)(i_s) = coeffs[axis] * X[axis];
        }
    }
}

void
refine_level_side_data(const int coarse_level_n,
                       const int fine_level_n,
                       const int dst_data_idx,
                       const int src_data_idx,
                       Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarse_level_n);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(fine_level_n);
    TBOX_ASSERT(coarse_level->getNumberOfPatches() == 1);
    TBOX_ASSERT(fine_level->getNumberOfPatches() == 1);

    const IntVector<NDIM> ratio = fine_level->getRatioToCoarserLevel();
    Pointer<Patch<NDIM>> coarse_patch = coarse_level->getPatch(0);
    Pointer<Patch<NDIM>> fine_patch = fine_level->getPatch(0);
    const Box<NDIM> fine_box = fine_patch->getBox();

    IBTK::CartSideDoubleRT0Refine refine_op;
    refine_op.refine(*fine_patch, *coarse_patch, dst_data_idx, src_data_idx, fine_box, ratio);
}

double
compute_cell_divergence(Pointer<SideData<NDIM, double>> u_data, const CellIndex<NDIM>& cell_idx, const double* dx)
{
    double div = 0.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const SideIndex<NDIM> lower(cell_idx, axis, SideIndex<NDIM>::Lower);
        SideIndex<NDIM> upper(cell_idx, axis, SideIndex<NDIM>::Lower);
        upper(axis) += 1;
        div += ((*u_data)(upper) - (*u_data)(lower)) / dx[axis];
    }
    return div;
}

double
compute_nonlinear_value(const VectorNd& X, const int axis)
{
    const double pi = 3.14159265358979323846;
    double transverse_sum = 0.0;
    double transverse_prod = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis) continue;
        transverse_sum += X[d];
        transverse_prod *= (X[d] + 0.1 * static_cast<double>(d + 1));
    }
    return std::sin(2.0 * pi * X[axis]) + 0.35 * std::cos(pi * transverse_sum) + 0.1 * X[axis] * X[axis] +
           0.05 * transverse_prod;
}

void
set_piecewise_rt0_side_field(Pointer<SideData<NDIM, double>> u_data, Pointer<Patch<NDIM>> patch)
{
    const Box<NDIM>& patch_box = patch->getBox();

    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const VectorNd X = IBTK::IndexUtilities::getSideCenter(*patch, i_s);
            (*u_data)(i_s) = compute_piecewise_rt0_value(X, axis);
        }
    }
}

void
set_nonlinear_side_field(Pointer<SideData<NDIM, double>> u_data, Pointer<Patch<NDIM>> patch)
{
    const Box<NDIM>& patch_box = patch->getBox();

    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const VectorNd X = IBTK::IndexUtilities::getSideCenter(*patch, i_s);
            (*u_data)(i_s) = compute_nonlinear_value(X, axis);
        }
    }
}

void
set_test_profile_side_field(Pointer<SideData<NDIM, double>> u_data,
                            Pointer<Patch<NDIM>> patch,
                            const ProfileType profile_type,
                            const std::array<double, NDIM>& coeffs)
{
    switch (profile_type)
    {
    case ProfileType::AFFINE:
        set_affine_side_field(u_data, patch, coeffs);
        break;
    case ProfileType::PIECEWISE_RT0:
        set_piecewise_rt0_side_field(u_data, patch);
        break;
    case ProfileType::NONLINEAR:
        set_nonlinear_side_field(u_data, patch);
        break;
    case ProfileType::UNKNOWN:
    default:
        TBOX_ERROR("Unknown ProfileType encountered.\n");
    }
}

int
run_affine_exactness_test(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                          const int u_sc_idx,
                          const int exact_sc_idx,
                          const std::array<double, NDIM>& coeffs,
                          const ProfileType profile_type,
                          const double tol)
{
    const int coarse_level_n = 0;
    const int fine_level_n = 1;
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarse_level_n);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(fine_level_n);
    TBOX_ASSERT(coarse_level->getNumberOfPatches() == 1);
    TBOX_ASSERT(fine_level->getNumberOfPatches() == 1);

    Pointer<Patch<NDIM>> coarse_patch = coarse_level->getPatch(0);
    Pointer<Patch<NDIM>> fine_patch = fine_level->getPatch(0);
    Pointer<SideData<NDIM, double>> coarse_u_data = coarse_patch->getPatchData(u_sc_idx);
    Pointer<SideData<NDIM, double>> fine_u_data = fine_patch->getPatchData(u_sc_idx);
    Pointer<SideData<NDIM, double>> fine_exact_data = fine_patch->getPatchData(exact_sc_idx);

    coarse_u_data->fillAll(0.0);
    fine_u_data->fillAll(0.0);
    fine_exact_data->fillAll(0.0);
    set_test_profile_side_field(coarse_u_data, coarse_patch, profile_type, coeffs);
    set_test_profile_side_field(fine_exact_data, fine_patch, profile_type, coeffs);
    refine_level_side_data(coarse_level_n, fine_level_n, u_sc_idx, u_sc_idx, patch_hierarchy);

    const Box<NDIM>& fine_box = fine_patch->getBox();
    double max_interior_error = 0.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(fine_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            if (i_s(axis) == side_box.lower(axis) || i_s(axis) == side_box.upper(axis)) continue;
            const double err = std::abs((*fine_u_data)(i_s) - (*fine_exact_data)(i_s));
            max_interior_error = std::max(max_interior_error, err);
        }
    }

    pout << "rt0 affine exactness max interior side error: " << max_interior_error << '\n';
    const int test_failures = IBTK::abs_equal_eps(max_interior_error, 0.0, tol) ? 0 : 1;
    pout << "rt0 affine exactness test_failures: " << test_failures << '\n';
    return test_failures;
}

int
run_divergence_preservation_test(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                 const int u_sc_idx,
                                 const std::array<double, NDIM>& coeffs,
                                 const ProfileType profile_type,
                                 const double tol)
{
    const int coarse_level_n = 0;
    const int fine_level_n = 1;
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarse_level_n);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(fine_level_n);
    TBOX_ASSERT(coarse_level->getNumberOfPatches() == 1);
    TBOX_ASSERT(fine_level->getNumberOfPatches() == 1);

    Pointer<Patch<NDIM>> coarse_patch = coarse_level->getPatch(0);
    Pointer<Patch<NDIM>> fine_patch = fine_level->getPatch(0);
    Pointer<SideData<NDIM, double>> coarse_u_data = coarse_patch->getPatchData(u_sc_idx);
    Pointer<SideData<NDIM, double>> fine_u_data = fine_patch->getPatchData(u_sc_idx);
    coarse_u_data->fillAll(0.0);
    fine_u_data->fillAll(0.0);
    set_test_profile_side_field(coarse_u_data, coarse_patch, profile_type, coeffs);
    refine_level_side_data(coarse_level_n, fine_level_n, u_sc_idx, u_sc_idx, patch_hierarchy);

    Pointer<CartesianPatchGeometry<NDIM>> coarse_pgeom = coarse_patch->getPatchGeometry();
    Pointer<CartesianPatchGeometry<NDIM>> fine_pgeom = fine_patch->getPatchGeometry();
    const double* coarse_dx = coarse_pgeom->getDx();
    const double* fine_dx = fine_pgeom->getDx();

    const IntVector<NDIM> ratio = fine_level->getRatioToCoarserLevel();
    const Box<NDIM>& coarse_box = coarse_patch->getBox();
    double max_div_error = 0.0;
    int bad_count = 0;
    for (Box<NDIM>::Iterator b(coarse_box); b; b++)
    {
        const CellIndex<NDIM> coarse_idx = b();
        const double coarse_div = compute_cell_divergence(coarse_u_data, coarse_idx, coarse_dx);

        CellIndex<NDIM> fine_lower, fine_upper;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            fine_lower(axis) = ratio(axis) * coarse_idx(axis);
            fine_upper(axis) = fine_lower(axis) + ratio(axis) - 1;
        }
        const Box<NDIM> fine_children_box(fine_lower, fine_upper);
        for (Box<NDIM>::Iterator b_child(fine_children_box); b_child; b_child++)
        {
            const CellIndex<NDIM> fine_idx = b_child();
            const double fine_div = compute_cell_divergence(fine_u_data, fine_idx, fine_dx);
            const double err = std::abs(fine_div - coarse_div);
            max_div_error = std::max(max_div_error, err);
            if (!IBTK::abs_equal_eps(fine_div, coarse_div, tol)) ++bad_count;
        }
    }

    pout << "rt0 divergence preservation max error: " << max_div_error << '\n';
    pout << "rt0 divergence preservation bad_count: " << bad_count << '\n';
    const int test_failures = bad_count > 0 ? 1 : 0;
    pout << "rt0 divergence preservation test_failures: " << test_failures << '\n';
    return test_failures;
}
} // namespace

// Verify that we can correctly refine a piecewise linear solution with the
// RT0 refinement class. Since the RT0 element is a vector-valued element that
// is, on Cartesian grids,
//
//     RT0_K = (P^1(x) * P^0(y), P^0(x) * P^1(y))
//
// we expect refining a vector field which is in that space to have zero error.

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> test_db =
            input_db->keyExists("test") ? input_db->getDatabase("test") : Pointer<Database>(input_db, false);
        const std::string test_mode = test_db->getStringWithDefault("test_mode", "legacy");
        const std::string profile_string = test_db->getStringWithDefault("profile_type", "affine");
        const ProfileType profile_type = string_to_profile_type(profile_string);
        const double tol = test_db->getDoubleWithDefault("tol", 1.0e-12);
        std::array<double, NDIM> affine_coeffs;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            affine_coeffs[axis] = test_db->getDoubleWithDefault("affine_coefficient_" + std::to_string(axis), 1.0);
        }
        plog << "Input database:\n";
        input_db->printClassData(plog);

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

        if (test_mode == "legacy")
        {
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
                const Box<NDIM> patch_box_1 = level_1->getPatch(0)->getBox();

                const IntVector<NDIM> ratio = level_1->getRatioToCoarserLevel();
                IBTK::CartSideDoubleRT0Refine refine_op;
                refine_op.refine(*level_1->getPatch(0), *level_0->getPatch(0), u_sc_idx, u_sc_idx, patch_box_1, ratio);

                solv::SAMRAIVectorReal<NDIM, double> exact_vec("e", patch_hierarchy, coarse_level_n, fine_level_n);
                exact_vec.addComponent(exact_sc_var, exact_sc_idx);
                exact_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false),
                                   Pointer<SAMRAIVectorReal<NDIM, double>>(&exact_vec, false));

                pout << "test results for " << db_u_fcn_name << '\n';
                pout << "max norm of u_sc: " << u_vec.maxNorm() << '\n';
                pout << "max norm of exact - refined: " << std::abs(exact_vec.maxNorm()) << '\n';
            };

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
        }
        else if (test_mode == "affine_exactness")
        {
            const int test_failures =
                run_affine_exactness_test(patch_hierarchy, u_sc_idx, exact_sc_idx, affine_coeffs, profile_type, tol);
            if (test_failures != 0) return test_failures;
        }
        else if (test_mode == "divergence_preservation")
        {
            const int test_failures =
                run_divergence_preservation_test(patch_hierarchy, u_sc_idx, affine_coeffs, profile_type, tol);
            if (test_failures != 0) return test_failures;
        }
        else
        {
            TBOX_ERROR("Unknown test_mode = " << test_mode << "\n");
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
