// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/StaggeredStokesFACPreconditioner.h>
#include <ibamr/StaggeredStokesLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Logger.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include <ibamr/app_namespaces.h>

namespace
{
// std::max keeps its first argument when the second is NaN, so a nonfinite
// difference must make the accumulated error infinite instead.
double
accumulate_error(const double error, const double difference)
{
    return std::isfinite(difference) ? std::max(error, difference) : std::numeric_limits<double>::infinity();
}

// Print a difference that must vanish and abort if it does not.
void
report_vanishing_difference(const std::string& name, const double difference)
{
    plog << name << " = " << difference << '\n';
    if (!(difference <= 1.0e-12))
    {
        TBOX_ERROR("FAC restriction test:\n  " << name << " = " << difference << " does not vanish.\n");
    }
}

// Check the FAC restriction of the RT0 velocity between the two levels against the RT0 prolongation matrix P, for fine
// data whose ghost cells hold garbage. On the faces with a complete stencil it is the adjoint of the prolongation with
// respect to the cell-volume inner products, (1 / prod(ratio)) P^T. On the faces of a coarse-fine interface, whose
// stencil is truncated, it is the same adjoint; there the fine faces that do not exist contribute nothing. On the
// faces of a physical boundary it normalizes the weights of the fine faces inside the domain to sum to one, as the
// matrix restriction of PETScMatUtilities::constructRestrictionScalingOp() does everywhere. The coarse right-hand side
// is zero, so that the restricted level holds only the transferred fine data.
void
check_rt0_restriction_matrix(Pointer<PatchHierarchy<NDIM>> hierarchy,
                             StaggeredStokesFACPreconditioner& fac,
                             StaggeredStokesLevelRelaxationFACOperator& fac_op,
                             SAMRAIVectorReal<NDIM, double>& u_vec,
                             SAMRAIVectorReal<NDIM, double>& f_vec,
                             SAMRAIVectorReal<NDIM, double>& r_vec,
                             const bool dirichlet_boundary)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("adjoint_check");
    Pointer<SideVariable<NDIM, int>> dof_var = new SideVariable<NDIM, int>("adjoint_check_dof");
    const int dof_idx = var_db->registerVariableAndContext(dof_var, ctx, IntVector<NDIM>(1));
    Pointer<SideVariable<NDIM, double>> mask_var = new SideVariable<NDIM, double>("adjoint_check_mask");
    const int mask_idx = var_db->registerVariableAndContext(mask_var, ctx, IntVector<NDIM>(0));
    Pointer<PatchLevel<NDIM>> coarse = hierarchy->getPatchLevel(0);
    Pointer<PatchLevel<NDIM>> fine = hierarchy->getPatchLevel(1);
    std::vector<int> counts[2];
    for (int ln = 0; ln < 2; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(dof_idx);
        PETScVecUtilities::constructPatchLevelDOFIndices(counts[ln], dof_idx, level);
    }
    AO ordering = nullptr;
    PETScVecUtilities::constructPatchLevelAO(ordering, counts[0], dof_idx, coarse, 0);
    Mat P = nullptr;
    PETScMatUtilities::constructProlongationOp(P, "RT0", dof_idx, counts[1], counts[0], fine, coarse, ordering, 0);
    Vec fine_vec = nullptr, coarse_vec = nullptr, expected = nullptr, restricted = nullptr, matrix_expected = nullptr,
        weight = nullptr, mask = nullptr, scale = nullptr;
    int ierr = MatCreateVecs(P, &coarse_vec, &fine_vec);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &expected, &restricted, &matrix_expected, &weight, &mask })
    {
        ierr = VecDuplicate(coarse_vec, v);
        IBTK_CHKERRQ(ierr);
    }
    PETScMatUtilities::constructRestrictionScalingOp(P, scale);

    const int f_sc_idx = f_vec.getComponentDescriptorIndex(0);
    const int r_sc_idx = r_vec.getComponentDescriptorIndex(0);
    f_vec.setToScalar(0.0, false);
    r_vec.setToScalar(0.0, false);
    // With Dirichlet conditions the residual vanishes on the faces of a physical boundary, where the homogeneous
    // ghost fill of the restriction sets it to zero. Other conditions leave these values alone, so they stay nonzero.
    const Box<NDIM>& fine_domain = fine->getPhysicalDomain()[0];
    const IntVector<NDIM> fine_periodic_shift = fine->getGridGeometry()->getPeriodicShift(fine->getRatio());
    for (PatchLevel<NDIM>::Iterator p(fine); p; p++)
    {
        Pointer<SideData<NDIM, double>> data = fine->getPatch(p())->getPatchData(f_sc_idx);
        const Box<NDIM>& patch_box = fine->getPatch(p())->getBox();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(data->getGhostBox(), axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                const bool interior = SideGeometry<NDIM>::toSideBox(patch_box, axis).contains(b());
                const bool physical_boundary =
                    fine_periodic_shift(axis) == 0 &&
                    (b()(axis) == fine_domain.lower(axis) || b()(axis) == fine_domain.upper(axis) + 1);
                (*data)(i) = !interior                               ? 1.0e3 :
                             physical_boundary && dirichlet_boundary ? 0.0 :
                                                                       1.0 + (7 * b()(0) + 13 * b()(1) + 5 * axis) % 11;
            }
        }
    }
    PETScVecUtilities::copyToPatchLevelVec(fine_vec, f_sc_idx, dof_idx, fine);
    ierr = MatMultTranspose(P, fine_vec, expected);
    IBTK_CHKERRQ(ierr);
    double ratio_product = 1.0;
    for (int d = 0; d < NDIM; ++d) ratio_product *= fine->getRatioToCoarserLevel()(d);
    ierr = VecScale(expected, 1.0 / ratio_product);
    IBTK_CHKERRQ(ierr);

    // The expected restriction is the matrix restriction, except on the faces of a coarse-fine interface, where the
    // matrix restriction normalizes the truncated stencil and the adjoint does not. Those faces have a stencil weight
    // P^T 1 / prod(ratio) below one and are not on a physical boundary.
    ierr = VecPointwiseMult(matrix_expected, expected, scale);
    IBTK_CHKERRQ(ierr);
    ierr = VecScale(matrix_expected, ratio_product);
    IBTK_CHKERRQ(ierr);
    Vec ones = nullptr;
    ierr = VecDuplicate(fine_vec, &ones);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(ones, 1.0);
    IBTK_CHKERRQ(ierr);
    ierr = MatMultTranspose(P, ones, weight);
    IBTK_CHKERRQ(ierr);
    ierr = VecScale(weight, 1.0 / ratio_product);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&ones);
    IBTK_CHKERRQ(ierr);
    coarse->allocatePatchData(mask_idx);
    PETScVecUtilities::copyFromPatchLevelVec(weight, mask_idx, dof_idx, coarse, nullptr, nullptr);
    const Box<NDIM>& domain = coarse->getPhysicalDomain()[0];
    const IntVector<NDIM> periodic_shift = coarse->getGridGeometry()->getPeriodicShift(coarse->getRatio());
    for (PatchLevel<NDIM>::Iterator p(coarse); p; p++)
    {
        Pointer<Patch<NDIM>> patch = coarse->getPatch(p());
        Pointer<SideData<NDIM, double>> data = patch->getPatchData(mask_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
            {
                const bool physical_boundary = periodic_shift(axis) == 0 && (si()(axis) == domain.lower(axis) ||
                                                                             si()(axis) == domain.upper(axis) + 1);
                const double w = (*data)(si());
                (*data)(si()) = (w > 1.0e-12 && w < 1.0 - 1.0e-12 && !physical_boundary) ? 1.0 : 0.0;
            }
        }
    }
    PETScVecUtilities::copyToPatchLevelVec(mask, mask_idx, dof_idx, coarse);
    coarse->deallocatePatchData(mask_idx);
    // expected = matrix_expected + mask * (adjoint - matrix_expected)
    ierr = VecAXPY(expected, -1.0, matrix_expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecPointwiseMult(expected, expected, mask);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(expected, 1.0, matrix_expected);
    IBTK_CHKERRQ(ierr);

    fac.initializeSolverState(u_vec, f_vec);
    fac_op.restrictResidual(f_vec, r_vec, 0);
    fac.deallocateSolverState();
    PETScVecUtilities::copyToPatchLevelVec(restricted, r_sc_idx, dof_idx, coarse);
    ierr = VecAXPY(restricted, -1.0, expected);
    IBTK_CHKERRQ(ierr);
    double error = 0.0, norm = 0.0;
    ierr = VecNorm(restricted, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(expected, NORM_INFINITY, &norm);
    IBTK_CHKERRQ(ierr);
    if (!std::isfinite(error) || !std::isfinite(norm) || norm < 1.0 || error > 1.0e-12 * norm)
    {
        TBOX_ERROR("FAC RT0 restriction differs from the prolongation matrix restriction: error = "
                   << error << ", norm = " << norm << "\n");
    }
    plog << "restriction_matrix_check = passed\n";

    for (Vec* v : { &fine_vec, &coarse_vec, &expected, &restricted, &matrix_expected, &weight, &mask, &scale })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&P);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&ordering);
    IBTK_CHKERRQ(ierr);
    for (int ln = 0; ln < 2; ++ln) hierarchy->getPatchLevel(ln)->deallocatePatchData(dof_idx);
}

// The largest difference between two side-centered and two cell-centered data sets on a level,
// over their full ghost boxes.
double
max_data_difference(Pointer<PatchHierarchy<NDIM>> hierarchy,
                    const int level_number,
                    const int sc_idx_a,
                    const int sc_idx_b,
                    const int cc_idx_a,
                    const int cc_idx_b)
{
    double difference = 0.0;
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> sc_a = patch->getPatchData(sc_idx_a);
        Pointer<SideData<NDIM, double>> sc_b = patch->getPatchData(sc_idx_b);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(sc_a->getGhostBox(), axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                difference = accumulate_error(difference, std::abs((*sc_a)(i) - (*sc_b)(i)));
            }
        }
        Pointer<CellData<NDIM, double>> cc_a = patch->getPatchData(cc_idx_a);
        Pointer<CellData<NDIM, double>> cc_b = patch->getPatchData(cc_idx_b);
        for (Box<NDIM>::Iterator b(cc_a->getGhostBox()); b; b++)
        {
            difference = accumulate_error(difference, std::abs((*cc_a)(b(), 0) - (*cc_b)(b(), 0)));
        }
    }
    return IBTK_MPI::maxReduction(difference);
}
} // namespace

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/

// The following program tests the FAC restriction of the RT0 velocity and the invariants of the FAC solve on a
// two-level hierarchy.
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "sc_poisson.log");
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

        // Setup the Boundary Conditions
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // State variables: Velocity and pressure. Need to get both cell sides and centers
        // since we are using the MAC scheme.
        Pointer<SideVariable<NDIM, double>> u_sc_var = new SideVariable<NDIM, double>("u_sc");
        Pointer<CellVariable<NDIM, double>> p_cc_var = new CellVariable<NDIM, double>("p_cc");

        // Results of operator "forces" and "divergence"
        Pointer<SideVariable<NDIM, double>> f_sc_var = new SideVariable<NDIM, double>("f_sc");
        Pointer<CellVariable<NDIM, double>> f_cc_var = new CellVariable<NDIM, double>("f_cc");

        // Error terms.
        Pointer<SideVariable<NDIM, double>> e_sc_var = new SideVariable<NDIM, double>("e_sc");
        Pointer<CellVariable<NDIM, double>> e_cc_var = new CellVariable<NDIM, double>("e_cc");

        // Register patch data indices...
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(1));
        const int p_cc_idx = var_db->registerVariableAndContext(p_cc_var, ctx, IntVector<NDIM>(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVector<NDIM>(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));

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

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(f_sc_idx, 0.0);
            level->allocatePatchData(e_sc_idx, 0.0);
            level->allocatePatchData(p_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_sc_var, u_sc_idx, h_sc_idx);
        u_vec.addComponent(p_cc_var, p_cc_idx, h_cc_idx);
        f_vec.addComponent(f_sc_var, f_sc_idx, h_sc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_sc_var, e_sc_idx, h_sc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(0.0);

        // Setup velocity and pressures functions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction p_fcn("p", app_initializer->getComponentDatabase("p"), grid_geometry);

        // Setup exact solution functions
        muParserCartGridFunction f_u_fcn("f_u", app_initializer->getComponentDatabase("f_u"), grid_geometry);
        muParserCartGridFunction f_p_fcn("f_p", app_initializer->getComponentDatabase("f_p"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(p_cc_idx, p_cc_var, patch_hierarchy, 0.0);

        f_u_fcn.setDataOnPatchHierarchy(e_sc_idx, e_sc_var, patch_hierarchy, 0.0);
        f_p_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);

        // Setup stokes poisson specifications
        PoissonSpecifications poisson_spec("poisson_spec");
        const double D = input_db->getDouble("D");
        const double C = input_db->getDouble("C");
        poisson_spec.setDConstant(D);
        poisson_spec.setCConstant(C);
        // At refinement ratio two, cubic coarsening falls back to conservative coarsening
        // with a warning. Disable warning output to keep source-location text out of the
        // regression output.
        Logger::getInstance()->setWarning(false);
        Pointer<Database> fac_db = input_db->getDatabase("FACPreconditioner");
        Pointer<StaggeredStokesLevelRelaxationFACOperator> fac_op =
            new StaggeredStokesLevelRelaxationFACOperator("fac_op", fac_db, "fac_");
        StaggeredStokesFACPreconditioner fac("fac", fac_op, fac_db, "fac_");
        fac.setVelocityPoissonSpecifications(poisson_spec);
        fac.setComponentsHaveNullSpace(false, true);
        fac.setPhysicalBcCoefs(u_bc_coefs, nullptr);
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
        bc_helper->cacheBcCoefData(u_bc_coefs, 1.0, patch_hierarchy);
        fac.setPhysicalBoundaryHelper(bc_helper);
        fac.setTimeInterval(0.0, 1.0);
        fac.setSolutionTime(1.0);

        SAMRAIScopedVectorDuplicate<double> first_rhs_storage(f_vec);
        Pointer<SAMRAIVectorReal<NDIM, double>> first_rhs = first_rhs_storage;
        for (int trial = 0; trial < 2; ++trial)
        {
            // The two right-hand sides differ only in their fine velocity ghosts.
            f_vec.setToScalar(0.0, false);
            HierarchySideDataOpsReal<NDIM, double> fine_ops(patch_hierarchy, 1, 1);
            fine_ops.setToScalar(f_sc_idx, trial == 0 ? -7.0 : 11.0, false);
            f_u_fcn.setDataOnPatchHierarchy(f_sc_idx, f_sc_var, patch_hierarchy, 0.0);
            f_p_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
            // Copy the fine-level right-hand side with its ghost cells.
            SAMRAIScopedVectorDuplicate<double> rhs_before_storage(f_vec);
            Pointer<SAMRAIVectorReal<NDIM, double>> rhs_before = rhs_before_storage;
            static const bool interior_only = false;
            fine_ops.copyData(rhs_before->getComponentDescriptorIndex(0), f_sc_idx, interior_only);
            HierarchyCellDataOpsReal<NDIM, double> fine_cc_ops(patch_hierarchy, 1, 1);
            fine_cc_ops.copyData(rhs_before->getComponentDescriptorIndex(1), f_cc_idx, interior_only);
            const bool solved = fac.solveSystem(u_vec, f_vec);
            // The solve must leave the caller's fine-level right-hand side unchanged, including its ghost cells.
            report_vanishing_difference("fine_rhs_change",
                                        max_data_difference(patch_hierarchy,
                                                            1,
                                                            rhs_before->getComponentDescriptorIndex(0),
                                                            f_sc_idx,
                                                            rhs_before->getComponentDescriptorIndex(1),
                                                            f_cc_idx));
            const double correction_norm = u_vec.L2Norm();
            if (!solved || !std::isfinite(correction_norm) || correction_norm <= 1.0e-12)
            {
                TBOX_ERROR("FAC produced an invalid correction\n");
            }
            if (trial == 0)
            {
                first_rhs->copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false));
                e_vec.copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false));
            }
        }

        // FAC without presmoothing overwrites the coarse RHS. Compare coarse interiors
        // without hierarchy weights, so that differences in regions covered by the fine
        // level are included.
        if (fac.getNumPreSmoothingSweeps() == 0)
        {
            first_rhs->subtract(first_rhs, Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false));
            HierarchySideDataOpsReal<NDIM, double> coarse_ops(patch_hierarchy, 0, 0);
            report_vanishing_difference("coarse_rhs_difference",
                                        coarse_ops.maxNorm(first_rhs->getComponentDescriptorIndex(0)));
        }
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false));
        report_vanishing_difference("correction_difference", e_vec.maxNorm());
        if (input_db->getBoolWithDefault("check_restriction_matrix", false))
        {
            check_rt0_restriction_matrix(patch_hierarchy,
                                         fac,
                                         *fac_op,
                                         u_vec,
                                         f_vec,
                                         e_vec,
                                         input_db->getBoolWithDefault("dirichlet_boundary", true));
        }

        // Deallocate level data
        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(u_sc_idx);
            level->deallocatePatchData(f_sc_idx);
            level->deallocatePatchData(e_sc_idx);
            level->deallocatePatchData(p_cc_idx);
            level->deallocatePatchData(f_cc_idx);
            level->deallocatePatchData(e_cc_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
