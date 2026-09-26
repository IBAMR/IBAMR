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

// Check that resetting only the finest level of a StaggeredStokesIBLevelRelaxationFACOperator keeps the coarse solver
// on the retained coarsest level, i.e., that solving on the coarsest level after a partial reset works and gives the
// same answer as before the reset.

#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScMatUtilities.h>

#include <petscmat.h>
#include <petscvec.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>
#include <SAMRAI_config.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <ibamr/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "stokes_ib_fac_partial_reset.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Two-level periodic hierarchy.
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

        // Patch data: solution, right-hand side, and the DOF numbering used to build the IB operators. These must be
        // registered before the first level is created, since that fixes the maximum ghost width.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("ctx");
        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
        Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
        Pointer<CellVariable<NDIM, double>> g_var = new CellVariable<NDIM, double>("g");
        Pointer<SideVariable<NDIM, int>> u_dof_var = new SideVariable<NDIM, int>("u_dof");
        Pointer<CellVariable<NDIM, int>> p_dof_var = new CellVariable<NDIM, int>("p_dof");
        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
        const int f_idx = var_db->registerVariableAndContext(f_var, ctx, IntVector<NDIM>(1));
        const int g_idx = var_db->registerVariableAndContext(g_var, ctx, IntVector<NDIM>(1));
        // The IB_4 interpolation stencil reaches two cells beyond the patch.
        const int u_dof_idx = var_db->registerVariableAndContext(u_dof_var, ctx, IntVector<NDIM>(3));
        const int p_dof_idx = var_db->registerVariableAndContext(p_dof_var, ctx, IntVector<NDIM>(1));

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, /*tag_buffer*/ 1);
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        if (finest_ln != 1) TBOX_ERROR("expected a two-level patch hierarchy\n");
        Pointer<PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(finest_ln);

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int idx : { u_idx, p_idx, f_idx, g_idx, u_dof_idx, p_dof_idx })
            {
                level->allocatePatchData(idx, 0.0);
            }
        }

        // Interpolation operator J and force Jacobian A for a single Lagrangian point on the finest level.
        std::vector<int> num_dofs_per_proc;
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
            num_dofs_per_proc, u_dof_idx, p_dof_idx, finest_level);

        const int num_local_lag_dofs = IBTK_MPI::getRank() == 0 ? NDIM : 0;
        Vec X_vec;
        int ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_lag_dofs, PETSC_DETERMINE, &X_vec);
        IBTK_CHKERRQ(ierr);
        for (int d = 0; d < num_local_lag_dofs; ++d)
        {
            ierr = VecSetValue(X_vec, d, 0.43 + 0.05 * d, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(X_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(X_vec);
        IBTK_CHKERRQ(ierr);

        Mat J = nullptr;
        PETScMatUtilities::constructPatchLevelSCInterpOp(
            J, PETScMatUtilities::ib_4_interp_fcn, 4, X_vec, num_dofs_per_proc, u_dof_idx, finest_level);

        Mat A;
        ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                            num_local_lag_dofs,
                            num_local_lag_dofs,
                            PETSC_DETERMINE,
                            PETSC_DETERMINE,
                            1,
                            nullptr,
                            0,
                            nullptr,
                            &A);
        IBTK_CHKERRQ(ierr);
        for (int d = 0; d < num_local_lag_dofs; ++d)
        {
            ierr = MatSetValue(A, d, d, -input_db->getDouble("SPRING_STIFFNESS"), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);

        // Vectors.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        SAMRAIVectorReal<NDIM, double> sol_vec("sol", patch_hierarchy, 0, finest_ln);
        sol_vec.addComponent(u_var, u_idx, wgt_sc_idx);
        sol_vec.addComponent(p_var, p_idx, wgt_cc_idx);
        SAMRAIVectorReal<NDIM, double> rhs_vec("rhs", patch_hierarchy, 0, finest_ln);
        rhs_vec.addComponent(f_var, f_idx, wgt_sc_idx);
        rhs_vec.addComponent(g_var, g_idx, wgt_cc_idx);
        sol_vec.setToScalar(0.0);
        rhs_vec.setToScalar(0.0);

        // A smooth but non-constant velocity right-hand side on every level.
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(f_data->getBox(), axis)); b; b++)
                    {
                        const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                        (*f_data)(i) = std::sin(0.7 * i(0) + 1.3 * i(1) + axis);
                    }
                }
            }
        }

        // Set up the FAC operator.
        const double dt = input_db->getDouble("DT");
        PoissonSpecifications U_problem_coefs("U_problem_coefs");
        U_problem_coefs.setCConstant(input_db->getDouble("RHO") / dt);
        U_problem_coefs.setDConstant(-input_db->getDouble("MU"));
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();

        StaggeredStokesIBLevelRelaxationFACOperator fac_op(
            "fac_op", input_db->getDatabase("stokes_ib_precond_db"), "stokes_ib_pc_");
        fac_op.setVelocityPoissonSpecifications(U_problem_coefs);
        fac_op.setPhysicalBcCoefs(u_bc_coefs, nullptr);
        fac_op.setPhysicalBoundaryHelper(bc_helper);
        fac_op.setTimeInterval(0.0, dt);
        fac_op.setSolutionTime(dt);
        fac_op.setComponentsHaveNullSpace(false, true);
        fac_op.setIBTimeSteppingType(BACKWARD_EULER);
        fac_op.setIBForceJacobian(A);
        fac_op.setIBInterpOp(J);
        fac_op.initializeOperatorState(sol_vec, rhs_vec);

        // Solve on the coarsest level with the fully initialized operator.
        Pointer<SAMRAIVectorReal<NDIM, double>> error_before = sol_vec.cloneVector("error_before");
        error_before->allocateVectorData();
        error_before->setToScalar(0.0);
        fac_op.solveCoarsestLevel(*error_before, rhs_vec, 0);

        // Reset only the finest level and solve on the coarsest level again.
        fac_op.setResetLevels(finest_ln, finest_ln);
        fac_op.initializeOperatorState(sol_vec, rhs_vec);
        if (!fac_op.getStaggeredStokesPETScLevelSolver(0))
        {
            TBOX_ERROR("the coarsest level lost its solver after resetting only the finest level\n");
        }

        Pointer<SAMRAIVectorReal<NDIM, double>> error_after = sol_vec.cloneVector("error_after");
        error_after->allocateVectorData();
        error_after->setToScalar(0.0);
        fac_op.solveCoarsestLevel(*error_after, rhs_vec, 0);

        // Compare the level 0 solutions without control volume weights, which vanish where level 1 covers level 0.
        HierarchySideDataOpsReal<NDIM, double> u_ops(patch_hierarchy, 0, 0);
        HierarchyCellDataOpsReal<NDIM, double> p_ops(patch_hierarchy, 0, 0);
        const int u_before_idx = error_before->getComponentDescriptorIndex(0);
        const int p_before_idx = error_before->getComponentDescriptorIndex(1);
        const int u_after_idx = error_after->getComponentDescriptorIndex(0);
        const int p_after_idx = error_after->getComponentDescriptorIndex(1);
        const double norm_before = std::max(u_ops.maxNorm(u_before_idx, -1), p_ops.maxNorm(p_before_idx, -1));
        u_ops.axpy(u_after_idx, -1.0, u_before_idx, u_after_idx, false);
        p_ops.axpy(p_after_idx, -1.0, p_before_idx, p_after_idx, false);
        const double norm_diff = std::max(u_ops.maxNorm(u_after_idx, -1), p_ops.maxNorm(p_after_idx, -1));
        if (!(norm_before > 0.0))
        {
            TBOX_ERROR("the coarsest level solve returned a trivial solution\n");
        }
        if (!(norm_diff <= 1.0e-10 * norm_before))
        {
            TBOX_ERROR("resetting only the finest level changed the coarsest level solve\n");
        }

        error_before->freeVectorComponents();
        error_after->freeVectorComponents();
        fac_op.deallocateOperatorState();
        ierr = VecDestroy(&X_vec);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&A);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&J);
        IBTK_CHKERRQ(ierr);
    } // cleanup dynamically allocated objects prior to shutdown
} // main
