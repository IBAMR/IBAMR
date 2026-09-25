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

// PETScMFFDJacobianOperator::formJacobian() works both inside a Newton solve (where PETSc's SNES
// supplies the base state and function value) and standalone (where the operator must evaluate the
// function itself). This test exercises only the standalone path, since implicit_stokes_ib_solver_
// components_01's Jacobian fixtures already exercise the SNES-driven path.

#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StokesSpecifications.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
        Pointer<Database> input = app->getInputDatabase();
        const std::string test_case = input->getStringWithDefault("test_case", "success");
        if (test_case != "success")
        {
            Pointer<Logger::Appender> abort_appender = new TestAppender();
            Logger::getInstance()->setAbortAppender(abort_appender);
        }

        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tagging = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), tagging, boxes, load_balancer);
        gridding->makeCoarsestLevel(hierarchy, 0.0);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");
        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
        Pointer<SideVariable<NDIM, double>> f_u_var = new SideVariable<NDIM, double>("f_u");
        Pointer<CellVariable<NDIM, double>> f_p_var = new CellVariable<NDIM, double>("f_p");
        Pointer<SideVariable<NDIM, double>> e_u_var = new SideVariable<NDIM, double>("e_u");
        Pointer<CellVariable<NDIM, double>> e_p_var = new CellVariable<NDIM, double>("e_p");
        Pointer<SideVariable<NDIM, double>> exp_u_var = new SideVariable<NDIM, double>("exp_u");
        Pointer<CellVariable<NDIM, double>> exp_p_var = new CellVariable<NDIM, double>("exp_p");

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
        const int f_u_idx = var_db->registerVariableAndContext(f_u_var, ctx, IntVector<NDIM>(1));
        const int f_p_idx = var_db->registerVariableAndContext(f_p_var, ctx, IntVector<NDIM>(1));
        const int e_u_idx = var_db->registerVariableAndContext(e_u_var, ctx, IntVector<NDIM>(1));
        const int e_p_idx = var_db->registerVariableAndContext(e_p_var, ctx, IntVector<NDIM>(1));
        const int exp_u_idx = var_db->registerVariableAndContext(exp_u_var, ctx, IntVector<NDIM>(1));
        const int exp_p_idx = var_db->registerVariableAndContext(exp_p_var, ctx, IntVector<NDIM>(1));

        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        level->allocatePatchData(u_idx, 0.0);
        level->allocatePatchData(p_idx, 0.0);
        level->allocatePatchData(f_u_idx, 0.0);
        level->allocatePatchData(f_p_idx, 0.0);
        level->allocatePatchData(e_u_idx, 0.0);
        level->allocatePatchData(e_p_idx, 0.0);
        level->allocatePatchData(exp_u_idx, 0.0);
        level->allocatePatchData(exp_p_idx, 0.0);

        HierarchyMathOps hier_math_ops("hier_math_ops", hierarchy);
        const int h_u_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int h_p_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> f_vec("f", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> e_vec("e", hierarchy, 0, 0);
        u_vec.addComponent(u_var, u_idx, h_u_idx);
        u_vec.addComponent(p_var, p_idx, h_p_idx);
        f_vec.addComponent(f_u_var, f_u_idx, h_u_idx);
        f_vec.addComponent(f_p_var, f_p_idx, h_p_idx);
        e_vec.addComponent(e_u_var, e_u_idx, h_u_idx);
        e_vec.addComponent(e_p_var, e_p_idx, h_p_idx);

        u_vec.setToScalar(0.0);
        muParserCartGridFunction u_fcn("u", app->getComponentDatabase("u"), geometry);
        u_fcn.setDataOnPatchHierarchy(u_idx, u_var, hierarchy, 0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        Pointer<StaggeredStokesOperator> stokes_op = new StaggeredStokesOperator("stokes_op", true);
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setDConstant(input->getDouble("D"));
        poisson_spec.setCConstant(input->getDouble("C"));
        stokes_op->setVelocityPoissonSpecifications(poisson_spec);
        stokes_op->initializeOperatorState(u_vec, f_vec);

        // stokes_op is linear (no boundary offset, since the domain is periodic), so its exact
        // Jacobian-vector product at any base state equals the operator applied directly to the
        // direction: F'(u)v = F(v). A finite-difference approximation of F'(u)v that is close to
        // F(v) is therefore evidence both that the base function value formJacobian() computed
        // standalone is correct (a wrong base value would throw off the difference quotient) and
        // that the resulting action is sane.
        SAMRAIVectorReal<NDIM, double> expected_vec("expected", hierarchy, 0, 0);
        expected_vec.addComponent(exp_u_var, exp_u_idx, h_u_idx);
        expected_vec.addComponent(exp_p_var, exp_p_idx, h_p_idx);
        stokes_op->apply(u_vec, expected_vec);

        PETScMFFDJacobianOperator mffd("mffd");
        if (test_case != "no_operator")
        {
            mffd.setOperator(stokes_op);
        }
        if (test_case != "uninitialized")
        {
            mffd.initializeOperatorState(u_vec, f_vec);
        }
        // Neither a base state nor a nonlinear solver is set: this is the standalone path.
        mffd.formJacobian(u_vec);
        mffd.apply(u_vec, f_vec);

        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double>>(&expected_vec, false));
        const double action_error = e_vec.maxNorm();
        const double action_scale = expected_vec.maxNorm();
        pout << "jacobian action relative error = " << action_error / action_scale << "\n";
        if (!std::isfinite(action_error) || !(action_error < 1.0e-5 * action_scale))
        {
            TBOX_ERROR(
                "The finite-difference Jacobian action computed outside a Newton solve does not match a "
                "direct application of the (linear) wrapped operator.\n");
        }

        mffd.deallocateOperatorState();
        stokes_op->deallocateOperatorState();

        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(p_idx);
        level->deallocatePatchData(f_u_idx);
        level->deallocatePatchData(f_p_idx);
        level->deallocatePatchData(e_u_idx);
        level->deallocatePatchData(e_p_idx);
        level->deallocatePatchData(exp_u_idx);
        level->deallocatePatchData(exp_p_idx);
    }
    return 0;
} // main
