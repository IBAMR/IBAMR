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

// Check PETScMFFDJacobianOperator::formJacobian() and apply() when the operator is not attached to
// a nonlinear solver.

#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StokesSpecifications.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/muParserCartGridFunction.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
        Pointer<Database> input = app->getInputDatabase();

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

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
        const int f_u_idx = var_db->registerVariableAndContext(f_u_var, ctx, IntVector<NDIM>(1));
        const int f_p_idx = var_db->registerVariableAndContext(f_p_var, ctx, IntVector<NDIM>(1));

        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        level->allocatePatchData(u_idx, 0.0);
        level->allocatePatchData(p_idx, 0.0);
        level->allocatePatchData(f_u_idx, 0.0);
        level->allocatePatchData(f_p_idx, 0.0);

        HierarchyMathOps hier_math_ops("hier_math_ops", hierarchy);
        const int h_u_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int h_p_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", hierarchy, 0, 0);
        SAMRAIVectorReal<NDIM, double> f_vec("f", hierarchy, 0, 0);
        u_vec.addComponent(u_var, u_idx, h_u_idx);
        u_vec.addComponent(p_var, p_idx, h_p_idx);
        f_vec.addComponent(f_u_var, f_u_idx, h_u_idx);
        f_vec.addComponent(f_p_var, f_p_idx, h_p_idx);

        u_vec.setToScalar(0.0);
        muParserCartGridFunction u_fcn("u", app->getComponentDatabase("u"), geometry);
        u_fcn.setDataOnPatchHierarchy(u_idx, u_var, hierarchy, 0.0);
        f_vec.setToScalar(0.0);

        Pointer<StaggeredStokesOperator> stokes_op = new StaggeredStokesOperator("stokes_op", true);
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setDConstant(input->getDouble("D"));
        poisson_spec.setCConstant(input->getDouble("C"));
        stokes_op->setVelocityPoissonSpecifications(poisson_spec);
        stokes_op->initializeOperatorState(u_vec, f_vec);

        // The Stokes operator F is linear on a periodic domain, so F'(u) v = F(v). Compare the
        // finite-difference Jacobian action on v = u with F(u).
        Pointer<SAMRAIVectorReal<NDIM, double>> expected_vec = f_vec.cloneVector("expected");
        expected_vec->allocateVectorData();
        stokes_op->apply(u_vec, *expected_vec);

        PETScMFFDJacobianOperator mffd("mffd");
        mffd.setOperator(stokes_op);
        mffd.initializeOperatorState(u_vec, f_vec);
        mffd.formJacobian(u_vec);
        mffd.apply(u_vec, f_vec);

        const double action_scale = expected_vec->maxNorm();
        plog << "max norm of F(u) = " << action_scale << "\n";
        plog << "max norm of the finite-difference Jacobian action = " << f_vec.maxNorm() << "\n";

        // The two differ by rounding error, which the comparison of the output does not resolve.
        f_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false), expected_vec);
        if (!(f_vec.maxNorm() < 1.0e-5 * action_scale))
        {
            TBOX_ERROR("The finite-difference Jacobian action differs from F(u).\n");
        }

        mffd.deallocateOperatorState();
        stokes_op->deallocateOperatorState();
        expected_vec->freeVectorComponents();

        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(p_idx);
        level->deallocatePatchData(f_u_idx);
        level->deallocatePatchData(f_p_idx);
    }
    return 0;
} // main
