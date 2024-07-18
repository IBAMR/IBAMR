// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/PETScKrylovStaggeredStokesSolver.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StaggeredStokesSolverManager.h>
#include <ibamr/StokesSpecifications.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/

// The following program is used to test out whether or not the StokesOperator is performing correctly.
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "sc_poisson.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.

        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);
        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
                                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                    error_detector,
                                                    box_generator,
                                                    load_balancer);

        // Setup the Boundary Conditions
        const IntVectorNd& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategyNd*> u_bc_coefs(NDIM);
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
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        // State variables: Velocity and pressure. Need to get both cell sides and centers
        // since we are using the MAC scheme.
        SAMRAIPointer<SideVariableNd<double> > u_sc_var = make_samrai_shared<SideVariableNd<double> >("u_sc");
        SAMRAIPointer<CellVariableNd<double> > p_cc_var = make_samrai_shared<CellVariableNd<double> >("p_cc");

        // Results of operator "forces" and "divergence"
        SAMRAIPointer<SideVariableNd<double> > f_sc_var = make_samrai_shared<SideVariableNd<double> >("f_sc");
        SAMRAIPointer<CellVariableNd<double> > f_cc_var = make_samrai_shared<CellVariableNd<double> >("f_cc");

        // Error terms.
        SAMRAIPointer<SideVariableNd<double> > e_sc_var = make_samrai_shared<SideVariableNd<double> >("e_sc");
        SAMRAIPointer<CellVariableNd<double> > e_cc_var = make_samrai_shared<CellVariableNd<double> >("e_cc");

        // Register patch data indices...
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVectorNd(1));
        const int p_cc_idx = var_db->registerVariableAndContext(p_cc_var, ctx, IntVectorNd(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVectorNd(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVectorNd(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVectorNd(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVectorNd(1));

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
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
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

        SAMRAIVectorRealNd<double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorRealNd<double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorRealNd<double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorRealNd<double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

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
        if (periodic_shift.min() > 0)
        {
            StaggeredStokesOperator stokes_op("stokes_op", true);

            stokes_op.setVelocityPoissonSpecifications(poisson_spec);
            stokes_op.initializeOperatorState(u_vec, f_vec);
            // Apply the operator
            stokes_op.apply(u_vec, f_vec);
        }
        else
        {
            auto bc_helper = make_samrai_shared<StaggeredStokesPhysicalBoundaryHelper>();
            bc_helper->cacheBcCoefData(u_bc_coefs, 0.0, patch_hierarchy);
            bc_helper->copyDataAtDirichletBoundaries(e_sc_idx, u_sc_idx);
            // Setup the stokes operator
            StaggeredStokesOperator stokes_op("stokes_op", false, input_db);
            stokes_op.setPhysicalBcCoefs(u_bc_coefs, nullptr);
            stokes_op.setVelocityPoissonSpecifications(poisson_spec);
            stokes_op.setPhysicalBoundaryHelper(bc_helper);
            stokes_op.initializeOperatorState(u_vec, f_vec);

            // Apply the operator
            stokes_op.apply(u_vec, f_vec);
        }

        // Compute error and print error norms.
        e_vec.subtract(SAMRAIPointer<SAMRAIVectorRealNd<double> >(&f_vec, false),
                       SAMRAIPointer<SAMRAIVectorRealNd<double> >(&e_vec, false));
        // print out the errors in each norm
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Deallocate level data
        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(u_sc_idx);
            level->deallocatePatchData(f_sc_idx);
            level->deallocatePatchData(e_sc_idx);
            level->deallocatePatchData(p_cc_idx);
            level->deallocatePatchData(f_cc_idx);
            level->deallocatePatchData(e_cc_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
