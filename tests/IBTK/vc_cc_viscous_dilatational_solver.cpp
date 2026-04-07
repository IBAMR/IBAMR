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

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScKrylovPoissonSolver.h>
#include <ibtk/ProblemSpecification.h>
#include <ibtk/VCCCViscousDilatationalOperator.h>
#include <ibtk/VCCCViscousDilatationalPETScLevelSolver.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

Pointer<PoissonSolver>
allocate_vc_cc_velocity_krylov_solver(const std::string& solver_object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                      const std::string& solver_default_options_prefix)
{
    Pointer<PETScKrylovPoissonSolver> krylov_solver =
        new PETScKrylovPoissonSolver(solver_object_name, solver_input_db, solver_default_options_prefix);
    krylov_solver->setOperator(new VCCCViscousDilatationalOperator(solver_object_name + "::vc_cc_viscous_operator"));
    return krylov_solver;
}

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

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "vc_cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        Pointer<CellVariable<NDIM, double> > f_cc_var = new CellVariable<NDIM, double>("f_cc", NDIM);
        Pointer<CellVariable<NDIM, double> > e_cc_var = new CellVariable<NDIM, double>("e_cc", NDIM);
        Pointer<CellVariable<NDIM, double> > r_cc_var = new CellVariable<NDIM, double>("r_cc", NDIM);

        Pointer<CellVariable<NDIM, double> > mu_cc_var = new CellVariable<NDIM, double>("mu_cc");
        const int mu_cc_idx = var_db->registerVariableAndContext(mu_cc_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > lambda_cc_var = new CellVariable<NDIM, double>("lambda_cc");
        const int lambda_cc_idx = var_db->registerVariableAndContext(lambda_cc_var, ctx, IntVector<NDIM>(1));

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));
        const int r_cc_idx = var_db->registerVariableAndContext(r_cc_var, ctx, IntVector<NDIM>(1));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + std::to_string(d), "SCALAR", u_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "VECTOR", f_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(f_cc_var->getName() + std::to_string(d), "SCALAR", f_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(e_cc_var->getName() + std::to_string(d), "SCALAR", e_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(r_cc_var->getName(), "VECTOR", r_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(r_cc_var->getName() + std::to_string(d), "SCALAR", r_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(mu_cc_var->getName(), "SCALAR", mu_cc_idx);
        visit_data_writer->registerPlotQuantity(lambda_cc_var->getName(), "SCALAR", lambda_cc_idx);

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

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(r_cc_idx, 0.0);
            level->allocatePatchData(mu_cc_idx, 0.0);
            level->allocatePatchData(lambda_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_cc_var, u_cc_idx, h_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);
        r_vec.addComponent(r_cc_var, r_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);
        muParserCartGridFunction lambda_fcn("lambda", app_initializer->getComponentDatabase("lambda"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
        mu_fcn.setDataOnPatchHierarchy(mu_cc_idx, mu_cc_var, patch_hierarchy, 0.0);
        lambda_fcn.setDataOnPatchHierarchy(lambda_cc_idx, lambda_cc_var, patch_hierarchy, 0.0);

        // Fill ghost cells of shear and bulk viscosity.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comp(2);
        transaction_comp[0] = InterpolationTransactionComponent(mu_cc_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CONSERVATIVE_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                /*mu_bc_coef*/ nullptr,
                                                                Pointer<VariableFillPattern<NDIM> >(nullptr));
        transaction_comp[1] = InterpolationTransactionComponent(lambda_cc_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CONSERVATIVE_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                /*mu_bc_coef*/ nullptr,
                                                                Pointer<VariableFillPattern<NDIM> >(nullptr));

        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
        hier_bdry_fill->setHomogeneousBc(false);
        hier_bdry_fill->fillData(/*time*/ 0.0);

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > U_nul_vecs(NDIM);
        const bool has_velocity_nullspace = periodic_shift.min() > 0;
        if (has_velocity_nullspace)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
            }

            for (unsigned int k = 0; k < NDIM; ++k)
            {
                U_nul_vecs[k] = f_vec.cloneVector("nul_vec_U_" + std::to_string(k));
                U_nul_vecs[k]->allocateVectorData(/*0.0*/);
                U_nul_vecs[k]->setToScalar(0.0);
                for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<CellData<NDIM, double> > U_nul_data =
                            patch->getPatchData(U_nul_vecs[k]->getComponentDescriptorIndex(0));
                        U_nul_data->getArrayData().fill(1.0, k);
                    }
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
            }
        }

        // Setup the implicit viscous and dilatational solver.
        VCViscousDilatationalOpSpec vc_op_spec;
        vc_op_spec.d_C_is_const = true;
        vc_op_spec.d_C_const = 0.0;
        vc_op_spec.d_D_is_const = false;
        vc_op_spec.d_D_idx = mu_cc_idx;
        vc_op_spec.d_L_is_const = false;
        vc_op_spec.d_L_idx = lambda_cc_idx;

        VCCCViscousDilatationalOperator viscous_dil_op("viscous_dil_op");
        viscous_dil_op.setProblemSpecification(&vc_op_spec);
        viscous_dil_op.setPhysicalBcCoefs(u_bc_coefs);
        viscous_dil_op.initializeOperatorState(u_vec, f_vec);
        viscous_dil_op.setSolutionTime(0.0);

        string solver_type = input_db->getString("solver_type");
        Pointer<Database> solver_db = input_db->getDatabase("solver_db");
        string precond_type = input_db->getString("precond_type");
        Pointer<Database> precond_db = input_db->getDatabase("precond_db");

        CCPoissonSolverManager* solver_manager = CCPoissonSolverManager::getManager();
        solver_manager->registerSolverFactoryFunction("VC_CC_VELOCITY_PETSC_KRYLOV_SOLVER",
                                                      allocate_vc_cc_velocity_krylov_solver);
        solver_manager->registerSolverFactoryFunction("VC_CC_VELOCITY_PETSC_LEVEL_SOLVER",
                                                      VCCCViscousDilatationalPETScLevelSolver::allocate_solver);
        Pointer<PoissonSolver> poisson_solver = solver_manager->allocateSolver(solver_type,
                                                                               "vc_cc_velocity_solver",
                                                                               solver_db,
                                                                               "vc_cc_velocity_",
                                                                               precond_type,
                                                                               "vc_cc_velocity_precond",
                                                                               precond_db,
                                                                               "vc_cc_velocity_pc_");

        poisson_solver->setProblemSpecification(&vc_op_spec);
        poisson_solver->setPhysicalBcCoefs(u_bc_coefs);
        poisson_solver->setSolutionTime(0.0);

        // Remove nullspace from RHS
        HierarchyCellDataOpsReal<NDIM, double> cc_data_ops(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        if (has_velocity_nullspace)
        {
            // Ensure that the right-hand-side vector has no components in the
            // nullspace of the operator.
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                r_vec.setToScalar(0.0);
                for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<CellData<NDIM, double> > r_data =
                            patch->getPatchData(r_vec.getComponentDescriptorIndex(0));
                        r_data->getArrayData().fill(1.0, k);
                    }
                }

                double fk_mean = f_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false)) /
                                 r_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));
                cc_data_ops.axpy(f_cc_idx, -fk_mean, r_cc_idx, f_cc_idx);
            }

            Pointer<PETScKrylovLinearSolver> krylov_solver = poisson_solver;
            if (krylov_solver)
            {
                krylov_solver->setNullSpace(false, U_nul_vecs);
            }
        }

        // Solve L*u = f.
        u_vec.setToScalar(0.0);
        poisson_solver->initializeSolverState(u_vec, f_vec);
        poisson_solver->solveSystem(u_vec, f_vec);

        // Subtract mean velocity from velocity components.
        if (has_velocity_nullspace)
        {
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                r_vec.setToScalar(0.0);
                for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<CellData<NDIM, double> > r_data =
                            patch->getPatchData(r_vec.getComponentDescriptorIndex(0));
                        r_data->getArrayData().fill(1.0, k);
                    }
                }

                const double uk_mean = u_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false)) /
                                       r_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));
                cc_data_ops.axpy(u_cc_idx, -uk_mean, r_cc_idx, u_cc_idx);
            }
        }

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&u_vec, false));
        const double e_max_norm = e_vec.maxNorm();
        const double e_l2_norm = e_vec.L2Norm();
        const double e_l1_norm = e_vec.L1Norm();
        pout << "|e|_oo = " << e_max_norm << "\n";
        pout << "|e|_2  = " << e_l2_norm << "\n";
        pout << "|e|_1  = " << e_l1_norm << "\n";

        // Compute the residual and print residual norms.
        viscous_dil_op.setHomogeneousBc(false);
        viscous_dil_op.apply(u_vec, r_vec);
        r_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&f_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&r_vec, false));

        const double r_max_norm = r_vec.maxNorm();
        const double r_l2_norm = r_vec.L2Norm();
        const double r_l1_norm = r_vec.L1Norm();
        pout << "|r|_oo = " << r_max_norm << "\n";
        pout << "|r|_2  = " << r_l2_norm << "\n";
        pout << "|r|_1  = " << r_l1_norm << "\n";

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "|e|_oo = " << e_max_norm << "\n";
            out << "|e|_2  = " << e_l2_norm << "\n";
            out << "|e|_1  = " << e_l1_norm << "\n";
        }

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > e_cc_data = patch->getPatchData(e_cc_idx);
                Pointer<CellData<NDIM, double> > r_cc_data = patch->getPatchData(r_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                        r_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example
