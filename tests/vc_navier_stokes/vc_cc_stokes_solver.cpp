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
#include <ibamr/VCCollocatedStokesOperator.h>
#include <ibamr/VCCollocatedStokesProjectionPreconditioner.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScKrylovPoissonSolver.h>
#include <ibtk/ProblemSpecification.h>
#include <ibtk/VCCCViscousDilatationalOpPointRelaxationFACOperator.h>
#include <ibtk/VCCCViscousDilatationalOperator.h>
#include <ibtk/VCCCViscousDilatationalPETScLevelSolver.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>
#include <ibtk/app_namespaces.h>

Pointer<PETScKrylovLinearSolver>
allocate_vc_collocated_stokes_krylov_solver(const std::string& solver_name,
                                            const std::string& solver_options_prefix,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db)
{
    Pointer<PETScKrylovLinearSolver> stokes_solver =
        new PETScKrylovLinearSolver(solver_name, solver_input_db, solver_options_prefix);
    stokes_solver->setOperator(new VCCollocatedStokesOperator(
        solver_name + "::vc_cc_collocated_stokes_operator", /*homogeneous_bc*/ true, solver_input_db));
    stokes_solver->setPreconditioner(new VCCollocatedStokesProjectionPreconditioner(
        solver_name + "::vc_collocated_stokes_precond", solver_input_db));
    return stokes_solver;
}

Pointer<PoissonSolver>
allocate_vc_cc_velocity_krylov_solver(const std::string& solver_name,
                                      const std::string& solver_options_prefix,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> pc_input_db)
{
    Pointer<PETScKrylovPoissonSolver> velocity_solver =
        new PETScKrylovPoissonSolver(solver_name, solver_input_db, solver_options_prefix);

    CCPoissonSolverManager* solver_manager = CCPoissonSolverManager::getManager();
    solver_manager->registerSolverFactoryFunction("VC_CC_VELOCITY_PETSC_LEVEL_SOLVER",
                                                  VCCCViscousDilatationalPETScLevelSolver::allocate_solver);
    Pointer<PoissonSolver> velocity_precond = VCCCViscousDilatationalOpPointRelaxationFACOperator::allocate_solver(
        solver_name + "::vc_cc_velocity_precond", pc_input_db, "vc_cc_velocity_pc_");

    velocity_solver->setOperator(new VCCCViscousDilatationalOperator(solver_name + "::vc_cc_viscous_operator"));
    velocity_solver->setPreconditioner(velocity_precond);

    return velocity_solver;
}

Pointer<PoissonSolver>
allocate_pressure_krylov_solver(const std::string& solver_name,
                                const std::string& solver_options_prefix,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> pc_input_db)
{
    Pointer<PETScKrylovPoissonSolver> pressure_solver =
        new PETScKrylovPoissonSolver(solver_name, solver_input_db, solver_options_prefix);
    Pointer<PoissonSolver> pressure_precond = CCPoissonPointRelaxationFACOperator::allocate_solver(
        solver_name + "::pressure_precond", pc_input_db, "vc_cc_pressure_pc_");

    pressure_solver->setOperator(new CCLaplaceOperator(solver_name + "::cc_laplace_operator"));
    pressure_solver->setPreconditioner(pressure_precond);

    return pressure_solver;
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

        Pointer<CellVariable<NDIM, double> > u_var = new CellVariable<NDIM, double>("u", NDIM);
        Pointer<CellVariable<NDIM, double> > fu_var = new CellVariable<NDIM, double>("fu", NDIM);
        Pointer<CellVariable<NDIM, double> > p_var = new CellVariable<NDIM, double>("p");
        Pointer<CellVariable<NDIM, double> > fp_var = new CellVariable<NDIM, double>("fp");
        Pointer<CellVariable<NDIM, double> > eu_var = new CellVariable<NDIM, double>("eu", NDIM);
        Pointer<CellVariable<NDIM, double> > ru_var = new CellVariable<NDIM, double>("ru", NDIM);
        Pointer<CellVariable<NDIM, double> > ep_var = new CellVariable<NDIM, double>("ep");
        Pointer<CellVariable<NDIM, double> > rp_var = new CellVariable<NDIM, double>("rp");

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        const int mu_idx = var_db->registerVariableAndContext(mu_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > lambda_var = new CellVariable<NDIM, double>("lambda");
        const int lambda_idx = var_db->registerVariableAndContext(lambda_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > C_var = new CellVariable<NDIM, double>("C", NDIM);
        const int C_idx = var_db->registerVariableAndContext(C_var, ctx, IntVector<NDIM>(1));

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
        const int fu_idx = var_db->registerVariableAndContext(fu_var, ctx, IntVector<NDIM>(1));
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
        const int fp_idx = var_db->registerVariableAndContext(fp_var, ctx, IntVector<NDIM>(1));
        const int eu_idx = var_db->registerVariableAndContext(eu_var, ctx, IntVector<NDIM>(1));
        const int ep_idx = var_db->registerVariableAndContext(ep_var, ctx, IntVector<NDIM>(1));
        const int ru_idx = var_db->registerVariableAndContext(ru_var, ctx, IntVector<NDIM>(1));
        const int rp_idx = var_db->registerVariableAndContext(rp_var, ctx, IntVector<NDIM>(1));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_var->getName(), "VECTOR", u_idx);
        visit_data_writer->registerPlotQuantity(fu_var->getName(), "VECTOR", fu_idx);
        visit_data_writer->registerPlotQuantity(eu_var->getName(), "VECTOR", eu_idx);
        visit_data_writer->registerPlotQuantity(ru_var->getName(), "VECTOR", ru_idx);
        visit_data_writer->registerPlotQuantity(C_var->getName(), "VECTOR", C_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_var->getName() + std::to_string(d), "SCALAR", u_idx, d);
            visit_data_writer->registerPlotQuantity(fu_var->getName() + std::to_string(d), "SCALAR", fu_idx, d);
            visit_data_writer->registerPlotQuantity(eu_var->getName() + std::to_string(d), "SCALAR", eu_idx, d);
            visit_data_writer->registerPlotQuantity(ru_var->getName() + std::to_string(d), "SCALAR", ru_idx, d);
            visit_data_writer->registerPlotQuantity(C_var->getName() + std::to_string(d), "SCALAR", C_idx, d);
        }

        visit_data_writer->registerPlotQuantity(p_var->getName(), "SCALAR", p_idx);
        visit_data_writer->registerPlotQuantity(fp_var->getName(), "SCALAR", fp_idx);
        visit_data_writer->registerPlotQuantity(ep_var->getName(), "SCALAR", ep_idx);
        visit_data_writer->registerPlotQuantity(rp_var->getName(), "SCALAR", rp_idx);
        visit_data_writer->registerPlotQuantity(mu_var->getName(), "SCALAR", mu_idx);
        visit_data_writer->registerPlotQuantity(lambda_var->getName(), "SCALAR", lambda_idx);

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

            level->allocatePatchData(u_idx, 0.0);
            level->allocatePatchData(fu_idx, 0.0);
            level->allocatePatchData(eu_idx, 0.0);
            level->allocatePatchData(ru_idx, 0.0);
            level->allocatePatchData(p_idx, 0.0);
            level->allocatePatchData(fp_idx, 0.0);
            level->allocatePatchData(ep_idx, 0.0);
            level->allocatePatchData(rp_idx, 0.0);
            level->allocatePatchData(mu_idx, 0.0);
            level->allocatePatchData(lambda_idx, 0.0);
            level->allocatePatchData(C_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(
            patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        const int h_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> p_vec("p", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> fu_vec("fu", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> fp_vec("fp", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> ru_vec("ru", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> rp_vec("rp", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        SAMRAIVectorReal<NDIM, double> x_vec("up", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        // Individual vecs
        u_vec.addComponent(u_var, u_idx, h_idx);
        p_vec.addComponent(p_var, p_idx, h_idx);
        fu_vec.addComponent(fu_var, fu_idx, h_idx);
        fp_vec.addComponent(fp_var, fp_idx, h_idx);
        ru_vec.addComponent(ru_var, ru_idx, h_idx);
        rp_vec.addComponent(rp_var, rp_idx, h_idx);

        // Composite vecs
        x_vec.addComponent(u_var, u_idx, h_idx);
        x_vec.addComponent(p_var, p_idx, h_idx);
        f_vec.addComponent(fu_var, fu_idx, h_idx);
        f_vec.addComponent(fp_var, fp_idx, h_idx);
        e_vec.addComponent(eu_var, eu_idx, h_idx);
        e_vec.addComponent(ep_var, ep_idx, h_idx);
        r_vec.addComponent(ru_var, ru_idx, h_idx);
        r_vec.addComponent(rp_var, rp_idx, h_idx);

        // Initialize composite vecs. This will also initialize individual vecs.
        x_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction fu_fcn("fu", app_initializer->getComponentDatabase("fu"), grid_geometry);
        muParserCartGridFunction p_fcn("p", app_initializer->getComponentDatabase("p"), grid_geometry);
        muParserCartGridFunction fp_fcn("fp", app_initializer->getComponentDatabase("fp"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);
        muParserCartGridFunction lambda_fcn("lambda", app_initializer->getComponentDatabase("lambda"), grid_geometry);
        muParserCartGridFunction C_fcn("C", app_initializer->getComponentDatabase("C"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(eu_idx, eu_var, patch_hierarchy, 0.0);
        fu_fcn.setDataOnPatchHierarchy(fu_idx, fu_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(ep_idx, ep_var, patch_hierarchy, 0.0);
        fp_fcn.setDataOnPatchHierarchy(fp_idx, fp_var, patch_hierarchy, 0.0);
        mu_fcn.setDataOnPatchHierarchy(mu_idx, mu_var, patch_hierarchy, 0.0);
        lambda_fcn.setDataOnPatchHierarchy(lambda_idx, lambda_var, patch_hierarchy, 0.0);
        C_fcn.setDataOnPatchHierarchy(C_idx, C_var, patch_hierarchy, 0.0);

        // Make mu and lambda negative
        hier_cc_data_ops.scale(mu_idx, -1.0, mu_idx);
        hier_cc_data_ops.scale(lambda_idx, -1.0, lambda_idx);

        // Fill ghost cells of shear and bulk viscosity.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comp(3);
        transaction_comp[0] = InterpolationTransactionComponent(mu_idx,
                                                                /*DATA_REFINE_TYPE*/ "LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                /*mu_bc_coef*/ nullptr,
                                                                Pointer<VariableFillPattern<NDIM> >(nullptr));
        transaction_comp[1] = InterpolationTransactionComponent(lambda_idx,
                                                                /*DATA_REFINE_TYPE*/ "LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                /*mu_bc_coef*/ nullptr,
                                                                Pointer<VariableFillPattern<NDIM> >(nullptr));
        transaction_comp[2] = InterpolationTransactionComponent(C_idx,
                                                                /*DATA_REFINE_TYPE*/ "LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
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
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

            const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

            Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
            u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
        }
        RobinBcCoefStrategy<NDIM>* P_bc_coef = nullptr;
        {
            const std::string bc_coefs_name = "p_bc_coefs";
            const std::string bc_coefs_db_name = "PressureBcCoefs";
            Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
            P_bc_coef = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
        }

        // Setup the implicit viscous and dilatational solver.
        VCViscousDilatationalOpSpec vc_op_spec;
        vc_op_spec.d_C_is_const = false;
        vc_op_spec.d_C_idx = C_idx;
        vc_op_spec.d_D_is_const = false;
        vc_op_spec.d_D_idx = mu_idx;
        vc_op_spec.d_L_is_const = false;
        vc_op_spec.d_L_idx = lambda_idx;

        string stokes_solver_name = "vc_cc_stokes_solver";
        string stokes_solver_prefix = "vc_cc_stokes_";
        string velocity_solver_name = "vc_cc_velocity_solver";
        string velocity_solver_prefix = "vc_cc_velocity_";
        string pressure_solver_name = "pressure_solver";
        string pressure_solver_prefix = "pressure_";
        Pointer<Database> stokes_solver_db = input_db->getDatabase("vc_cc_stokes_db");
        Pointer<Database> velocity_solver_db = stokes_solver_db->getDatabase("velocity_solver_db");
        Pointer<Database> velocity_pc_db = stokes_solver_db->getDatabase("velocity_precond_db");
        Pointer<Database> pressure_solver_db = stokes_solver_db->getDatabase("pressure_solver_db");
        Pointer<Database> pressure_pc_db = stokes_solver_db->getDatabase("pressure_precond_db");

        Pointer<PETScKrylovLinearSolver> stokes_solver =
            allocate_vc_collocated_stokes_krylov_solver(stokes_solver_name, stokes_solver_prefix, stokes_solver_db);
        Pointer<PoissonSolver> velocity_solver = allocate_vc_cc_velocity_krylov_solver(
            velocity_solver_name, velocity_solver_prefix, velocity_solver_db, velocity_pc_db);
        velocity_solver->setProblemSpecification(&vc_op_spec);
        velocity_solver->setPhysicalBcCoefs(u_bc_coefs);
        Pointer<PoissonSolver> pressure_solver = allocate_pressure_krylov_solver(
            pressure_solver_name, pressure_solver_prefix, pressure_solver_db, pressure_pc_db);
        PoissonSpecifications pressure_poisson_spec("pressure_solver_spec");
        pressure_poisson_spec.setCZero();
        pressure_poisson_spec.setDConstant(-1.0);
        pressure_solver->setPoissonSpecifications(pressure_poisson_spec);
        pressure_solver->setPhysicalBcCoef(P_bc_coef);

        Pointer<VCCollocatedStokesOperator> stokes_op = stokes_solver->getOperator();
        Pointer<VCCollocatedStokesProjectionPreconditioner> stokes_pc = stokes_solver->getPreconditioner();
        stokes_pc->setVelocitySubdomainSolver(velocity_solver);
        stokes_pc->setPressureSubdomainSolver(pressure_solver);

        stokes_op->setProblemSpecification(&vc_op_spec);
        stokes_op->setPhysicalBcCoefs(u_bc_coefs, P_bc_coef);
        stokes_pc->setProblemSpecification(&vc_op_spec);
        stokes_pc->setPhysicalBcCoefs(u_bc_coefs, P_bc_coef);

        // Set nullspace vectors
        std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > stokes_null_vec(1);
        auto& x_null_vec = stokes_null_vec[0];
        x_null_vec = x_vec.cloneVector("stokes_null_vec");
        x_null_vec->allocateVectorData(0.0);
        hier_cc_data_ops.setToScalar(x_null_vec->getComponentDescriptorIndex(0), 0.0);
        hier_cc_data_ops.setToScalar(x_null_vec->getComponentDescriptorIndex(1), 1.0);
        stokes_solver->setNullSpace(false, stokes_null_vec);
        auto p_pressure_solver = dynamic_cast<LinearSolver*>(pressure_solver.getPointer());
        p_pressure_solver->setNullSpace(true);

        // Initialize velocity and pressure solvers
        velocity_solver->initializeSolverState(u_vec, fu_vec);
        pressure_solver->initializeSolverState(p_vec, fp_vec);

        // Remove the nullspace of the operator from the RHS vector
        for (unsigned int k = 0; k < NDIM; ++k)
        {
            ru_vec.setToScalar(0.0);
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > ru_data =
                        patch->getPatchData(ru_vec.getComponentDescriptorIndex(0));
                    ru_data->getArrayData().fill(1.0, k);
                }
            }

            double fk_mean = fu_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&ru_vec, false)) /
                             ru_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&ru_vec, false));
            hier_cc_data_ops.axpy(fu_idx, -fk_mean, ru_idx, fu_idx);
        }
        {
            rp_vec.setToScalar(1.0);
            double fp_mean = fp_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&rp_vec, false)) /
                             rp_vec.dot(Pointer<SAMRAIVectorReal<NDIM, double> >(&rp_vec, false));
            hier_cc_data_ops.axpy(fp_idx, -fp_mean, rp_idx, fp_idx);
        }

        // Solve L*u = f.
        x_vec.setToScalar(0.0);
        stokes_solver->initializeSolverState(x_vec, f_vec);
        stokes_solver->solveSystem(x_vec, f_vec);

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&x_vec, false));
        const double e_max_norm = e_vec.maxNorm();
        const double e_l2_norm = e_vec.L2Norm();
        const double e_l1_norm = e_vec.L1Norm();
        pout << "|e|_oo = " << e_max_norm << "\n";
        pout << "|e|_2  = " << e_l2_norm << "\n";
        pout << "|e|_1  = " << e_l1_norm << "\n";

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
                Pointer<CellData<NDIM, double> > eu_data = patch->getPatchData(eu_idx);
                Pointer<CellData<NDIM, double> > ep_data = patch->getPatchData(ep_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        eu_data->fillAll(0.0, intersection);
                        ep_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example
