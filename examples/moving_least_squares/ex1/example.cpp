// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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
#include <IBAMR_config.h>
#include <IBTK_config.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/KrylovMobilitySolver.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

//////////////////////////////////////////////////////////////////////////////

// Center of mass velocity
void
ConstrainedCOMVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    U_com[1] = 1.0;
    W_com[2] = 1.0;

    return;
} // ConstrainedCOMVel

void
ConstrainedNodalVel(Vec /*U_k*/, const RigidDOFVector& /*U*/, const Eigen::Vector3d& /*X_com*/, void* /*ctx*/)
{
    // intentionally left blank.
    return;
} // ConstrainedNodalVel

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    SAMRAIManager::setMaxNumberPatchDataEntries(2054);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "CIB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Read default Petsc options
        if (input_db->keyExists("petsc_options_file"))
        {
            std::string petsc_options_file = input_db->getString("petsc_options_file");
            PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, petsc_options_file.c_str(), PETSC_TRUE);
        }

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.

        // INS integrator
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        // CIB method
        const unsigned int num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        Pointer<CIBMethod> ib_method_ops =
            new CIBMethod("CIBMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);

        // Krylov solver for INS integrator that solves for [u,p,U,L]
        Pointer<CIBStaggeredStokesSolver> CIBSolver =
            new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver",
                                         input_db->getDatabase("CIBStaggeredStokesSolver"),
                                         navier_stokes_integrator,
                                         ib_method_ops,
                                         "SP_");

        // Register the Krylov solver with INS integrator
        navier_stokes_integrator->setStokesSolver(CIBSolver);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Specify structure kinematics
        FreeRigidDOFVector solve_dofs;
        solve_dofs.setZero();
        ib_method_ops->setSolveRigidBodyVelocity(0, solve_dofs);
        ib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedCOMVel, NULL, 0);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        Pointer<CartGridFunction> sc_mask = new muParserCartGridFunction(
            "sc_mask", app_initializer->getComponentDatabase("VelocityMask"), grid_geometry);
        Pointer<CartGridFunction> cc_mask = new muParserCartGridFunction(
            "cc_mask", app_initializer->getComponentDatabase("PressureMask"), grid_geometry);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
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
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Set physical boundary operator used in spreading.
        ib_method_ops->setVelocityPhysBdryOp(time_integrator->getVelocityPhysBdryOp());

        // Register mobility matrices (if needed)
        std::string mobility_solver_type = input_db->getString("MOBILITY_SOLVER_TYPE");
        if (mobility_solver_type == "DIRECT")
        {
            std::string mat_name = "MOB-MATRIX-1";
            const int managing_proc = 0;
            std::vector<std::vector<unsigned> > struct_ids;
            std::vector<unsigned> prototype_structs;
            prototype_structs.push_back(0);
            struct_ids.push_back(prototype_structs);

            DirectMobilitySolver* direct_solvers = NULL;
            CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(NULL, &direct_solvers, NULL);

            direct_solvers->registerMobilityMat(
                mat_name, prototype_structs, RPY, std::make_pair(LAPACK_LU, LAPACK_LU), managing_proc);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name, struct_ids);
        }
        navier_stokes_integrator->setStokesSolverNeedsInit();

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Setup data used to determine the accuracy of the computed solution.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> main_ctx = var_db->getContext("main_context");

        const std::string weighting_fcn = input_db->getString("DELTA_FUNCTION");
        const int n_ghosts = LEInteractor::getMinimumGhostWidth(weighting_fcn);

        const Pointer<Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const int sc_masked_idx = var_db->registerVariableAndContext(u_var, main_ctx, IntVector<NDIM>(n_ghosts));

        const Pointer<Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        const int cc_masked_idx = var_db->registerVariableAndContext(p_var, main_ctx, IntVector<NDIM>(n_ghosts));
        visit_data_writer->registerPlotQuantity("cc_mask", "SCALAR", cc_masked_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(sc_masked_idx);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(cc_masked_idx);
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        sc_mask->setDataOnPatchHierarchy(sc_masked_idx, u_var, patch_hierarchy, loop_time);
        cc_mask->setDataOnPatchHierarchy(cc_masked_idx, p_var, patch_hierarchy, loop_time);

        // Fill ghost data.
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_components(2);
        ghost_cell_components[0] = ITC(sc_masked_idx,
                                       "CONSERVATIVE_LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // u_bc_coefs
                                       nullptr);
        ghost_cell_components[1] = ITC(cc_masked_idx,
                                       "LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // cc_bc_coefs
                                       nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        // Register sc mask patch data index with CIBMethod
        ib_method_ops->registerMaskingPatchDataIndex(sc_masked_idx);

        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }
        if (dump_postproc_data)
        {
            output_data(patch_hierarchy,
                        ib_method_ops->getLDataManager(),
                        iteration_num,
                        loop_time,
                        postproc_data_dump_dirname);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();

            pout << "Advancing hierarchy by timestep size dt = " << dt << "\n";

            if (time_integrator->atRegridPoint()) navier_stokes_integrator->setStokesSolverNeedsInit();
            if (ib_method_ops->flagRegrid())
            {
                time_integrator->regridHierarchy();
                navier_stokes_integrator->setStokesSolverNeedsInit();
            }
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            RigidDOFVector U0;
            ib_method_ops->getNewRigidBodyVelocity(0, U0);
            pout << "\nRigid body velocity of the structure is : \n" << U0 << "\n";

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
            LDataManager* /*l_data_manager*/,
            const int iteration_num,
            const double loop_time,
            const string& /*data_dump_dirname*/)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    return;
} // output_data
