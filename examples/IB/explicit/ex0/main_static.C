// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_prefix_config.h>
#include <IBTK_prefix_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

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
main(
    int argc,
    char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

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
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator", app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n" <<
                       "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBHierarchyIntegrator> time_integrator = new IBHierarchyIntegrator(
            "IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        time_integrator->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        time_integrator->registerIBLagrangianForceFunction(ib_force_fcn);

        // Create Eulerian initial condition specification objects.  These
        // objects also are used to specify exact solution values for error
        // analysis.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> u_bc_coefs;
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
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        time_integrator->freeLInitStrategy();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Setup data used to determine the accuracy of the computed solution.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();

        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        const Pointer<Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = navier_stokes_integrator->getCurrentContext();

        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
        visit_data_writer->registerPlotQuantity("P error", "SCALAR", p_cloned_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx);
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Make a copy of the initial positions of the structure.
        LDataManager* l_data_manager = time_integrator->getLDataManager();
        const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
        Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
        Vec X_petsc_vec = X_data->getVec();
        Vec X_lag_vec;
        VecDuplicate(X_petsc_vec, &X_lag_vec);
        l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout <<                                                    "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";

            dt = time_integrator->getTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout <<                                                    "\n";
            pout << "At end       of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout <<                                                    "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num%viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            // Compute velocity and pressure error norms.
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(u_cloned_idx)) level->allocatePatchData(u_cloned_idx);
                if (!level->checkAllocated(p_cloned_idx)) level->allocatePatchData(p_cloned_idx);
            }

            pout << "\n"
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                 << "Computing error norms.\n\n";

            u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
            p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time-0.5*dt);

            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

            Pointer<CellVariable<NDIM,double> > u_cc_var = u_var;
            if (!u_cc_var.isNull())
            {
                HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
                hier_cc_data_ops.subtract(u_cloned_idx, u_idx, u_cloned_idx);
                pout << "Error in u at time " << loop_time << ":\n"
                     << "  L1-norm:  " << hier_cc_data_ops. L1Norm(u_cloned_idx,wgt_cc_idx)  << "\n"
                     << "  L2-norm:  " << hier_cc_data_ops. L2Norm(u_cloned_idx,wgt_cc_idx)  << "\n"
                     << "  max-norm: " << hier_cc_data_ops.maxNorm(u_cloned_idx,wgt_cc_idx) << "\n";
            }

            Pointer<SideVariable<NDIM,double> > u_sc_var = u_var;
            if (!u_sc_var.isNull())
            {
                HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
                hier_sc_data_ops.subtract(u_cloned_idx, u_idx, u_cloned_idx);
                pout << "Error in u at time " << loop_time << ":\n"
                     << "  L1-norm:  " << hier_sc_data_ops. L1Norm(u_cloned_idx,wgt_sc_idx)  << "\n"
                     << "  L2-norm:  " << hier_sc_data_ops. L2Norm(u_cloned_idx,wgt_sc_idx)  << "\n"
                     << "  max-norm: " << hier_sc_data_ops.maxNorm(u_cloned_idx,wgt_sc_idx) << "\n";
            }

            HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            hier_cc_data_ops.subtract(p_cloned_idx, p_idx, p_cloned_idx);
            pout << "Error in p at time " << loop_time-0.5*dt << ":\n"
                 << "  L1-norm:  " << hier_cc_data_ops. L1Norm(p_cloned_idx,wgt_cc_idx)  << "\n"
                 << "  L2-norm:  " << hier_cc_data_ops. L2Norm(p_cloned_idx,wgt_cc_idx)  << "\n"
                 << "  max-norm: " << hier_cc_data_ops.maxNorm(p_cloned_idx,wgt_cc_idx) << "\n"
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        }

        // Estimate the rate of volume change.
        Pointer<LData> U_data = l_data_manager->getLData("U", finest_hier_level);
        Vec U_petsc_vec = U_data->getVec();
        Vec U_lag_vec;
        VecDuplicate(U_petsc_vec, &U_lag_vec);
        l_data_manager->scatterPETScToLagrangian(U_petsc_vec, U_lag_vec, finest_hier_level);
        double* X_arr;
        VecGetArray(X_lag_vec, &X_arr);
        double* U_arr;
        VecGetArray(U_lag_vec, &U_arr);
        const int n_layers = input_db->getInteger("n_layers");
        vector<double> dV(n_layers > 1 ? 2 : 1, 0.0);
        const int n_nodes = X_data->getLocalNodeCount();
        const int n_nodes_per_layer = n_nodes/n_layers;
        for (int side = 0; side <= (n_layers > 1 ? 1 : 0); ++side)
        {
            for (int k = 0; k < n_nodes_per_layer; ++k)
            {
                int k_prev = (k-1)%n_nodes_per_layer;
                int k_next = (k+1)%n_nodes_per_layer;
                double* X = &X_arr[NDIM*k + side*(n_layers-1)*n_nodes_per_layer];
                double* U = &U_arr[NDIM*k + side*(n_layers-1)*n_nodes_per_layer];
                double n[NDIM];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    n[d] = X[d] - 0.5;
                }
                const double n_norm = sqrt(inner_product(n,n+NDIM,n,0.0));
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    n[d] /= n_norm;
                }

                double* X_prev = &X_arr[NDIM*k_prev + side*(n_layers-1)*n_nodes_per_layer];
                double* X_next = &X_arr[NDIM*k_next + side*(n_layers-1)*n_nodes_per_layer];
                double dX[NDIM];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    dX[d] = 0.5*(X_next[d]-X_prev[d]);
                }
                const double dx = sqrt(inner_product(dX,dX+NDIM,dX,0.0));

                dV[side] += inner_product(U,U+NDIM,n,0.0)*dx;
            }
        }
        if (dV.size() == 1) pout << "dV = " << dV[0] << "\n";
        if (dV.size() == 2) pout << "dV_i = " << -dV[0] << "\ndV_o = " << dV[1] << "\n";
        VecRestoreArray(X_lag_vec, &X_arr);
        VecRestoreArray(U_lag_vec, &U_arr);
        VecDestroy(X_lag_vec);
        VecDestroy(U_lag_vec);

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}// main
