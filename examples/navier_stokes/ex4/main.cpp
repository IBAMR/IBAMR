// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
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
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> time_integrator;
        const string ins_solver_type = main_db->getStringWithDefault("ins_solver_type", "STAGGERED");
        if (ins_solver_type == "STAGGERED")
        {
            time_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (ins_solver_type == "COLLOCATED")
        {
            time_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << ins_solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = main_db->getStringWithDefault("adv_diff_solver_type", "GODUNOV");
        if (adv_diff_solver_type == "GODUNOV")
        {
            Pointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps(
                "AdvectorExplicitPredictorPatchOps",
                app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
            adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator(
                "AdvDiffPredictorCorrectorHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"),
                predictor);
        }
        else if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: GODUNOV, SEMI_IMPLICIT");
        }
        time_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
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

        // Set up the fluid solver.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);

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
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();
                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Setup the advected and diffused quantity.
        Pointer<CellVariable<NDIM, double> > U_adv_diff_var = new CellVariable<NDIM, double>("U_adv_diff", NDIM);
        adv_diff_integrator->registerTransportedQuantity(U_adv_diff_var);
        adv_diff_integrator->setDiffusionCoefficient(U_adv_diff_var,
                                                     input_db->getDouble("MU") / input_db->getDouble("RHO"));
        adv_diff_integrator->setInitialConditions(U_adv_diff_var, u_init);

        Pointer<FaceVariable<NDIM, double> > u_adv_var = time_integrator->getAdvectionVelocityVariable();
        adv_diff_integrator->setAdvectionVelocity(U_adv_diff_var, u_adv_var);

        vector<RobinBcCoefStrategy<NDIM>*> U_adv_diff_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_adv_diff_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "U_adv_diff_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();
                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "AdvDiffBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
                U_adv_diff_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            adv_diff_integrator->setPhysicalBcCoefs(U_adv_diff_var, U_adv_diff_bc_coefs);
        }

        if (input_db->keyExists("AdvDiffForcingFunction"))
        {
            Pointer<CellVariable<NDIM, double> > F_adv_diff_var = new CellVariable<NDIM, double>("F_adv_diff_", NDIM);
            Pointer<CartGridFunction> F_adv_diff_fcn = new muParserCartGridFunction(
                "F_adv_diff_fcn", app_initializer->getComponentDatabase("AdvDiffForcingFunction"), grid_geometry);
            adv_diff_integrator->registerSourceTerm(F_adv_diff_var);
            adv_diff_integrator->setSourceTermFunction(F_adv_diff_var, F_adv_diff_fcn);
            adv_diff_integrator->setSourceTerm(U_adv_diff_var, F_adv_diff_var);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
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
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

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
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
        const Pointer<Variable<NDIM> > p_var = time_integrator->getPressureVariable();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, time_integrator->getCurrentContext());
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, time_integrator->getCurrentContext());
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
        const int U_adv_diff_idx =
            var_db->mapVariableAndContextToIndex(U_adv_diff_var, adv_diff_integrator->getCurrentContext());
        const int U_adv_diff_cloned_idx = var_db->registerClonedPatchDataIndex(U_adv_diff_var, U_adv_diff_idx);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(U_adv_diff_cloned_idx, loop_time);
        }
        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);
        u_init->setDataOnPatchHierarchy(U_adv_diff_cloned_idx, U_adv_diff_var, patch_hierarchy, loop_time);

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        Pointer<CellVariable<NDIM, double> > u_cc_var = u_var;
        if (u_cc_var)
        {
            hier_cc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
            pout << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << hier_cc_data_ops.L1Norm(u_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << hier_cc_data_ops.L2Norm(u_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << hier_cc_data_ops.maxNorm(u_idx, wgt_cc_idx) << "\n";
        }

        Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
        if (u_sc_var)
        {
            hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
            pout << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n";
        }

        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        pout << "Error in p at time " << loop_time - 0.5 * dt << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n";

        hier_cc_data_ops.subtract(U_adv_diff_idx, U_adv_diff_idx, U_adv_diff_cloned_idx);
        pout << "Error in U_adv_diff at time " << loop_time << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(U_adv_diff_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(U_adv_diff_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(U_adv_diff_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data && uses_visit)
        {
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        for (unsigned int d = 0; d < NDIM; ++d) delete U_adv_diff_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main
