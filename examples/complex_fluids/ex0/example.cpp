// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/CFINSForcing.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> ins_integrator,
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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

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

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> time_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
        if (solver_type == "STAGGERED")
        {
            time_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            time_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }

        // Create the advection diffusion integrator for the extra stress.
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
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
        time_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);

        // Create boundary condition specification objects (when necessary).
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

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Generate the extra stress component and set up necessary components. The constructor registers the extra
        // stress with the advection diffusion integrator.
        Pointer<CFINSForcing> polymericStressForcing;
        bool using_exact_u = input_db->getBool("USING_EXACT_U");
        if (input_db->keyExists("ComplexFluid"))
        {
            if (!using_exact_u)
            {
                // If we are not using the exact velocity, we need to register the CFINSForcing object with the fluid
                // solver.
                polymericStressForcing = new CFINSForcing("PolymericStressForcing",
                                                          app_initializer->getComponentDatabase("ComplexFluid"),
                                                          time_integrator,
                                                          grid_geometry,
                                                          adv_diff_integrator,
                                                          visit_data_writer);
                time_integrator->registerBodyForceFunction(polymericStressForcing);
            }
            else
            {
                polymericStressForcing = new CFINSForcing("PolymericStressForcing",
                                                          app_initializer->getComponentDatabase("ComplexFluid"),
                                                          u_init,
                                                          grid_geometry,
                                                          adv_diff_integrator,
                                                          visit_data_writer);
            }
        }
        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Create a CartGridFunction that can evaluate the exact solution for the extra stress.
        Pointer<muParserCartGridFunction> s_init = new muParserCartGridFunction(
            "S_exact",
            app_initializer->getComponentDatabase("ComplexFluid")->getDatabase("InitialConditions"),
            grid_geometry);

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
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
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
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy, time_integrator, iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }

        // Compute errors at the final time.
        Pointer<CellVariable<NDIM, double> > s_var = polymericStressForcing->getVariable();
        Pointer<CellVariable<NDIM, double> > sxx_var = new CellVariable<NDIM, double>("Sxx");
        Pointer<CellVariable<NDIM, double> > syy_var = new CellVariable<NDIM, double>("Syy");
        Pointer<CellVariable<NDIM, double> > sxy_var = new CellVariable<NDIM, double>("Sxy");
        const Pointer<VariableContext> s_ctx = adv_diff_integrator->getCurrentContext();

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int s_idx = var_db->mapVariableAndContextToIndex(s_var, s_ctx);
        const int s_cloned_idx = var_db->registerClonedPatchDataIndex(s_var, s_idx);
        const int sxx_idx = var_db->registerVariableAndContext(sxx_var, s_ctx);
        const int syy_idx = var_db->registerVariableAndContext(syy_var, s_ctx);
        const int sxy_idx = var_db->registerVariableAndContext(sxy_var, s_ctx);

        const Pointer<Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();

        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(s_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(sxx_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(syy_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(sxy_idx, loop_time);
            if (!using_exact_u) patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
        }

        s_init->setDataOnPatchHierarchy(s_cloned_idx, s_var, patch_hierarchy, loop_time);
        if (!using_exact_u) u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        if (!using_exact_u) hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
        hier_cc_data_ops.subtract(s_idx, s_idx, s_cloned_idx);

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > s_data = patch->getPatchData(s_idx);
                Pointer<CellData<NDIM, double> > sxx_data = patch->getPatchData(sxx_idx);
                Pointer<CellData<NDIM, double> > syy_data = patch->getPatchData(syy_idx);
                Pointer<CellData<NDIM, double> > sxy_data = patch->getPatchData(sxy_idx);

                sxx_data->copyDepth(0, *s_data, 0);
                syy_data->copyDepth(0, *s_data, 1);
                sxy_data->copyDepth(0, *s_data, 2);
            }
        }
        if (!using_exact_u)
        {
            pout << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << std::setprecision(10) << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << std::setprecision(10) << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n";
        }

        pout << "Error in sxx at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(sxx_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(sxx_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(sxx_idx, wgt_cc_idx) << "\n";

        pout << "Error in syy at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(syy_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(syy_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(syy_idx, wgt_cc_idx) << "\n";

        pout << "Error in sxy at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(sxy_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << std::setprecision(10) << hier_cc_data_ops.L2Norm(sxy_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << std::setprecision(10) << hier_cc_data_ops.maxNorm(sxy_idx, wgt_cc_idx) << "\n";

        time_integrator->setupPlotData();
        visit_data_writer->registerPlotQuantity("SXX_Err", "SCALAR", sxx_idx);
        visit_data_writer->registerPlotQuantity("SYY_Err", "SCALAR", syy_idx);
        visit_data_writer->registerPlotQuantity("SXY_Err", "SCALAR", sxy_idx);
        visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);

        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (!using_exact_u) level->deallocatePatchData(u_cloned_idx);
            level->deallocatePatchData(s_cloned_idx);
            level->deallocatePatchData(sxx_idx);
            level->deallocatePatchData(syy_idx);
            level->deallocatePatchData(sxy_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> ins_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getVelocityVariable(),
                                                           ins_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getPressureVariable(),
                                                           ins_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data
