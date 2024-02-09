// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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
#include <ibamr/AcousticStreamingHierarchyIntegrator.h>
#include <ibamr/AcousticStreamingPETScMatUtilities.h>
#include <ibamr/AcousticStreamingPETScVecUtilities.h>
#include <ibamr/FOAcousticStreamingPETScLevelSolver.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <sstream>

#include <ibamr/namespaces.h>

// Routines to reset fluid properties
void
callSetFluidDensityCallbackFunction(int rho_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int /*cycle_num*/,
                                    const double time,
                                    const double /*current_time*/,
                                    const double /*new_time*/,
                                    void* ctx)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    CartGridFunction* rho_fcn = static_cast<CartGridFunction*>(ctx);
    rho_fcn->setDataOnPatchHierarchy(rho_idx, rho_var, patch_hierarchy, time);
    return;
} // callSetFluidDensityCallbackFunction

void
callSetFluidShearViscosityCallbackFunction(int mu_idx,
                                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           const int /*cycle_num*/,
                                           const double time,
                                           const double /*current_time*/,
                                           const double /*new_time*/,
                                           void* ctx)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    CartGridFunction* mu_fcn = static_cast<CartGridFunction*>(ctx);
    mu_fcn->setDataOnPatchHierarchy(mu_idx, mu_var, patch_hierarchy, time);
    return;
} // callSetFluidShearViscosityCallbackFunction

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

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer =
            new AppInitializer(argc, argv, "acoustic_streaming_hier_integrator.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const std::string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create acoustic hierarchy integrator
        Pointer<AcousticStreamingHierarchyIntegrator> time_integrator = new AcousticStreamingHierarchyIntegrator(
            "AcousticStreamingHierarchyIntegrator",
            app_initializer->getComponentDatabase("AcousticStreamingHierarchyIntegrator"),
            /*register_for_restart*/ true);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
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

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u1_init = new muParserCartGridFunction(
            "u1_init", app_initializer->getComponentDatabase("FOVelocityInitialConditions"), grid_geometry);
        time_integrator->registerFirstOrderVelocityInitialConditions(u1_init);
        Pointer<CartGridFunction> p1_init = new muParserCartGridFunction(
            "p1_init", app_initializer->getComponentDatabase("FOPressureInitialConditions"), grid_geometry);
        time_integrator->registerFirstOrderPressureInitialConditions(p1_init);

        Pointer<CartGridFunction> u2_init = new muParserCartGridFunction(
            "u2_init", app_initializer->getComponentDatabase("SOVelocityInitialConditions"), grid_geometry);
        time_integrator->registerSecondOrderVelocityInitialConditions(u2_init);
        Pointer<CartGridFunction> p2_init = new muParserCartGridFunction(
            "p2_init", app_initializer->getComponentDatabase("SOPressureInitialConditions"), grid_geometry);
        time_integrator->registerSecondOrderPressureInitialConditions(p2_init);

        Pointer<CartGridFunction> rho_init =
            new muParserCartGridFunction("rho_init", app_initializer->getComponentDatabase("rho"), grid_geometry);
        time_integrator->registerMassDensityInitialConditions(rho_init);
        Pointer<CartGridFunction> mu_init =
            new muParserCartGridFunction("mu_init", app_initializer->getComponentDatabase("mu"), grid_geometry);
        time_integrator->registerShearViscosityInitialConditions(mu_init);

        // Create boundary condition specification objects for the first-order system.
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2> u1_bc_coefs;
        u1_bc_coefs[0].resize(NDIM);
        u1_bc_coefs[1].resize(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (int comp = 0; comp < 2; ++comp)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    u1_bc_coefs[comp][d] = nullptr;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u1_real_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "FOVelocityRealBcCoefs_" + std::to_string(d);

                u1_bc_coefs[0][d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u1_imag_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "FOVelocityImagBcCoefs_" + std::to_string(d);

                u1_bc_coefs[1][d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }
        time_integrator->registerFirstOrderPhysicalBoundaryConditions(u1_bc_coefs);

        // Create boundary condition specification objects for the second-order system.
        std::vector<RobinBcCoefStrategy<NDIM>*> u2_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u2_bc_coefs[d] = nullptr;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u2_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "SOVelocityBcCoefs_" + std::to_string(d);

                u2_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }
        time_integrator->registerSecondOrderPhysicalBoundaryConditions(u2_bc_coefs);

        // Boundary conditions for density and viscosity
        std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs(NDIM, nullptr);
        if (!(periodic_shift.min() > 0))
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                rho_bc_coefs[d] = new muParserRobinBcCoefs(
                    "rho_bc_coef",
                    app_initializer->getComponentDatabase("DensityBcCoefs_" + std::to_string(d)),
                    grid_geometry);
            }
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ShearViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ShearViscosityBcCoefs"), grid_geometry);
        }

        // Setup the integrator maintained material properties.
        Pointer<SideVariable<NDIM, double> > rho_var = new SideVariable<NDIM, double>("rho_var");
        time_integrator->registerMassDensityVariable(rho_var);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerShearViscosityVariable(mu_var);

        // Callback fncs to reset fluid material properties
        time_integrator->registerResetFluidDensityFcn(&callSetFluidDensityCallbackFunction,
                                                      static_cast<void*>(rho_init.getPointer()));
        time_integrator->registerResetFluidShearViscosityFcn(&callSetFluidShearViscosityCallbackFunction,
                                                             static_cast<void*>(mu_init.getPointer()));

        // Register body force and velocity diveregnce functions for the first and second order systems
        Pointer<CartGridFunction> F1_fcn =
            new muParserCartGridFunction("F1_fcn", app_initializer->getComponentDatabase("F1_fcn"), grid_geometry);
        Pointer<CartGridFunction> F2_fcn =
            new muParserCartGridFunction("F2_fcn", app_initializer->getComponentDatabase("F2_fcn"), grid_geometry);
        Pointer<CartGridFunction> Q1_fcn =
            new muParserCartGridFunction("Q1_fcn", app_initializer->getComponentDatabase("Q1_fcn"), grid_geometry);
        Pointer<CartGridFunction> Q2_fcn =
            new muParserCartGridFunction("Q2_fcn", app_initializer->getComponentDatabase("Q2_fcn"), grid_geometry);

        time_integrator->registerFirstOrderBodyForceFunction(F1_fcn);
        time_integrator->registerSecondOrderBodyForceFunction(F2_fcn);
        time_integrator->registerFirstOrderVelocityDivergenceFunction(Q1_fcn);
        time_integrator->registerSecondOrderVelocityDivergenceFunction(Q2_fcn);

        // Get exact solutions to compute errors
        muParserCartGridFunction u1_exact(
            "u1_exact", app_initializer->getComponentDatabase("FOVelocityExact"), grid_geometry);
        muParserCartGridFunction p1_exact(
            "p1_exact", app_initializer->getComponentDatabase("FOPressureExact"), grid_geometry);
        muParserCartGridFunction u2_exact(
            "u2_exact", app_initializer->getComponentDatabase("SOVelocityExact"), grid_geometry);
        muParserCartGridFunction p2_exact(
            "p2_exact", app_initializer->getComponentDatabase("SOPressureExact"), grid_geometry);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Remove the AppInitializer
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
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<Variable<NDIM> > u1_var = time_integrator->getFirstOrderVelocityVariable();
        const Pointer<Variable<NDIM> > p1_var = time_integrator->getFirstOrderPressureVariable();
        const Pointer<Variable<NDIM> > u2_var = time_integrator->getSecondOrderVelocityVariable();
        const Pointer<Variable<NDIM> > p2_var = time_integrator->getSecondOrderPressureVariable();
        const Pointer<VariableContext> ctx = time_integrator->getCurrentContext();

        const int u1_idx = var_db->mapVariableAndContextToIndex(u1_var, ctx);
        const int u1_cloned_idx = var_db->registerClonedPatchDataIndex(u1_var, u1_idx);
        const int p1_idx = var_db->mapVariableAndContextToIndex(p1_var, ctx);
        const int p1_cloned_idx = var_db->registerClonedPatchDataIndex(p1_var, p1_idx);

        const int u2_idx = var_db->mapVariableAndContextToIndex(u2_var, ctx);
        const int u2_cloned_idx = var_db->registerClonedPatchDataIndex(u2_var, u2_idx);
        const int p2_idx = var_db->mapVariableAndContextToIndex(p2_var, ctx);
        const int p2_cloned_idx = var_db->registerClonedPatchDataIndex(p2_var, p2_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u1_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p1_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u2_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p2_cloned_idx, loop_time);
        }

        u1_exact.setDataOnPatchHierarchy(u1_cloned_idx, u1_var, patch_hierarchy, loop_time);
        p1_exact.setDataOnPatchHierarchy(p1_cloned_idx, p1_var, patch_hierarchy, loop_time);
        u2_exact.setDataOnPatchHierarchy(u2_cloned_idx, u2_var, patch_hierarchy, loop_time);
        p2_exact.setDataOnPatchHierarchy(p2_cloned_idx, p2_var, patch_hierarchy, loop_time);

        Pointer<HierarchyMathOps> hier_math_ops = time_integrator->getHierarchyMathOps();
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();

        std::vector<double> u1_err(3), p1_err(3), u2_err(3), p2_err(3);

        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_sc_data_ops.subtract(u1_idx, u1_idx, u1_cloned_idx);
        u1_err[0] = hier_sc_data_ops.L1Norm(u1_idx, wgt_sc_idx);
        u1_err[1] = hier_sc_data_ops.L2Norm(u1_idx, wgt_sc_idx);
        u1_err[2] = hier_sc_data_ops.maxNorm(u1_idx, wgt_sc_idx);
        pout << "Error in u1 at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << u1_err[0] << "\n"
             << "  L2-norm:  " << u1_err[1] << "\n"
             << "  max-norm: " << u1_err[2] << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        hier_sc_data_ops.subtract(u2_idx, u2_idx, u2_cloned_idx);
        u2_err[0] = hier_sc_data_ops.L1Norm(u2_idx, wgt_sc_idx);
        u2_err[1] = hier_sc_data_ops.L2Norm(u2_idx, wgt_sc_idx);
        u2_err[2] = hier_sc_data_ops.maxNorm(u2_idx, wgt_sc_idx);
        pout << "Error in u2 at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << u2_err[0] << "\n"
             << "  L2-norm:  " << u2_err[1] << "\n"
             << "  max-norm: " << u2_err[2] << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        hier_cc_data_ops.subtract(p1_idx, p1_idx, p1_cloned_idx);
        p1_err[0] = hier_cc_data_ops.L1Norm(p1_idx, wgt_cc_idx);
        p1_err[1] = hier_cc_data_ops.L2Norm(p1_idx, wgt_cc_idx);
        p1_err[2] = hier_cc_data_ops.maxNorm(p1_idx, wgt_cc_idx);

        pout << "Error in p1 at time " << loop_time << ":\n"
             << "  L1-norm:  " << p1_err[0] << "\n"
             << "  L2-norm:  " << p1_err[1] << "\n"
             << "  max-norm: " << p1_err[2] << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        hier_cc_data_ops.subtract(p2_idx, p2_idx, p2_cloned_idx);
        p2_err[0] = hier_cc_data_ops.L1Norm(p2_idx, wgt_cc_idx);
        p2_err[1] = hier_cc_data_ops.L2Norm(p2_idx, wgt_cc_idx);
        p2_err[2] = hier_cc_data_ops.maxNorm(p2_idx, wgt_cc_idx);

        pout << "Error in p2 at time " << loop_time << ":\n"
             << "  L1-norm:  " << p2_err[0] << "\n"
             << "  L2-norm:  " << p2_err[1] << "\n"
             << "  max-norm: " << p2_err[2] << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data && uses_visit)
        {
            pout << "Printing velocity and pressure differences as last visit output" << std::endl;
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

        // Cleanup dumb pointers
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u1_bc_coefs[0][d];
            delete u1_bc_coefs[1][d];

            delete u2_bc_coefs[d];

            delete rho_bc_coefs[d];
        }

        delete mu_bc_coef;
    } // cleanup dynamically allocated objects prior to shutdown
} // main
