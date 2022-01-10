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

// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/BrinkmanAdvDiffBcHelper.h>
#include <ibamr/BrinkmanAdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBInterpolantHierarchyIntegrator.h>

#include "ibtk/IndexUtilities.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Application specific includes
#include "LevelSetInitialConditionCircle.h"
#include "LevelSetInitialConditionHexagram.h"

// Struct to specify the function required for inhomogeneous Neumann conditions for Brinkman penalization
struct BrinkmanPenalizationCtx
{
    // intentionally blank
};

void
evaluate_brinkman_bc_callback_fcn(int B_idx,
                                  Pointer<CellVariable<NDIM, double> > /*ls_var*/,
                                  Pointer<HierarchyMathOps> hier_math_ops,
                                  double /*time*/,
                                  void* /*ctx*/)
{
#if (!NDEBUG)
    TBOX_ASSERT(NDIM == 2);
#endif

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    IntVector<NDIM> ratio = finest_level->getRatio();
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_hierarchy->getGridGeometry();
        Pointer<SideData<NDIM, double> > B_data = patch->getPatchData(B_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

                // Get physical coordinates
                IBTK::Vector coord = IndexUtilities::getSideCenter(grid_geom, ratio, si);
                if (axis == 0)
                {
                    (*B_data)(si) = cos(coord[0]) * sin(coord[1]);
                }
                else if (axis == 1)
                {
                    (*B_data)(si) = sin(coord[0]) * cos(coord[1]);
                }
            }
        }
    }
}

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

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
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
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        //
        Pointer<AdvDiffHierarchyIntegrator> time_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "PREDICTOR_CORRECTOR");
        time_integrator = new BrinkmanAdvDiffSemiImplicitHierarchyIntegrator(
            "BrinkmanAdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("BrinkmanAdvDiffSemiImplicitHierarchyIntegrator"));

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

        // Set up the advected and diffused quantity.
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();

        Pointer<CellVariable<NDIM, double> > phi_outer_solid_var =
            new CellVariable<NDIM, double>("level_set_outer_solid");
        time_integrator->registerTransportedQuantity(phi_outer_solid_var, true);
        time_integrator->setDiffusionCoefficient(phi_outer_solid_var, 0.0);

        // Origin of the Hexagram.
        IBTK::Vector2d origin1(M_PI, M_PI);
        Pointer<CartGridFunction> phi_outer_solid_init =
            new LevelSetInitialConditionHexagram("ls_outer_solid_init", grid_geometry, origin1);
        time_integrator->setInitialConditions(phi_outer_solid_var, phi_outer_solid_init);
        std::vector<Pointer<CellVariable<NDIM, double> > > ls_vars;
        ls_vars.push_back(phi_outer_solid_var);

        const string& ls_name = "level_set_inner_solid";
        Pointer<CellVariable<NDIM, double> > phi_inner_solid_var = new CellVariable<NDIM, double>(ls_name);
        time_integrator->registerTransportedQuantity(phi_inner_solid_var, true);
        time_integrator->setDiffusionCoefficient(phi_inner_solid_var, 0.0);

        // Origin of the circle1.
        IBTK::Vector2d origin2(M_PI, M_PI);
        const double radius = input_db->getDouble("RADIUS");
        Pointer<CartGridFunction> phi_inner_solid_init =
            new LevelSetInitialConditionCircle("ls_inner_solid_init", grid_geometry, radius, origin2);
        time_integrator->setInitialConditions(phi_inner_solid_var, phi_inner_solid_init);
        ls_vars.push_back(phi_inner_solid_var);

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
            time_integrator->setPhysicalBcCoef(phi_outer_solid_var, phi_bc_coef);
            time_integrator->setPhysicalBcCoef(phi_inner_solid_var, phi_bc_coef);
        }
        Pointer<CellVariable<NDIM, double> > q_var = new CellVariable<NDIM, double>("q");
        time_integrator->registerTransportedQuantity(q_var, true);
        time_integrator->setDiffusionCoefficient(q_var, input_db->getDouble("KAPPA"));

        if (input_db->keyExists("TransportedQuantityInitialConditions"))
        {
            Pointer<CartGridFunction> q_init = new muParserCartGridFunction(
                "q_init", app_initializer->getComponentDatabase("TransportedQuantityInitialConditions"), grid_geometry);
            time_integrator->setInitialConditions(q_var, q_init);
        }

        RobinBcCoefStrategy<NDIM>* q_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TransportedQuantityBcCoefs"))
        {
            q_bc_coef = new muParserRobinBcCoefs(
                "q_bc_coef", app_initializer->getComponentDatabase("TransportedQuantityBcCoefs"), grid_geometry);
        }
        time_integrator->setPhysicalBcCoef(q_var, q_bc_coef);

        // Brinkman penalization for Trsanported quantity.
        const string indicator_function_type = input_db->getString("INDICATOR_FUNCTION_TYPE");
        const double eta = input_db->getDouble("ETA");
        const double num_of_interface_cells = input_db->getDouble("NUMBER_OF_INTERFACE_CELLS");
        Pointer<BrinkmanAdvDiffBcHelper> brinkman_adv_diff =
            new BrinkmanAdvDiffBcHelper("BrinkmanAdvDiffBcHelper", time_integrator);

        // setting inhomogeneous Neumann at the lower cylinder surface.
        brinkman_adv_diff->registerInhomogeneousBC(q_var,
                                                   phi_outer_solid_var,
                                                   "NEUMANN",
                                                   &evaluate_brinkman_bc_callback_fcn,
                                                   nullptr,
                                                   indicator_function_type,
                                                   num_of_interface_cells,
                                                   eta);

        brinkman_adv_diff->registerInhomogeneousBC(q_var,
                                                   phi_inner_solid_var,
                                                   "NEUMANN",
                                                   &evaluate_brinkman_bc_callback_fcn,
                                                   nullptr,
                                                   indicator_function_type,
                                                   num_of_interface_cells,
                                                   eta);

        Pointer<BrinkmanAdvDiffSemiImplicitHierarchyIntegrator> bp_adv_diff_hier_integrator = time_integrator;
        bp_adv_diff_hier_integrator->registerBrinkmanAdvDiffBcHelper(brinkman_adv_diff);

        if (input_db->keyExists("TransportedQuantityForcingFunction"))
        {
            Pointer<CellVariable<NDIM, double> > F_var = new CellVariable<NDIM, double>("F");
            Pointer<CartGridFunction> q_forcing_fcn = new muParserCartGridFunction(
                "q_forcing_fcn",
                app_initializer->getComponentDatabase("TransportedQuantityForcingFunction"),
                grid_geometry);
            time_integrator->registerSourceTerm(F_var, true);
            time_integrator->setSourceTermFunction(F_var, q_forcing_fcn);
            time_integrator->setSourceTerm(q_var, F_var);
        }

        Pointer<CartGridFunction> q_ex = new muParserCartGridFunction(
            "q_ex", app_initializer->getComponentDatabase("TransportedQuantityExactSolutions"), grid_geometry);

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
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int q_idx = var_db->mapVariableAndContextToIndex(q_var, time_integrator->getCurrentContext());
        const int q_ex_cloned_idx = var_db->registerClonedPatchDataIndex(q_var, q_idx);
        const int q_err_cloned_idx = var_db->registerClonedPatchDataIndex(q_var, q_idx);
        const int phi_outer_idx =
            var_db->mapVariableAndContextToIndex(phi_outer_solid_var, time_integrator->getCurrentContext());
        const int phi_cloned_idx = var_db->registerClonedPatchDataIndex(phi_outer_solid_var, phi_outer_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(q_ex_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(q_err_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(phi_cloned_idx, loop_time);
        }
        visit_data_writer->registerPlotQuantity("q_exact", "SCALAR", q_ex_cloned_idx);
        visit_data_writer->registerPlotQuantity("q_error", "SCALAR", q_err_cloned_idx);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);

        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
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

            // Determine the accuracy of the computed solution.
            pout << "\n"
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                 << "Computing error norms.\n\n";

            // Calculate Heaviside function and mask the error indices.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(phi_cloned_idx);

                    for (Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        CellIndex<NDIM> ci(it());
                        double chi = 0.0;
                        for (const auto& phi_solid_var : ls_vars)
                        {
                            const int phi_idx = var_db->mapVariableAndContextToIndex(
                                phi_solid_var, time_integrator->getCurrentContext());
                            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
                            const double phi = (*phi_data)(ci);
                            double Hphi;
                            if (phi > 0.0)
                            {
                                Hphi = 1.0;
                            }
                            else
                            {
                                Hphi = 0.0;
                            }
                            chi += (1.0 - Hphi);
                        }
                        (*H_data)(ci) = 1.0 - chi;
                    }
                }
            }

            double vol = hier_cc_data_ops.integral(phi_cloned_idx, wgt_cc_idx);

            // Get volume weights in the region
            hier_cc_data_ops.multiply(phi_cloned_idx, phi_cloned_idx, wgt_cc_idx);
            double q_integral = hier_cc_data_ops.integral(q_idx, phi_cloned_idx);
            hier_cc_data_ops.addScalar(q_idx, q_idx, -q_integral / vol);
            q_integral = hier_cc_data_ops.integral(q_idx, phi_cloned_idx);
            pout << "q integral" << std::setprecision(10) << q_integral << "\n";

            q_ex->setDataOnPatchHierarchy(q_ex_cloned_idx, q_var, patch_hierarchy, loop_time);
            hier_cc_data_ops.subtract(q_err_cloned_idx, q_idx, q_ex_cloned_idx);

            pout << "Error in q at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(q_err_cloned_idx, phi_cloned_idx)
                 << "\n"
                 << "  L2-norm:  " << hier_cc_data_ops.L2Norm(q_err_cloned_idx, phi_cloned_idx) << "\n"
                 << "  max-norm: " << hier_cc_data_ops.maxNorm(q_err_cloned_idx, phi_cloned_idx) << "\n";

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
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(q_ex_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(q_err_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(phi_cloned_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main
