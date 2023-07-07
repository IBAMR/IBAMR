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
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateInterface.h"

// This function uses a spatially varying Lagrange multiplier.
void
fix_mass_loss_ls_callback_fcn(double /*current_time*/, double /*new_time*/, int /*cycle_num*/, void* ctx)
{
    LevelSetUtilities::LevelSetMassLossFixer* mass_fixer = static_cast<LevelSetUtilities::LevelSetMassLossFixer*>(ctx);
    const double v0 = mass_fixer->getInitialVolume();

    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = mass_fixer->getAdvDiffHierarchyIntegrator();
    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx =
        var_db->mapVariableAndContextToIndex(mass_fixer->getLevelSetVariable(), adv_diff_integrator->getNewContext());
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_integral = 0.0;
    double vol_interface_integral = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(ls_idx);
            Pointer<CellData<NDIM, double> > wgt_data = patch->getPatchData(wgt_cc_idx);

            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d) cell_size *= patch_dx[d];
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double num_cells = 1.0;
            const double alpha = num_cells * cell_size;

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci);
                const double dv = (*wgt_data)(ci);

                double h_prime, h_phi; // smoothed delta and Heaviside functions
                if (phi < -alpha)
                {
                    h_phi = 0.0;
                    h_prime = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    h_phi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                    h_prime = 0.5 / alpha + 1.0 / (2.0 * alpha) * std::cos(M_PI * phi / alpha);
                }
                else
                {
                    h_phi = 1.0;
                    h_prime = 0.0;
                }

                vol_integral += (1.0 - h_phi) * dv;
                vol_interface_integral += h_prime * h_phi * (1.0 - h_phi) * dv;
            }
        }
    }
    const double nmr = IBTK_MPI::sumReduction(vol_integral) - v0;
    const double dnr = IBTK_MPI::sumReduction(vol_interface_integral);
    const double q = nmr / dnr;

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d) cell_size *= patch_dx[d];
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double num_cells = 1.0;
            const double alpha = num_cells * cell_size;

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(ls_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci);

                double h_phi; // smoothed delta
                if (phi < -alpha)
                {
                    h_phi = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    h_phi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                }
                else
                {
                    h_phi = 1.0;
                }

                (*phi_data)(ci) = (*phi_data)(ci) + q * h_phi * (1.0 - h_phi);
            }
        }
    }

    return;
} // fix_mass_loss_ls_callback_fcn

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<AdvDiffHierarchyIntegrator> time_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
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

        // Setup the advection velocity.
        Pointer<FaceVariable<NDIM, double> > u_var = new FaceVariable<NDIM, double>("u");
        Pointer<CartGridFunction> u_fcn =
            new muParserCartGridFunction("u", app_initializer->getComponentDatabase("UFunction"), grid_geometry);
        const bool u_is_div_free = true;
        time_integrator->registerAdvectionVelocity(u_var);
        time_integrator->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        time_integrator->setAdvectionVelocityFunction(u_var, u_fcn);

        // Setup the level set variable Q_var.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        pout << "solving the advection-diffusion equation in "
             << IBAMR::enum_to_string<ConvectiveDifferencingType>(difference_form) << " form.\n";
        Pointer<CellVariable<NDIM, double> > Q_var = new CellVariable<NDIM, double>("Q");
        Pointer<CartGridFunction> Q_init =
            new muParserCartGridFunction("Q", app_initializer->getComponentDatabase("Q_init"), grid_geometry);

        // Boundary condition for the level set variable
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        RobinBcCoefStrategy<NDIM>* Q_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("QBcCoefs"))
        {
            Q_bc_coef =
                new muParserRobinBcCoefs("Q_bc_coef", app_initializer->getComponentDatabase("QBcCoefs"), grid_geometry);
        }

        time_integrator->registerTransportedQuantity(Q_var);
        time_integrator->setAdvectionVelocity(Q_var, u_var);
        time_integrator->setDiffusionCoefficient(Q_var, 0.0);
        time_integrator->setConvectiveDifferencingType(Q_var, difference_form);
        time_integrator->setInitialConditions(Q_var, Q_init);
        time_integrator->setPhysicalBcCoef(Q_var, Q_bc_coef);

        // Reinitialize level set variable
        Pointer<RelaxationLSMethod> level_set_ops = new RelaxationLSMethod(
            "Reinit Level Set", app_initializer->getComponentDatabase("LevelSetReinitialization"));
        LSLocateInterface* ptr_LSLocateInterface = new LSLocateInterface("LSLocateInterface", time_integrator, Q_var);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateInterface));
        level_set_ops->registerPhysicalBoundaryCondition(Q_bc_coef);

        LevelSetUtilities::SetLSProperties* ptr_SetLSProperties =
            new LevelSetUtilities::SetLSProperties("SetLSProperties", level_set_ops);
        time_integrator->registerResetFunction(
            Q_var, &LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(ptr_SetLSProperties));

        // Fix the volume loss due to advection
        NULL_USE(fix_mass_loss_ls_callback_fcn);
        std::vector<Pointer<CellVariable<NDIM, double> > > ls_vars{ Q_var };
        LevelSetUtilities::LevelSetMassLossFixer level_set_fixer(
            "LevelSet  MassLossFixer",
            time_integrator,
            ls_vars,
            app_initializer->getComponentDatabase("LevelSetMassFixer"),
            /*restart*/ true);
        time_integrator->registerPostprocessIntegrateHierarchyCallback(&LevelSetUtilities::fixMassLoss2PhaseFlows,
                                                                       static_cast<void*>(&level_set_fixer));

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Close the restart manager.
        RestartManager::getManager()->closeRestartFile();

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

        // Set the target volume as the gas domain volume (where Q is negative)
        std::vector<double> H_integrals = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(&level_set_fixer);
        level_set_fixer.setInitialVolume(H_integrals[0]);

        // Open stream to save the volume of the two phase and the Lagrange multiplier.
        ofstream vol_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            vol_stream.open("vol.curve", ios_base::out);
            vol_stream.precision(16);
            vol_stream.setf(ios::fixed, ios::floatfield);
            vol_stream << 0.0 << "\t" << H_integrals[0] << "\t" << H_integrals[1] << "\t" << 0.0 << std::endl;
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

            std::vector<double> H_integrals = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(&level_set_fixer);
            if (SAMRAI_MPI::getRank() == 0)
            {
                vol_stream.precision(16);
                vol_stream.setf(ios::fixed, ios::floatfield);
                vol_stream << loop_time << "\t" << H_integrals[0] << "\t" << H_integrals[1] << "\t"
                           << level_set_fixer.getLagrangeMultiplier() << std::endl;
            }

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
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
                pout << "\nWriting restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            vol_stream.close();
        }

        delete ptr_SetLSProperties;
        delete ptr_LSLocateInterface;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
