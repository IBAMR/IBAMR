// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/RelaxationLSMethod.h>

#include "ibtk/IndexUtilities.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateCircularInterface.h"
#include "SetFluidProperties.h"
#include "SetLSProperties.h"
#include "VelocityInitialCondition.h"

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

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
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
        Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            time_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "SEMI_IMPLICIT");
        if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: SEMI_IMPLICIT");
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

        // Setup level set information
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        const string& ls_name = "level_set";
        Pointer<CellVariable<NDIM, double> > phi_var = new CellVariable<NDIM, double>(ls_name);
        adv_diff_integrator->registerTransportedQuantity(phi_var);
        adv_diff_integrator->setDiffusionCoefficient(phi_var, 0.0);
        // Set the advection velocity of the bubble.
        adv_diff_integrator->setAdvectionVelocity(phi_var, time_integrator->getAdvectionVelocityVariable());

        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateCircularInterface* ptr_LSLocateCircularInterface =
            new LSLocateCircularInterface("LSLocateCircularInterface", adv_diff_integrator, phi_var, &circle);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateCircularInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateCircularInterface));
        SetLSProperties* ptr_SetLSProperties = new SetLSProperties("SetLSProperties", NULL, level_set_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var, &callSetGasLSCallbackFunction, static_cast<void*>(ptr_SetLSProperties));

        // LS initial conditions
        if (input_db->keyExists("LevelSetInitialConditions"))
        {
            Pointer<CartGridFunction> phi_init = new muParserCartGridFunction(
                "phi_init", app_initializer->getComponentDatabase("LevelSetInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(phi_var, phi_init);
        }

        // Setup the INS maintained material properties.
        Pointer<Variable<NDIM> > rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariable<NDIM, double>("rho");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho");
        }

        time_integrator->registerMassDensityVariable(rho_var);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_inside = input_db->getDouble("RHO_I");
        const double rho_outside = input_db->getDouble("RHO_O");
        const double mu_inside = input_db->getDouble("MU_I");
        const double mu_outside = input_db->getDouble("MU_O");
        const int ls_reinit_interval = input_db->getInteger("LS_REINIT_INTERVAL");
        const double num_interface_cells = input_db->getDouble("NUM_INTERFACE_CELLS");
        std::vector<double> inside_velocity(NDIM), outside_velocity(NDIM);
        input_db->getDoubleArray("INSIDE_VELOCITY", &inside_velocity[0], NDIM);
        input_db->getDoubleArray("OUTSIDE_VELOCITY", &outside_velocity[0], NDIM);

        // Callback functions can either be registered with the NS integrator, or the advection-diffusion integrator
        SetFluidProperties* ptr_SetFluidProperties = new SetFluidProperties("SetFluidProperties",
                                                                            adv_diff_integrator,
                                                                            phi_var,
                                                                            rho_outside,
                                                                            rho_inside,
                                                                            mu_outside,
                                                                            mu_inside,
                                                                            ls_reinit_interval,
                                                                            num_interface_cells);
        time_integrator->registerResetFluidDensityFcn(&callSetFluidDensityCallbackFunction,
                                                      static_cast<void*>(ptr_SetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&callSetFluidViscosityCallbackFunction,
                                                        static_cast<void*>(ptr_SetFluidProperties));

        // Create Eulerian initial condition specification objects.
        Pointer<CartGridFunction> u_init =
            new VelocityInitialCondition("u_init", num_interface_cells, inside_velocity, outside_velocity, circle);
        time_integrator->registerVelocityInitialConditions(u_init);

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            time_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
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

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(phi_var, phi_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("RhoBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("RhoBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("MuBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("MuBcCoefs"), grid_geometry);
            time_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

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

        // Get the probe points from the input file
        Pointer<Database> probe_db = input_db->getDatabase("ProbePoints");
        const int num_probes = (probe_db->getAllKeys()).getSize();
        std::vector<std::vector<double> > probe_points;
        probe_points.resize(num_probes);
        for (int i = 0; i < num_probes; ++i)
        {
            std::string probe_name = "probe_" + Utilities::intToString(i);
            probe_points[i].resize(NDIM);
            probe_db->getDoubleArray(probe_name, &probe_points[i][0], NDIM);
        }
        std::ofstream out("output");
        out << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        // File to write to for fluid mass data
        ofstream mass_file;
        mass_file.open("mass_fluid.txt");
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

            // Compute the fluid mass in the domain from interpolated density
            const int rho_ins_idx = time_integrator->getLinearOperatorRhoPatchDataIndex();
#if !defined(NDEBUG)
            TBOX_ASSERT(rho_ins_idx >= 0);
#endif
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            HierarchySideDataOpsReal<NDIM, double> hier_rho_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const double mass_fluid = hier_rho_data_ops.integral(rho_ins_idx, wgt_sc_idx);

            // Write to file
            if (!IBTK_MPI::getRank())
            {
                mass_file << std::setprecision(13) << loop_time << "\t" << mass_fluid << std::endl;
            }

            // Print out the x-velocity values at probe locations
            // Note that it will print the nearest cell center
            // Max reduction over u_val array to ensure that only processor 0 prints
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            int u_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                             time_integrator->getCurrentContext());
            std::vector<double> u_val(num_probes);
            for (int i = 0; i < num_probes; ++i)
            {
                u_val[i] = -std::numeric_limits<double>::max();
                bool found_point_in_patch = false;
                for (int ln = finest_ln; ln >= coarsest_ln; --ln)
                {
                    // Get the cell index for this point
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    const CellIndex<NDIM> cell_idx =
                        IndexUtilities::getCellIndex(&probe_points[i][0], level->getGridGeometry(), level->getRatio());
                    const int axis = 0;
                    const SideIndex<NDIM> side_idx(cell_idx, axis, SideIndex<NDIM>::Lower);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        const Box<NDIM>& patch_box = patch->getBox();
                        const bool contains_probe = patch_box.contains(side_idx);
                        if (!contains_probe) continue;
                        // Get the level set value at this particular cell and print to stream
                        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                        u_val[i] = (*u_data)(side_idx);
                        found_point_in_patch = true;
                        if (found_point_in_patch) break;
                    }
                    if (found_point_in_patch) break;
                }
            }
            IBTK_MPI::maxReduction(&u_val[0], num_probes);
            if (!IBTK_MPI::getRank())
            {
                for (int i = 0; i < num_probes; ++i)
                {
                    out << "x-velocity near point (" << probe_points[i][0] << "," << probe_points[i][1]
#if (NDIM == 3)
                        << "," << probe_points[i][2]
#endif
                        << ") = " << std::setprecision(5) << u_val[i] << std::endl;
                }
            }

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

        out << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        // Close file
        mass_file.close();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        // Cleanup other dumb pointers
        delete ptr_SetLSProperties;
        delete ptr_SetFluidProperties;
        delete ptr_LSLocateCircularInterface;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
