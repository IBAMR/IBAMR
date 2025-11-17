// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
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
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/ConstraintIBMethod.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBHydrodynamicForceEvaluator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Application objects
#include "IBEELKinematics.h"

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

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
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        // Create odor transport integrator (advection-diffusion solver)
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffHierarchyIntegrator(
            "AdvDiffHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

        // Register advection velocity (provided by Navier-Stokes solver)
        adv_diff_integrator->setAdvectionVelocity(navier_stokes_integrator->getAdvectionVelocityVariable());
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        const int num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        Pointer<ConstraintIBMethod> ib_method_ops = new ConstraintIBMethod(
            "ConstraintIBMethod", app_initializer->getComponentDatabase("ConstraintIBMethod"), num_structures);
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
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create odor concentration initial condition (Gaussian source)
        if (input_db->keyExists("OdorInitialConditions"))
        {
            Pointer<CartGridFunction> C_init = new muParserCartGridFunction(
                "C_init", app_initializer->getComponentDatabase("OdorInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(C_init);
        }

        // Create odor source term (continuous release from spherical source)
        if (input_db->keyExists("OdorSourceTerm"))
        {
            Pointer<CartGridFunction> C_source = new muParserCartGridFunction(
                "C_source", app_initializer->getComponentDatabase("OdorSourceTerm"), grid_geometry);
            adv_diff_integrator->setSourceTerm(C_source);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
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
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create odor concentration boundary conditions
        RobinBcCoefStrategy<NDIM>* C_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("OdorBcCoefs"))
        {
            C_bc_coef = new muParserRobinBcCoefs(
                "C_bc_coef", app_initializer->getComponentDatabase("OdorBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(C_bc_coef);
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
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            adv_diff_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // =============================================================================
        // CREATE KINEMATICS OBJECTS FOR 4 FISH
        // Following official IBAMR eel2d example pattern
        // =============================================================================
        vector<Pointer<ConstraintIBKinematics> > ibkinematics_ops_vec;
        Pointer<ConstraintIBKinematics> ib_kinematics_op;
        
        // struct_0 - FISH-1 (bottom-left, y = -0.3)
        ib_kinematics_op =
            new IBEELKinematics("eel2d_1",
                                app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase("eel2d_1"),
                                ib_method_ops->getLDataManager(),
                                patch_hierarchy);
        ibkinematics_ops_vec.push_back(ib_kinematics_op);
        
        // struct_1 - FISH-2 (bottom-right, y = -0.3)
        ib_kinematics_op =
            new IBEELKinematics("eel2d_2",
                                app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase("eel2d_2"),
                                ib_method_ops->getLDataManager(),
                                patch_hierarchy);
        ibkinematics_ops_vec.push_back(ib_kinematics_op);

        // struct_2 - FISH-3 (top-left, y = +0.3)
        ib_kinematics_op =
            new IBEELKinematics("eel2d_3",
                                app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase("eel2d_3"),
                                ib_method_ops->getLDataManager(),
                                patch_hierarchy);
        ibkinematics_ops_vec.push_back(ib_kinematics_op);

        // struct_3 - FISH-4 (top-right, y = +0.3)
        ib_kinematics_op =
            new IBEELKinematics("eel2d_4",
                                app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase("eel2d_4"),
                                ib_method_ops->getLDataManager(),
                                patch_hierarchy);
        ibkinematics_ops_vec.push_back(ib_kinematics_op);

        // register ConstraintIBKinematics objects with ConstraintIBMethod.
        ib_method_ops->registerConstraintIBKinematics(ibkinematics_ops_vec);
        ib_method_ops->initializeHierarchyOperatorsandData();

        // =============================================================================
        // CREATE HYDRODYNAMIC FORCE EVALUATOR
        // =============================================================================
        double rho_fluid = input_db->getDouble("RHO");
        double mu_fluid = input_db->getDouble("MU");
        double start_time = time_integrator->getIntegratorTime();
        Pointer<IBHydrodynamicForceEvaluator> hydro_force =
            new IBHydrodynamicForceEvaluator("IBHydrodynamicForce", rho_fluid, mu_fluid, start_time, true);

        // =============================================================================
        // REGISTER CONTROL VOLUMES FOR 4 FISH
        // =============================================================================
        
        // Get the initial box position and velocity from input - FISH 1
        const string init_hydro_force_box_db_name = "InitHydroForceBox_0";
        IBTK::Vector3d box_X_lower, box_X_upper, box_init_vel;
        input_db->getDatabase(init_hydro_force_box_db_name)->getDoubleArray("lower_left_corner", &box_X_lower[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name)->getDoubleArray("upper_right_corner", &box_X_upper[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name)->getDoubleArray("init_velocity", &box_init_vel[0], 3);
        hydro_force->registerStructure(box_X_lower, box_X_upper, patch_hierarchy, box_init_vel, 0);

        // Get the initial box position and velocity from input - FISH 2
        const string init_hydro_force_box_db_name_2 = "InitHydroForceBox_1";
        IBTK::Vector3d box_X_lower_2, box_X_upper_2, box_init_vel_2;
        input_db->getDatabase(init_hydro_force_box_db_name_2)->getDoubleArray("lower_left_corner", &box_X_lower_2[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name_2)->getDoubleArray("upper_right_corner", &box_X_upper_2[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name_2)->getDoubleArray("init_velocity", &box_init_vel_2[0], 3);
        hydro_force->registerStructure(box_X_lower_2, box_X_upper_2, patch_hierarchy, box_init_vel_2, 1);

        // Get the initial box position and velocity from input - FISH 3
        const string init_hydro_force_box_db_name_3 = "InitHydroForceBox_2";
        IBTK::Vector3d box_X_lower_3, box_X_upper_3, box_init_vel_3;
        input_db->getDatabase(init_hydro_force_box_db_name_3)->getDoubleArray("lower_left_corner", &box_X_lower_3[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name_3)->getDoubleArray("upper_right_corner", &box_X_upper_3[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name_3)->getDoubleArray("init_velocity", &box_init_vel_3[0], 3);
        hydro_force->registerStructure(box_X_lower_3, box_X_upper_3, patch_hierarchy, box_init_vel_3, 2);

        // Get the initial box position and velocity from input - FISH 4
        const string init_hydro_force_box_db_name_4 = "InitHydroForceBox_3";
        IBTK::Vector3d box_X_lower_4, box_X_upper_4, box_init_vel_4;
        input_db->getDatabase(init_hydro_force_box_db_name_4)->getDoubleArray("lower_left_corner", &box_X_lower_4[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name_4)->getDoubleArray("upper_right_corner", &box_X_upper_4[0], 3);
        input_db->getDatabase(init_hydro_force_box_db_name_4)->getDoubleArray("init_velocity", &box_init_vel_4[0], 3);
        hydro_force->registerStructure(box_X_lower_4, box_X_upper_4, patch_hierarchy, box_init_vel_4, 3);

        // =============================================================================
        // CREATE COM VARIABLES AND SET TORQUE ORIGINS FOR 4 FISH
        // =============================================================================
        std::vector<std::vector<double> > structure_COM = ib_method_ops->getCurrentStructureCOM();
        
        // Create COM variable - FISH 1
        IBTK::Vector3d eel_COM;
        for (int d = 0; d < 3; ++d) eel_COM[d] = structure_COM[0][d];
        hydro_force->setTorqueOrigin(eel_COM, 0);
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 0);

        // Create COM variable - FISH 2
        IBTK::Vector3d eel_COM_2;
        for (int d = 0; d < 3; ++d) eel_COM_2[d] = structure_COM[1][d];
        hydro_force->setTorqueOrigin(eel_COM_2, 1);
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 1);

        // Create COM variable - FISH 3
        IBTK::Vector3d eel_COM_3;
        for (int d = 0; d < 3; ++d) eel_COM_3[d] = structure_COM[2][d];
        hydro_force->setTorqueOrigin(eel_COM_3, 2);
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 2);

        // Create COM variable - FISH 4
        IBTK::Vector3d eel_COM_4;
        for (int d = 0; d < 3; ++d) eel_COM_4[d] = structure_COM[3][d];
        hydro_force->setTorqueOrigin(eel_COM_4, 3);
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 3);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Get velocity and pressure variables from integrator
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);

        const Pointer<Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = navier_stokes_integrator->getCurrentContext();
        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);

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

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        double current_time, new_time;
        double box_disp = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
            current_time = loop_time;

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            loop_time += dt;
            new_time = loop_time;

            // Regrid the hierarchy if necessary.
            if (time_integrator->atRegridPoint()) time_integrator->regridHierarchy();

            // Set the box velocity to nonzero only if the eel has moved sufficiently far.
            IBTK::Vector3d box_vel;
            box_vel.setZero();

            // Velocity due to free-swimming
            std::vector<std::vector<double> > COM_vel = ib_method_ops->getCurrentCOMVelocity();
            for (int d = 0; d < NDIM; ++d) box_vel(d) = COM_vel[0][d];

            int coarsest_ln = 0;
            Pointer<PatchLevel<NDIM> > coarsest_level = patch_hierarchy->getPatchLevel(coarsest_ln);
            const Pointer<CartesianGridGeometry<NDIM> > coarsest_grid_geom = coarsest_level->getGridGeometry();
            const double* const DX = coarsest_grid_geom->getDx();

            // Set the box velocity to ensure that the immersed body remains inside the control volume at all times.
            // If the body's COM has moved 0.9 coarse mesh widths in the x-direction, set the CV velocity such that
            // the CV will translate by 1 coarse mesh width in the direction of swimming (negative x-direction).
            // Otherwise, keep the CV in place by setting its velocity to zero.

            box_disp += box_vel[0] * dt;
            if (abs(box_disp) >= abs(0.9 * DX[0]))
            {
                box_vel.setZero();
                box_vel[0] = -DX[0] / dt;

                box_disp = 0.0;
            }
            else
            {
                box_vel.setZero();
            }

            // Update the location of the box for time n + 1
            hydro_force->updateStructureDomain(box_vel, dt, patch_hierarchy, 0);

            // Compute the momentum of u^n in box n+1 on the newest hierarchy
            hydro_force->computeLaggedMomentumIntegral(
                u_idx, patch_hierarchy, navier_stokes_integrator->getVelocityBoundaryConditions());

            // Advance the hierarchy
            time_integrator->advanceHierarchy(dt);

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Get the momentum of the structures (CRITICAL: MUST BE DECLARED HERE)
            std::vector<std::vector<double> > structure_linear_momentum = ib_method_ops->getStructureMomentum();
            std::vector<std::vector<double> > structure_rotational_momentum =
                ib_method_ops->getStructureRotationalMomentum();

            // Get the momentum of the eel - FISH 1
            IBTK::Vector3d eel_mom, eel_rot_mom;
            eel_mom.setZero();
            eel_rot_mom.setZero();
            for (int d = 0; d < NDIM; ++d) eel_mom[d] = structure_linear_momentum[0][d];
            for (int d = 0; d < 3; ++d) eel_rot_mom[d] = structure_rotational_momentum[0][d];

            // Store the new momenta of the eel - FISH 1
            hydro_force->updateStructureMomentum(eel_mom, eel_rot_mom, 0);

            // Get the momentum of the eel - FISH 2
            IBTK::Vector3d eel_mom_2, eel_rot_mom_2;
            eel_mom_2.setZero();
            eel_rot_mom_2.setZero();
            for (int d = 0; d < NDIM; ++d) eel_mom_2[d] = structure_linear_momentum[1][d];
            for (int d = 0; d < 3; ++d) eel_rot_mom_2[d] = structure_rotational_momentum[1][d];

            // Store the new momenta of the eel - FISH 2
            hydro_force->updateStructureMomentum(eel_mom_2, eel_rot_mom_2, 1);

            // Get the momentum of the eel - FISH 3
            IBTK::Vector3d eel_mom_3, eel_rot_mom_3;
            eel_mom_3.setZero();
            eel_rot_mom_3.setZero();
            for (int d = 0; d < NDIM; ++d) eel_mom_3[d] = structure_linear_momentum[2][d];
            for (int d = 0; d < 3; ++d) eel_rot_mom_3[d] = structure_rotational_momentum[2][d];

            // Store the new momenta of the eel - FISH 3
            hydro_force->updateStructureMomentum(eel_mom_3, eel_rot_mom_3, 2);

            // Get the momentum of the eel - FISH 4
            IBTK::Vector3d eel_mom_4, eel_rot_mom_4;
            eel_mom_4.setZero();
            eel_rot_mom_4.setZero();
            for (int d = 0; d < NDIM; ++d) eel_mom_4[d] = structure_linear_momentum[3][d];
            for (int d = 0; d < 3; ++d) eel_rot_mom_4[d] = structure_rotational_momentum[3][d];

            // Store the new momenta of the eel - FISH 4
            hydro_force->updateStructureMomentum(eel_mom_4, eel_rot_mom_4, 3);

            // Evaluate hydrodynamic force on the eel.
            hydro_force->computeHydrodynamicForce(u_idx,
                                                  p_idx,
                                                  /*f_idx*/ -1,
                                                  patch_hierarchy,
                                                  dt,
                                                  navier_stokes_integrator->getVelocityBoundaryConditions(),
                                                  navier_stokes_integrator->getPressureBoundaryConditions());

            // Print the drag and torque
            hydro_force->postprocessIntegrateData(current_time, new_time);

            // Update CV plot data - FISH 1
            hydro_force->updateStructurePlotData(patch_hierarchy, 0);

            // Update CV plot data - FISH 2
            hydro_force->updateStructurePlotData(patch_hierarchy, 1);

            // Update CV plot data - FISH 3
            hydro_force->updateStructurePlotData(patch_hierarchy, 2);

            // Update CV plot data - FISH 4
            hydro_force->updateStructurePlotData(patch_hierarchy, 3);

            // Set the torque origin for the next time step - FISH 1
            structure_COM = ib_method_ops->getCurrentStructureCOM();
            for (int d = 0; d < 3; ++d) eel_COM[d] = structure_COM[0][d];

            // Set the torque evaluation axis to point from newest COM - FISH 1
            hydro_force->setTorqueOrigin(eel_COM, 0);

            // Set the torque origin for the next time step - FISH 2
            for (int d = 0; d < 3; ++d) eel_COM_2[d] = structure_COM[1][d];

            // Set the torque evaluation axis to point from newest COM - FISH 2
            hydro_force->setTorqueOrigin(eel_COM_2, 1);

            // Set the torque origin for the next time step - FISH 3
            for (int d = 0; d < 3; ++d) eel_COM_3[d] = structure_COM[2][d];

            // Set the torque evaluation axis to point from newest COM - FISH 3
            hydro_force->setTorqueOrigin(eel_COM_3, 2);

            // Set the torque origin for the next time step - FISH 4
            for (int d = 0; d < 3; ++d) eel_COM_4[d] = structure_COM[3][d];

            // Set the torque evaluation axis to point from newest COM - FISH 4
            hydro_force->setTorqueOrigin(eel_COM_4, 3);

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
                            navier_stokes_integrator,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete C_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            LDataManager* l_data_manager,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
    Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    VecView(X_lag_vec, viewer);
    PetscViewerDestroy(&viewer);
    VecDestroy(&X_lag_vec);
    return;
} // output_data