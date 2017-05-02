// Filename: example.cpp
// Modified by Amneet Bhalla from IBAMR examples.

// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
#include <ibamr/ConstraintIBMethod.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBHydrodynamicForceEvaluator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibtk/AppInitializer.h>
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
bool
run_example(int argc, char* argv[])
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
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Create ConstraintIBKinematics objects
        vector<Pointer<ConstraintIBKinematics> > ibkinematics_ops_vec;
        Pointer<ConstraintIBKinematics> ib_kinematics_op;
        // struct_0
        ib_kinematics_op =
            new IBEELKinematics("eel2d",
                                app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase("eel2d"),
                                ib_method_ops->getLDataManager(),
                                patch_hierarchy);
        ibkinematics_ops_vec.push_back(ib_kinematics_op);

        // register ConstraintIBKinematics objects with ConstraintIBMethod.
        ib_method_ops->registerConstraintIBKinematics(ibkinematics_ops_vec);
        ib_method_ops->initializeHierarchyOperatorsandData();

        // Create hydrodynamic force evaluator object.
        double rho_fluid = input_db->getDouble("RHO");
        double mu_fluid = input_db->getDouble("MU");
        double start_time = time_integrator->getIntegratorTime();
        Pointer<IBHydrodynamicForceEvaluator> hydro_force =
            new IBHydrodynamicForceEvaluator("IBHydrodynamicForce", rho_fluid, mu_fluid, start_time, true);

        // Get the initial box position and velocity from input
        const string init_hydro_force_box_db_name = "InitHydroForceBox_0";
        SAMRAI::tbox::Array<double> box_X_lower_array, box_X_upper_array, box_init_vel_array;
        IBTK::Vector3d box_X_lower, box_X_upper, box_init_vel;

        box_X_lower_array = input_db->getDatabase(init_hydro_force_box_db_name)->getDoubleArray("lower_left_corner");
        box_X_upper_array = input_db->getDatabase(init_hydro_force_box_db_name)->getDoubleArray("upper_right_corner");
        box_init_vel_array = input_db->getDatabase(init_hydro_force_box_db_name)->getDoubleArray("init_velocity");
        for (int d = 0; d < 3; ++d)
        {
            box_X_lower[d] = box_X_lower_array[d];
            box_X_upper[d] = box_X_upper_array[d];
            box_init_vel[d] = box_init_vel_array[d];
        }
        hydro_force->registerStructure(box_X_lower, box_X_upper, patch_hierarchy, box_init_vel, 0);

        // Create COM variable
        std::vector<std::vector<double> > structure_COM = ib_method_ops->getStructureCOM();
        IBTK::Vector3d eel_COM;
        for (int d = 0; d < 3; ++d) eel_COM[d] = structure_COM[0][d];

        // Set the torque evaluation axis to point from newest COM
        hydro_force->setTorqueOrigin(eel_COM, 0);

        // Register optional plot data
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 0);

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
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
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

            // Regrid the hierarchy if necessary
            if (time_integrator->atRegridPoint()) time_integrator->regridHierarchy();

            // Set the box velocity to nonzero only if the eel has moved sufficiently far
            IBTK::Vector3d box_vel;
            box_vel.setZero();
            // Velocity due to free-swimming
            std::vector<std::vector<double> > COM_vel = ib_method_ops->getCurrentCOMVelocity();
            for (int d = 0; d < NDIM; ++d) box_vel(d) = COM_vel[0][d];

            int coarsest_ln = 0;
            Pointer<PatchLevel<NDIM> > coarsest_level = patch_hierarchy->getPatchLevel(coarsest_ln);
            const Pointer<CartesianGridGeometry<NDIM> > coarsest_grid_geom = coarsest_level->getGridGeometry();
            const double* const DX = coarsest_grid_geom->getDx();
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

            // Get the momentum of the eel
            IBTK::Vector3d eel_mom, eel_rot_mom;
            eel_mom.setZero();
            eel_rot_mom.setZero();
            std::vector<std::vector<double> > structure_linear_momentum = ib_method_ops->getStructureMomentum();
            for (int d = 0; d < NDIM; ++d) eel_mom[d] = structure_linear_momentum[0][d];
            std::vector<std::vector<double> > structure_rotational_momentum =
                ib_method_ops->getStructureRotationalMomentum();
            for (int d = 0; d < 3; ++d) eel_rot_mom[d] = structure_rotational_momentum[0][d];

            // Store the new momenta of the eel
            hydro_force->updateStructureMomentum(eel_mom, eel_rot_mom, 0);

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

            // Update CV plot data
            hydro_force->updateStructurePlotData(patch_hierarchy, 0);

            // Set the torque origin for the next time step
            structure_COM = ib_method_ops->getStructureCOM();
            for (int d = 0; d < 3; ++d) eel_COM[d] = structure_COM[0][d];

            // Set the torque evaluation axis to point from newest COM
            hydro_force->setTorqueOrigin(eel_COM, 0);

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

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return true;
} // run_example

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
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
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
