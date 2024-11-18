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
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBHydrodynamicForceEvaluator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/KrylovMobilitySolver.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

#include <ibamr/app_namespaces.h>

//////////////////////////////////////////////////////////////////////////////

// Center of mass velocity
void
ConstrainedCOMOuterVel(double /*data_time*/, IBTK::Vector3d& U_com, IBTK::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    U_com[0] = 1.0;
    U_com[1] = 2.0;
    W_com[1] = 5.0;

    return;
} // ConstrainedCOMOuterVel

// Center of mass velocity
void
ConstrainedCOMInnerVel(double /*data_time*/, IBTK::Vector3d& U_com, IBTK::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // ConstrainedCOMInnerVel

void
ConstrainedNodalVel(Vec /*U_k*/, const RigidDOFVector& /*U*/, const IBTK::Vector3d& /*X_com*/, void* /*ctx*/)
{
    // intentionally left blank.
    return;
} // ConstrainedNodalVel

// These forces on the structure are computed when the "above" velocities are prescribed on them.
void
NetExternalForceTorqueOuter(double /*data_time*/, IBTK::Vector3d& F_ext, IBTK::Vector3d& T_ext, void* /*ctx*/)
{
    F_ext << 100.916, 1351.74, 0.000916801;
    T_ext << -2.23858e-09, 58206.1, -1.12567e-09;

    return;
} // NetExternalForceTorqueOuter

void
NetExternalForceTorqueInner(double /*data_time*/, IBTK::Vector3d& F_ext, IBTK::Vector3d& T_ext, void* /*ctx*/)
{
    F_ext << -116.762, -201.82, -0.00128781;
    T_ext << 1.46423e-09, -1393.33, 9.03505e-10;

    return;
} // NetExternalForceTorqueInner

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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Several parts of the code (such as LDataManager) expect mesh files,
    // specified in the input file, to exist in the current working
    // directory. Since tests are run in temporary directories we need to regenerate these input
    // to work. We also create a petsc options file for CIB solvers.
    if (IBTK_MPI::getRank() == 0)
    {
        std::ifstream shell_in_vertex_stream(SOURCE_DIR "/shell_3d_in.vertex");
        std::ofstream shell_in_cwd("shell_3d_in.vertex");
        shell_in_cwd << shell_in_vertex_stream.rdbuf();

        std::ifstream shell_out_vertex_stream(SOURCE_DIR "/shell_3d_out.vertex");
        std::ofstream shell_out_cwd("shell_3d_out.vertex");
        shell_out_cwd << shell_out_vertex_stream.rdbuf();

        std::ifstream petsc_options_stream(SOURCE_DIR "/petsc_options.dat");
        std::ofstream petsc_options_cwd("petsc_options.dat");
        petsc_options_cwd << petsc_options_stream.rdbuf();
    }

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "CIB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // some variants of this test use the SVD, some don't
        const auto inverse_type = IBAMR::string_to_enum<IBAMR::MobilityMatrixInverseType>(
            input_db->getStringWithDefault("mobility_inverse_type", "LAPACK_SVD"));

        // Read default Petsc options
        if (input_db->keyExists("petsc_options_file"))
        {
            std::string petsc_options_file = input_db->getString("petsc_options_file");
            PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, petsc_options_file.c_str(), PETSC_TRUE);
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
        FreeRigidDOFVector outer_free_dofs, inner_free_dofs;

        // Make some DOFs prescribed and some free.
        outer_free_dofs << 0, 0, 1, 1, 1, 1;
        inner_free_dofs << 1, 1, 0, 1, 1, 0;
        ib_method_ops->setSolveRigidBodyVelocity(0, outer_free_dofs);
        ib_method_ops->setSolveRigidBodyVelocity(1, inner_free_dofs);

        ib_method_ops->registerConstrainedVelocityFunction(nullptr, &ConstrainedCOMOuterVel, nullptr, 0);
        ib_method_ops->registerConstrainedVelocityFunction(nullptr, &ConstrainedCOMInnerVel, nullptr, 1);

        ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorqueOuter, nullptr, 0);
        ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorqueInner, nullptr, 1);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

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
                u_bc_coefs[d] = nullptr;
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
            std::string mat_name1 = "struct-1";
            std::string mat_name2 = "struct-2";
            std::vector<std::vector<unsigned> > struct_ids1;
            std::vector<std::vector<unsigned> > struct_ids2;
            std::vector<unsigned> prototype_structs1;
            std::vector<unsigned> prototype_structs2;

            // Dense matrix
            prototype_structs1.push_back(0);
            prototype_structs2.push_back(1);

            struct_ids1.push_back(prototype_structs1);
            struct_ids2.push_back(prototype_structs2);

            DirectMobilitySolver* direct_solvers = nullptr;
            CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(
                nullptr, &direct_solvers, nullptr);

            direct_solvers->registerMobilityMat(
                mat_name1, prototype_structs1, EMPIRICAL, std::make_pair(inverse_type, inverse_type), 0);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name1, struct_ids1);

            int next_proc = 0;
            if (IBTK_MPI::getNodes() > 1) next_proc = 1;
            direct_solvers->registerMobilityMat(
                mat_name2, prototype_structs2, EMPIRICAL, std::make_pair(inverse_type, inverse_type), next_proc);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name2, struct_ids2);
        }
        navier_stokes_integrator->setStokesSolverNeedsInit();

        // Set up the hydro force objects
        double rho_fluid = input_db->getDouble("RHO");
        double mu_fluid = input_db->getDouble("MU");
        double start_time = time_integrator->getIntegratorTime();
        Pointer<IBHydrodynamicForceEvaluator> hydro_force =
            new IBHydrodynamicForceEvaluator("IBHydrodynamicForce", rho_fluid, mu_fluid, start_time, true);

        // Get the initial box position and velocity from input
        const string init_hydro_force_box_out_db_name = "InitHydroForceBox_0";
        IBTK::Vector3d box_X_lower_out, box_X_upper_out, box_init_vel_out;

        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("lower_left_corner", &box_X_lower_out[0], 3);
        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("upper_right_corner", &box_X_upper_out[0], 3);
        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("init_velocity", &box_init_vel_out[0], 3);

        const string init_hydro_force_box_in_db_name = "InitHydroForceBox_1";
        IBTK::Vector3d box_X_lower_in, box_X_upper_in, box_init_vel_in;

        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("lower_left_corner", &box_X_lower_in[0], 3);
        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("upper_right_corner", &box_X_upper_in[0], 3);
        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("init_velocity", &box_init_vel_in[0], 3);

        // Register the control volumes
        hydro_force->registerStructure(box_X_lower_out, box_X_upper_out, patch_hierarchy, box_init_vel_out, 0);
        hydro_force->registerStructure(box_X_lower_in, box_X_upper_in, patch_hierarchy, box_init_vel_in, 1);

        // Set up optional visualization of select boxes
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 0);
        hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, 1);

        // Set the origin of the torque evaluation
        IBTK::Vector3d torque_origin_out, torque_origin_in;
        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("torque_origin", &torque_origin_out[0], 3);
        input_db->getDatabase(init_hydro_force_box_out_db_name)
            ->getDoubleArray("torque_origin", &torque_origin_in[0], 3);

        hydro_force->setTorqueOrigin(torque_origin_out, 0);
        hydro_force->setTorqueOrigin(torque_origin_in, 1);

        // Deallocate initialization objects.
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
        if (dump_postproc_data)
        {
            output_data(patch_hierarchy,
                        ib_method_ops->getLDataManager(),
                        iteration_num,
                        loop_time,
                        postproc_data_dump_dirname);
        }

        // File to write the output data
        ofstream output_file;
        if (!IBTK_MPI::getRank()) output_file.open("output");

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double current_time, new_time;
        double dt = 0.0;

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

            pout << "Advancing hierarchy by timestep size dt = " << dt << "\n";

            if (time_integrator->atRegridPoint()) navier_stokes_integrator->setStokesSolverNeedsInit();
            if (ib_method_ops->flagRegrid())
            {
                time_integrator->regridHierarchy();
                navier_stokes_integrator->setStokesSolverNeedsInit();
            }

            // Update the location of the box for time n + 1
            hydro_force->updateStructureDomain(IBTK::Vector3d::Zero(), dt, patch_hierarchy, 0);
            hydro_force->updateStructureDomain(IBTK::Vector3d::Zero(), dt, patch_hierarchy, 1);

            // Compute the momentum of u^n in box n+1 on the newest hierarchy
            hydro_force->computeLaggedMomentumIntegral(
                u_idx, patch_hierarchy, navier_stokes_integrator->getVelocityBoundaryConditions());

            // Advance the hierarchy
            time_integrator->advanceHierarchy(dt);

            pout << "\n\nNet rigid force and torque on structure 0 is : \n"
                 << ib_method_ops->getNetRigidGeneralizedForce(0) << "\n\n";
            pout << "\n\nNet rigid force and torque on structure 1 is : \n"
                 << ib_method_ops->getNetRigidGeneralizedForce(1) << "\n\n";

            RDV U;
            ib_method_ops->getNewRigidBodyVelocity(0, U);
            pout << "\n\nRigid body velocity of structure 0 is : \n" << U << "\n\n";
            ib_method_ops->getNewRigidBodyVelocity(1, U);
            pout << "\n\nRigid body velocity of structure 1 is : \n" << U << "\n\n";

            if (!IBTK_MPI::getRank())
            {
                output_file << std::setprecision(13) << "\n\nNet rigid force and torque on structure 0 is : \n"
                            << ib_method_ops->getNetRigidGeneralizedForce(0) << "\n\n"
                            << "\n\nNet rigid force and torque on structure 1 is : \n"
                            << ib_method_ops->getNetRigidGeneralizedForce(1) << "\n\n";

                ib_method_ops->getNewRigidBodyVelocity(0, U);
                output_file << "\n\nRigid body velocity of structure 0 is : \n" << U << "\n\n";
                ib_method_ops->getNewRigidBodyVelocity(1, U);
                output_file << "\n\nRigid body velocity of structure 1 is : \n" << U << "\n\n";
            }

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Set the linear and angular momentum of the spheres (zero in Stokes flow)
            IBTK::Vector3d Sphere_mom;
            IBTK::Vector3d Sphere_ang_mom;
            Sphere_mom.setZero();
            Sphere_ang_mom.setZero();

            // Store the new momenta of the sphere
            hydro_force->updateStructureMomentum(Sphere_mom, Sphere_ang_mom, 0);
            hydro_force->updateStructureMomentum(Sphere_mom, Sphere_ang_mom, 1);

            // Evaluate hydrodynamic force on the sphere.
            hydro_force->computeHydrodynamicForce(u_idx,
                                                  p_idx,
                                                  /*f_idx*/ -1,
                                                  patch_hierarchy,
                                                  dt,
                                                  navier_stokes_integrator->getVelocityBoundaryConditions(),
                                                  navier_stokes_integrator->getPressureBoundaryConditions());
            // Post processing call for writing data
            hydro_force->postprocessIntegrateData(current_time, new_time);

            // Update optional visualization of select boxes
            hydro_force->updateStructurePlotData(patch_hierarchy, 0);
            hydro_force->updateStructurePlotData(patch_hierarchy, 1);

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

        // Close file
        if (!IBTK_MPI::getRank()) output_file.close();

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
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
