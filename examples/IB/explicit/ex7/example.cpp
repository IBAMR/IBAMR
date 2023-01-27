// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------
// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <mpi.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

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

struct Lag_Coords
{
    int lag_idx; // The global lag index
    double X;
    double Y;
    double Z;

    // Default constructor
    Lag_Coords()
    : lag_idx(std::numeric_limits<int>::max())
    , X(std::numeric_limits<double>::quiet_NaN())
    , Y(std::numeric_limits<double>::quiet_NaN())
    , Z(std::numeric_limits<double>::quiet_NaN())
    {}

    Lag_Coords(int idx, double x1, double x2, double x3)
    {
        lag_idx = idx;
        X = x1;
        Y = x2;
        Z = x3;
        
    }
};


// NOTE: These values are read in below from the input file. The values listed here are the default values.
static int post_struct_id = 0; // post structure ID number [NOTE: This assumes that the post structures are the first to
                               // be listed in the input file]
static double post_length = 6.0e-3;                   // post length [cm]
static double post_deflection_radius = 2.95e-3;       // maximum post deflection [cm]
static double post_rotational_frequency = 1.0 / .006; // post rotational frequency [1/s]

void
update_target_points(const Pointer<PatchHierarchy<NDIM> >& hierarchy,
                     const LDataManager* const l_data_manager,
                     const double current_time,
                     const double dt)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const std ::pair<int, int>& post_lag_idxs =
        l_data_manager->getLagrangianStructureIndexRange(post_struct_id, finest_ln);
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);

    // update both the local and ghost nodes.
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    const std::vector<LNode*>& ghost_nodes = mesh->getGhostNodes();
    std::vector<LNode*> all_nodes = local_nodes;
    all_nodes.reserve(local_nodes.size() + ghost_nodes.size());
    all_nodes.insert(all_nodes.end(), ghost_nodes.begin(), ghost_nodes.end());
    for (auto* node : all_nodes)
    {
        auto* force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
        if (!force_spec) continue;
        const int lag_idx = node->getLagrangianIndex();
        Point& X_target = force_spec->getTargetPointPosition();
        if (post_lag_idxs.first <= lag_idx && lag_idx < post_lag_idxs.second)
        {
            // We are using the z component of the targeted position to figure out where we are on the post. This
            // fundamentally assumes that the posts are anchored at z=0.  If the posts are anchored at z != 0, then that
            // needs to be taken into account here.
            //
            // This works because we DO NOT change the value of X[2]!

            double h = X_target[2];
            if (h > sqrt(std::numeric_limits<double>::epsilon()))
            {
                // The "slanted post height" is the maximum z value in the slanted configuration. It is determined from
                // the post length and the prescribed deflection radius.
                //
                // Notice that we do not pre-compute this value above because doing so would make it harder to read in
                // user defined post lengths & deflection radii.
                double slanted_post_height =
                    sqrt(post_length * post_length - post_deflection_radius * post_deflection_radius);
                double deflection = post_deflection_radius * h / slanted_post_height;
                double previous_time = current_time - dt;
                X_target[0] = X_target[0] - deflection * cos(2.0 * M_PI * previous_time * post_rotational_frequency) +
                              deflection * cos(2.0 * M_PI * current_time * post_rotational_frequency);
                X_target[1] = X_target[1] - deflection * sin(2.0 * M_PI * previous_time * post_rotational_frequency) +
                              deflection * sin(2.0 * M_PI * current_time * post_rotational_frequency);
            }
        }
    }
}

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
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool is_from_restart = app_initializer->isFromRestart();
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

        // Determine whether to include an dye concentration field.
        const bool simulate_dye_concentration_field =
            input_db->getBoolWithDefault("SIMULATE_DYE_CONCENTRATION_FIELD", false);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        if (simulate_dye_concentration_field)
        {
            const string adv_diff_solver_type = input_db->getStringWithDefault("ADV_DIFF_SOLVER_TYPE", "SEMI_IMPLICIT");
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
            navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        }
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
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
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        if (simulate_dye_concentration_field)
        {
            // Setup the advected and diffused quantity.
            Pointer<CellVariable<NDIM, double> > Q_var = new CellVariable<NDIM, double>("Q");
            adv_diff_integrator->registerTransportedQuantity(Q_var);
            adv_diff_integrator->setDiffusionCoefficient(Q_var, input_db->getDouble("D"));
            if (input_db->keyExists("ConcentrationInitialConditions"))
            {
                Pointer<CartGridFunction> Q_init = new muParserCartGridFunction(
                    "Q_init", app_initializer->getComponentDatabase("ConcentrationInitialConditions"), grid_geometry);
                adv_diff_integrator->setInitialConditions(Q_var, Q_init);
            }
            Pointer<FaceVariable<NDIM, double> > u_adv_var = navier_stokes_integrator->getAdvectionVelocityVariable();
            adv_diff_integrator->setAdvectionVelocity(Q_var, u_adv_var);
            std::unique_ptr<RobinBcCoefStrategy<NDIM> > Q_bc_coefs = nullptr;
            if (periodic_shift.min() == 0)
            {
                Q_bc_coefs.reset(new muParserRobinBcCoefs(
                    "Q_bc_coefs", app_initializer->getComponentDatabase("ConcentrationBcCoefs"), grid_geometry));
                adv_diff_integrator->setPhysicalBcCoef(Q_var, Q_bc_coefs.get());
            }
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write restart data before starting main time integration loop.
        if (dump_restart_data && !is_from_restart)
        {
            pout << "\nWriting restart files...\n\n";
            RestartManager::getManager()->writeRestartFile(restart_dump_dirname, 0);
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

        // figure out how many MPI processes are running and initialize their rank
        const int n_proc = IBTK_MPI::getNodes();
        int rank = IBTK_MPI::getRank();

        // Create Text File which Stores the locations of each IB point.
        std::ofstream fout;

        // Define an MPI type for the Lag_Coords struct
        const int items = NDIM + 1;

        int blocklength[4] = { 1, 1, 1, 1 };
        MPI_Datatype types[4] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
        MPI_Datatype mpi_lag_coords_type;
        MPI_Aint offsets[4];
        offsets[0] = offsetof(Lag_Coords, lag_idx);
        offsets[1] = offsetof(Lag_Coords, X);
        offsets[2] = offsetof(Lag_Coords, Y);
        offsets[3] = offsetof(Lag_Coords, Z);

        MPI_Type_create_struct(items, blocklength, offsets, types, &mpi_lag_coords_type);
        MPI_Type_commit(&mpi_lag_coords_type);

        if (rank == 0)
        {
            fout.open("IBCoordinates.txt");
        }

        // Vector to store the number of nodes on each processor
        std::vector<int> nodes_per_proc(n_proc);

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

            update_target_points(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, dt);

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Print the Eulerian Coordinates of the IB points which constitute the Lagrangian structure.
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<LData> X_data = ib_method_ops->getLDataManager()->getLData("X", finest_ln);
            Vec X_vec = X_data->getVec();
            double* x_vals;
            int ierr = VecGetArray(X_vec, &x_vals);
            IBTK_CHKERRQ(ierr);
            Pointer<LMesh> l_mesh = ib_method_ops->getLDataManager()->getLMesh(finest_ln);
            const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
            int num_nodes = local_nodes.size();

            // Send the number of ib nodes on each processor to the root processor 0
            if (rank != 0)
            {
                MPI_Send(&num_nodes, 1, MPI_INT, 0, 0, IBTK_MPI::getCommunicator());
            }
            else
            {
                nodes_per_proc[0] = num_nodes;

                for (int r = 1; r < n_proc; ++r)
                {
                    MPI_Recv(&nodes_per_proc[r], 1, MPI_INT, r, 0, IBTK_MPI::getCommunicator(), MPI_STATUS_IGNORE);
                }
            }

            IBTK_MPI::barrier();
            // initialize vector to store the nodal coordinates on each processor
            std::vector<Lag_Coords> local_coordinates;
            for (const auto& node : local_nodes)
            {
                const int petsc_idx = node->getLocalPETScIndex();
                const int lag_idx = node->getLagrangianIndex();
                // put the lag idx as the first entry of the local_coordinates vector since the petsc indices may change
                // throughout a computation
                Eigen::Map<Vector3d> X(&x_vals[petsc_idx * NDIM]); // X is a vector containing the local coordinates
                Lag_Coords loc(lag_idx,X[0],X[1],X[2]); //define a temporary Lag_Coords struct 
                local_coordinates.push_back(loc); //add the temporary Lag_Coords struct into local_coordinates
            }

            IBTK_MPI::barrier();
            // Initialize a vector to recieve the coordinates to print to
            std::vector<Lag_Coords> global_coordinates;
            // Specify what's done on the root process
            int total_nodes = 0;
            if (rank == 0)
            {
                total_nodes = std::accumulate(nodes_per_proc.begin(), nodes_per_proc.end(), 0);
                global_coordinates.resize(total_nodes); // resize the receive buffer on the root process to receive the
                                                        // Lag indices and the coordinates

                // Setup the arrays for the counts and displacements
                std::vector<int> count(n_proc);
                std::vector<int> displacements(n_proc);
                displacements[0] = 0;
                count[0] = num_nodes;
                int sum_elem = 0;
                for (int j = 1; j < n_proc; j++)
                {
                    sum_elem = sum_elem + nodes_per_proc[j - 1];
                    count[j] = nodes_per_proc[j];
                    displacements[j] = sum_elem;
                }

                // Gather the data

                MPI_Gatherv(local_coordinates.data(),
                            local_coordinates.size(),
                            mpi_lag_coords_type,
                            global_coordinates.data(),
                            count.data(),
                            displacements.data(),
                            mpi_lag_coords_type,
                            0,
                            IBTK_MPI::getCommunicator());
            }
            else
            { // specify what's done on the non-root processes
                MPI_Gatherv(local_coordinates.data(),
                            local_coordinates.size(),
                            mpi_lag_coords_type,
                            NULL,
                            NULL,
                            NULL,
                            mpi_lag_coords_type,
                            0,
                            IBTK_MPI::getCommunicator());
            }
            IBTK_MPI::barrier();
            // Now if we are on the root process, let's go ahead and sort the global_coordinates vector according to the
            // lag index and print out the coordinates to "IBCoordinates.txt"
            if (rank == 0)
            {
                fout.unsetf(ios_base::showpos);
                fout.setf(ios_base::scientific);
                fout.precision(5);

                 std::sort(global_coordinates.begin(), global_coordinates.end(),
                [](const Lag_Coords &a, const Lag_Coords &b) {return a.lag_idx < b.lag_idx;});

                for (const Lag_Coords& ldata : global_coordinates)
                {
                    fout << loop_time << ", " << ldata.lag_idx << ", " << ldata.X << ", " << ldata.Y << ", " << ldata.Z
                         << "\n";
                }
            }

            ierr = VecRestoreArray(X_vec, &x_vals);
            IBTK_CHKERRQ(ierr);

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
        }

        MPI_Type_free(&mpi_lag_coords_type);
        
        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main
