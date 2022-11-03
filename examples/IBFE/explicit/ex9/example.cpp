// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
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

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
System *x_solid_system, *u_solid_system;
void
tether_force_function(VectorValue<double>& F,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the boundary mesh
                      const libMesh::Point& X_bndry, // X_bndry gives reference coordinates on the boundary mesh
                      Elem* const elem,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double /*time*/,
                      void* /*ctx*/)
{
    // tether_force_function() is called on elements of the boundary mesh.  Here
    // we look up the element in the solid mesh that the current boundary
    // element was extracted from.
    const Elem* const interior_parent = elem->interior_parent();

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.
    std::vector<double> x_solid(NDIM, 0.0), u_solid(NDIM, 0.0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_solid_system->point_value(d, X_bndry, interior_parent);
    }

    // Look up the velocity of the boundary mesh.
    const std::vector<double>& u_bndry = *var_data[0];

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s * (x_solid[d] - x_bndry(d)) + eta_s * (u_solid[d] - u_bndry[d]);
    }
    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

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
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

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
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string exodus_solid_filename = "solid_output.ex2"; // app_initializer->getExodusIIFilename();
        const string exodus_bndry_filename = "bndry_output.ex2"; // app_initializer->getExodusIIFilename();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 0.5;
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        // Ensure nodes on the surface are on the analytic boundary.
        MeshBase::element_iterator el_end = solid_mesh.elements_end();
        for (MeshBase::element_iterator el = solid_mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (!at_mesh_bdry) continue;
                for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                {
                    if (!elem->is_node_on_side(k, side)) continue;
                    Node& n = elem->node_ref(k);
                    n = R * n.unit();
                }
            }
        }
        solid_mesh.prepare_for_use();

        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &boundary_mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
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

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(ib_method_ops->getVelocitySystemName(), vars));
        IBFEMethod::LagBodyForceFcnData body_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagBodyForceFunction(body_fcn_data);
        EquationSystems* bndry_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        // Setup solid systems.
        std::unique_ptr<EquationSystems> solid_equation_systems(new EquationSystems(solid_mesh));
        x_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("position");
        u_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("velocity");
        Order order = SECOND;
        FEFamily family = LAGRANGE;
        for (int d = 0; d < NDIM; ++d)
        {
            x_solid_system->add_variable("X_" + std::to_string(d), order, family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            u_solid_system->add_variable("U_" + std::to_string(d), order, family);
        }
        solid_equation_systems->init();

        // Set up the position vector.
        //
        // \todo There needs to be an API so that this can be handled by the IBFEMethod class implementation.
        {
            MeshBase& mesh = solid_equation_systems->get_mesh();
            System& X_system = solid_equation_systems->get_system("position");
            const unsigned int X_sys_num = X_system.number();
            NumericVector<double>& X_coords = *X_system.solution;
            for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
            {
                Node* n = *it;
                if (n->n_vars(X_sys_num))
                {
                    TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
                    const libMesh::Point& X = *n;
                    libMesh::Point x = X;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int dof_index = n->dof_number(X_sys_num, d, 0);
                        X_coords.set(dof_index, x(d));
                    }
                }
            }
            X_coords.close();
            X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
            copy_and_synch(*X_system.solution, *X_system.current_local_solution);
        }

        x_solid_system->assemble_before_solve = false;
        x_solid_system->assemble();

        u_solid_system->assemble_before_solve = false;
        u_solid_system->assemble();

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
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_solid_io(uses_exodus ? new ExodusII_IO(solid_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_io(uses_exodus ? new ExodusII_IO(boundary_mesh) : NULL);

        // Check to see if this is a restarted run to append current exodus files
        if (uses_exodus)
        {
            const bool from_restart = RestartManager::getManager()->isFromRestart();
            exodus_solid_io->append(from_restart);
            exodus_bndry_io->append(from_restart);
        }

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_solid_io->write_timestep(
                    exodus_solid_filename, *solid_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_bndry_io->write_timestep(
                    exodus_bndry_filename, *bndry_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            // Setup the position and velocity vector.
            //
            // \todo There needs to be an API so that this can be handled by the IBFEMethod class implementation.
            {
                DenseVector<double> U(NDIM);
                U(1) = -1.0;
                MeshBase& mesh = solid_equation_systems->get_mesh();
                System& X_system = solid_equation_systems->get_system("position");
                const unsigned int X_sys_num = X_system.number();
                NumericVector<double>& X_coords = *X_system.solution;
                System& U_system = solid_equation_systems->get_system("velocity");
                const unsigned int U_sys_num = U_system.number();
                NumericVector<double>& U_coords = *U_system.solution;
                for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
                {
                    Node* n = *it;
                    if (n->n_vars(X_sys_num))
                    {
                        TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
                        const libMesh::Point& X = *n;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            const int dof_index = n->dof_number(U_sys_num, d, 0);
                            X_coords.set(dof_index, X(d) + loop_time * U(d));
                            U_coords.set(dof_index, U(d));
                        }
                    }
                }
                X_coords.close();
                X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
                copy_and_synch(X_coords, *X_system.current_local_solution);
                U_coords.close();
                U_system.get_dof_map().enforce_constraints_exactly(U_system, &U_coords);
                copy_and_synch(U_coords, *U_system.current_local_solution);
            }

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
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_solid_io->write_timestep(exodus_solid_filename,
                                                    *solid_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                    exodus_bndry_io->write_timestep(exodus_bndry_filename,
                                                    *bndry_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main
