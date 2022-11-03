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

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// The tether penalty functions each require some data that is set in the
// input file. This data is passed to each object through the void *ctx
// context data pointer. Here we collect all relevant tether data in a struct:
struct ElasticityData
{
    const double c1_s;
    const double kappa_s;
    const double mu_s;
    const double lambda_s;

    ElasticityData(Pointer<Database> input_db)
        : c1_s(input_db->getDouble("C1_S")),
          kappa_s(input_db->getDouble("KAPPA_S")),
          mu_s(input_db->getDouble("MU_S")),
          lambda_s(input_db->getDouble("LAMBDA_S"))
    {
    }
};

// Tether (penalty) force function for the solid block.
void
block_tether_force_function(VectorValue<double>& F,
                            const TensorValue<double>& /*FF*/,
                            const libMesh::Point& X,
                            const libMesh::Point& s,
                            Elem* const /*elem*/,
                            const std::vector<const std::vector<double>*>& /*var_data*/,
                            const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double /*time*/,
                            void* ctx)
{
    const ElasticityData* const elasticity_data = reinterpret_cast<ElasticityData*>(ctx);

    F = elasticity_data->kappa_s * (s - X);
    return;
} // block_tether_force_function

// Tether (penalty) force function for the thin beam.
void
beam_tether_force_function(VectorValue<double>& F,
                           const TensorValue<double>& /*FF*/,
                           const libMesh::Point& X,
                           const libMesh::Point& s,
                           Elem* const /*elem*/,
                           const std::vector<const std::vector<double>*>& /*var_data*/,
                           const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                           double /*time*/,
                           void* ctx)
{
    const double r = sqrt((s(0) - 0.2) * (s(0) - 0.2) + (s(1) - 0.2) * (s(1) - 0.2));
    if (r <= 0.05)
    {
        const ElasticityData* const elasticity_data = reinterpret_cast<ElasticityData*>(ctx);
        F = elasticity_data->kappa_s * (s - X);
    }
    else
    {
        F.zero();
    }
    return;
} // beam_tether_force_function

// (Penalty) stress tensor function for the solid block.
void
block_PK1_stress_function(TensorValue<double>& PP,
                          const TensorValue<double>& FF,
                          const libMesh::Point& /*X*/,
                          const libMesh::Point& /*s*/,
                          Elem* const /*elem*/,
                          const std::vector<const std::vector<double>*>& /*var_data*/,
                          const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                          double /*time*/,
                          void* ctx)
{
    const ElasticityData* const elasticity_data = reinterpret_cast<ElasticityData*>(ctx);

    PP = 2.0 * elasticity_data->c1_s * (FF - tensor_inverse_transpose(FF, NDIM));
    return;
} // block_PK1_stress_function

void
beam_PK1_stress_function(TensorValue<double>& PP,
                         const TensorValue<double>& FF,
                         const libMesh::Point& /*X*/,
                         const libMesh::Point& /*s*/,
                         Elem* const /*elem*/,
                         const std::vector<const std::vector<double>*>& /*var_data*/,
                         const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                         double /*time*/,
                         void* ctx)
{
    const ElasticityData* const elasticity_data = reinterpret_cast<ElasticityData*>(ctx);
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> EE = 0.5 * (CC - II);
    const TensorValue<double> SS = elasticity_data->lambda_s * EE.tr() * II + 2.0 * elasticity_data->mu_s * EE;
    PP = FF * SS;
    return;
} // beam_PK1_stress_function
} // namespace ModelData
using namespace ModelData;

// Function prototypes
static ofstream drag_stream, lift_stream, A_x_posn_stream, A_y_posn_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      const IBFEMethod* ib_method_ops,
                      Mesh& beam_mesh,
                      EquationSystems* beam_equation_systems,
                      Mesh& block_mesh,
                      EquationSystems* block_equation_systems,
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
        const string block_exodus_filename = app_initializer->getExodusIIFilename("block");
        const string beam_exodus_filename = app_initializer->getExodusIIFilename("beam");

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
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;

        Mesh block_mesh(init.comm(), NDIM);
        string block_elem_type = input_db->getString("BLOCK_ELEM_TYPE");
        const double R = 0.05;
        if (block_elem_type == "TRI3" || block_elem_type == "TRI6")
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                block_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(block_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(block_elem_type);
            triangle.desired_area() = sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
            block_mesh.prepare_for_use();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = ceil(2.0 * M_PI * R / ds);
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(block_mesh, R, r, Utility::string_to_enum<ElemType>(block_elem_type));
        }
        for (MeshBase::node_iterator n_it = block_mesh.nodes_begin(); n_it != block_mesh.nodes_end(); ++n_it)
        {
            Node& n = **n_it;
            n(0) += 0.2;
            n(1) += 0.2;
        }

        Mesh beam_mesh(init.comm(), NDIM);
        string beam_elem_type = input_db->getString("BEAM_ELEM_TYPE");
        MeshTools::Generation::build_square(beam_mesh,
                                            ceil(0.4 / ds),
                                            ceil(0.02 / ds),
                                            0.2,
                                            0.6,
                                            0.19,
                                            0.21,
                                            Utility::string_to_enum<ElemType>(beam_elem_type));
        beam_mesh.prepare_for_use();

        vector<MeshBase*> meshes(2);
        meshes[0] = &block_mesh;
        meshes[1] = &beam_mesh;

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
                           meshes,
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
        ElasticityData elasticity_data(input_db);
        void* const elasticity_data_ptr = reinterpret_cast<void*>(&elasticity_data);
        IBFEMethod::LagBodyForceFcnData block_tether_force_data(
            block_tether_force_function, std::vector<IBTK::SystemData>(), elasticity_data_ptr);
        IBFEMethod::PK1StressFcnData block_PK1_stress_data(
            block_PK1_stress_function, std::vector<IBTK::SystemData>(), elasticity_data_ptr);
        ib_method_ops->registerLagBodyForceFunction(block_tether_force_data, 0);
        ib_method_ops->registerPK1StressFunction(block_PK1_stress_data, 0);
        string block_kernel_fcn = input_db->getStringWithDefault("BLOCK_KERNEL_FUNCTION", "PIECEWISE_LINEAR");
        FEDataManager::InterpSpec block_interp_spec = ib_method_ops->getDefaultInterpSpec();
        block_interp_spec.kernel_fcn = block_kernel_fcn;
        ib_method_ops->setInterpSpec(block_interp_spec, 0);
        FEDataManager::SpreadSpec block_spread_spec = ib_method_ops->getDefaultSpreadSpec();
        block_spread_spec.kernel_fcn = block_kernel_fcn;
        ib_method_ops->setSpreadSpec(block_spread_spec, 0);

        IBFEMethod::LagBodyForceFcnData beam_tether_force_data(
            beam_tether_force_function, std::vector<IBTK::SystemData>(), elasticity_data_ptr);
        IBFEMethod::PK1StressFcnData beam_PK1_stress_data(
            beam_PK1_stress_function, std::vector<IBTK::SystemData>(), elasticity_data_ptr);
        ib_method_ops->registerLagBodyForceFunction(beam_tether_force_data, 1);
        ib_method_ops->registerPK1StressFunction(beam_PK1_stress_data, 1);
        string beam_kernel_fcn = input_db->getStringWithDefault("BEAM_KERNEL_FUNCTION", "IB_3");
        FEDataManager::InterpSpec beam_interp_spec = ib_method_ops->getDefaultInterpSpec();
        beam_interp_spec.kernel_fcn = beam_kernel_fcn;
        ib_method_ops->setInterpSpec(beam_interp_spec, 1);
        FEDataManager::SpreadSpec beam_spread_spec = ib_method_ops->getDefaultSpreadSpec();
        beam_spread_spec.kernel_fcn = beam_kernel_fcn;
        ib_method_ops->setSpreadSpec(beam_spread_spec, 1);

        ib_method_ops->initializeFEEquationSystems();
        EquationSystems* block_equation_systems = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* beam_equation_systems = ib_method_ops->getFEDataManager(1)->getEquationSystems();

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
        std::unique_ptr<ExodusII_IO> block_exodus_io(uses_exodus ? new ExodusII_IO(block_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> beam_exodus_io(uses_exodus ? new ExodusII_IO(beam_mesh) : NULL);

        // Check to see if this is a restarted run to append current exodus files
        if (uses_exodus)
        {
            const bool from_restart = RestartManager::getManager()->isFromRestart();
            block_exodus_io->append(from_restart);
            beam_exodus_io->append(from_restart);
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
                block_exodus_io->write_timestep(
                    block_exodus_filename, *block_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                beam_exodus_io->write_timestep(
                    beam_exodus_filename, *beam_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save lift and drag coefficients.
        if (IBTK_MPI::getRank() == 0)
        {
            drag_stream.open("C_D.curve", ios_base::out | ios_base::trunc);
            lift_stream.open("C_L.curve", ios_base::out | ios_base::trunc);
            A_x_posn_stream.open("A_x.curve", ios_base::out | ios_base::trunc);
            A_y_posn_stream.open("A_y.curve", ios_base::out | ios_base::trunc);
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
                    block_exodus_io->write_timestep(block_exodus_filename,
                                                    *block_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                    beam_exodus_io->write_timestep(
                        beam_exodus_filename, *beam_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
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
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting state data...\n\n";
                postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 ib_method_ops,
                                 beam_mesh,
                                 beam_equation_systems,
                                 block_mesh,
                                 block_equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }
        }

        // Close the logging streams.
        if (IBTK_MPI::getRank() == 0)
        {
            drag_stream.close();
            lift_stream.close();
            A_x_posn_stream.close();
            A_y_posn_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                 const IBFEMethod* const ib_method_ops,
                 Mesh& beam_mesh,
                 EquationSystems* beam_equation_systems,
                 Mesh& block_mesh,
                 EquationSystems* block_equation_systems,
                 const int /*iteration_num*/,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{
    double F_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
    Mesh* mesh[2] = { &beam_mesh, &block_mesh };
    EquationSystems* equation_systems[2] = { beam_equation_systems, block_equation_systems };
    for (unsigned int k = 0; k < 2; ++k)
    {
        System& F_system = equation_systems[k]->get_system<System>(ib_method_ops->getForceSystemName());
        NumericVector<double>* F_vec = F_system.solution.get();
        NumericVector<double>* F_ghost_vec = F_system.current_local_solution.get();
        copy_and_synch(*F_vec, *F_ghost_vec);
        DofMap& F_dof_map = F_system.get_dof_map();
        std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);
        std::unique_ptr<FEBase> fe(FEBase::build(NDIM, F_dof_map.variable_type(0)));
        std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi = fe->get_phi();
        const std::vector<double>& JxW = fe->get_JxW();
        boost::multi_array<double, 2> F_node;
        const auto el_begin = mesh[k]->active_local_elements_begin();
        const auto el_end = mesh[k]->active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            const auto elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_dof_map.dof_indices(elem, F_dof_indices[d], d);
            }
            const int n_qp = qrule->n_points();
            const int n_basis = static_cast<int>(F_dof_indices[0].size());
            get_values_for_interpolation(F_node, *F_ghost_vec, F_dof_indices);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                for (int k = 0; k < n_basis; ++k)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_integral[d] += F_node[k][d] * phi[k][qp] * JxW[qp];
                    }
                }
            }
        }
    }
    IBTK_MPI::sumReduction(F_integral, NDIM);
    if (IBTK_MPI::getRank() == 0)
    {
        drag_stream.precision(12);
        drag_stream.setf(ios::fixed, ios::floatfield);
        drag_stream << loop_time << " " << -F_integral[0] << endl;
        lift_stream.precision(12);
        lift_stream.setf(ios::fixed, ios::floatfield);
        lift_stream << loop_time << " " << -F_integral[1] << endl;
    }

    System& X_system = beam_equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
    NumericVector<double>* X_vec = X_system.solution.get();
    std::unique_ptr<NumericVector<Number> > X_serial_vec = NumericVector<Number>::build(X_vec->comm());
    X_serial_vec->init(X_vec->size(), true, SERIAL);
    X_vec->localize(*X_serial_vec);
    DofMap& X_dof_map = X_system.get_dof_map();
    vector<unsigned int> vars(2);
    vars[0] = 0;
    vars[1] = 1;
    MeshFunction X_fcn(*beam_equation_systems, *X_serial_vec, X_dof_map, vars);
    X_fcn.init();
    DenseVector<double> X_A(2);
    X_fcn(libMesh::Point(0.6, 0.2, 0), 0.0, X_A);
    if (IBTK_MPI::getRank() == 0)
    {
        A_x_posn_stream.precision(12);
        A_x_posn_stream.setf(ios::fixed, ios::floatfield);
        A_x_posn_stream << loop_time << " " << X_A(0) << endl;
        A_y_posn_stream.precision(12);
        A_y_posn_stream.setf(ios::fixed, ios::floatfield);
        A_y_posn_stream << loop_time << " " << X_A(1) << endl;
    }
    return;
} // postprocess_data
