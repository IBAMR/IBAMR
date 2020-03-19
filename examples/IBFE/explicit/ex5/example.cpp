// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBFESurfaceMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

inline double
kernel(double x)
{
    x += 4.;
    const double x2 = x * x;
    const double x3 = x * x2;
    const double x4 = x * x3;
    const double x5 = x * x4;
    const double x6 = x * x5;
    const double x7 = x * x6;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1984126984126984e-3 * x7;
    else if (x <= 2.)
        return .1111111111111111e-1 * x6 - .1388888888888889e-2 * x7 - .3333333333333333e-1 * x5 +
               .5555555555555556e-1 * x4 - .5555555555555556e-1 * x3 + .3333333333333333e-1 * x2 -
               .1111111111111111e-1 * x + .1587301587301587e-2;
    else if (x <= 3.)
        return .4333333333333333 * x5 - .6666666666666667e-1 * x6 + .4166666666666667e-2 * x7 - 1.500000000000000 * x4 +
               3.055555555555556 * x3 - 3.700000000000000 * x2 + 2.477777777777778 * x - .7095238095238095;
    else if (x <= 4.)
        return 9. * x4 - 1.666666666666667 * x5 + .1666666666666667 * x6 - .6944444444444444e-2 * x7 -
               28.44444444444444 * x3 + 53. * x2 - 54.22222222222222 * x + 23.59047619047619;
    else if (x <= 5.)
        return 96. * x3 - 22.11111111111111 * x4 + 3. * x5 - .2222222222222222 * x6 + .6944444444444444e-2 * x7 -
               245.6666666666667 * x2 + 344. * x - 203.9650793650794;
    else if (x <= 6.)
        return 483.5000000000000 * x2 - 147.0555555555556 * x3 + 26.50000000000000 * x4 - 2.833333333333333 * x5 +
               .1666666666666667 * x6 - .4166666666666667e-2 * x7 - 871.2777777777778 * x + 664.0904761904762;
    else if (x <= 7.)
        return 943.1222222222222 * x - 423.7000000000000 * x2 + 104.9444444444444 * x3 - 15.50000000000000 * x4 +
               1.366666666666667 * x5 - .6666666666666667e-1 * x6 + .1388888888888889e-2 * x7 - 891.1095238095238;
    else if (x <= 8.)
        return 416.1015873015873 - 364.0888888888889 * x + 136.5333333333333 * x2 - 28.44444444444444 * x3 +
               3.555555555555556 * x4 - .2666666666666667 * x5 + .1111111111111111e-1 * x6 - .1984126984126984e-3 * x7;
    else
        return 0.;
} // kernel

// Elasticity model data.
namespace ModelData
{
// The tether penalty functions each require some data that is set in the
// input file. This data is passed to each object through the void *ctx
// context data pointer. Here we collect all relevant tether data in a struct:
struct TetherData
{
    const double c1_s;
    const double kappa_s_body;
    const double eta_s_body;
    const double kappa_s_surface;
    const double eta_s_surface;

    TetherData(Pointer<Database> input_db)
        : c1_s(input_db->getDouble("C1_S")),
          kappa_s_body(input_db->getDouble("KAPPA_S_BODY")),
          eta_s_body(input_db->getDouble("ETA_S_BODY")),
          kappa_s_surface(input_db->getDouble("KAPPA_S_SURFACE")),
          eta_s_surface(input_db->getDouble("ETA_S_SURFACE"))
    {
    }
};

bool use_boundary_mesh = false;
bool compute_fluid_traction = false;

// Tether (penalty) stress function.
void
PK1_stress_function(TensorValue<double>& PP,
                    const TensorValue<double>& FF,
                    const libMesh::Point& /*x*/,
                    const libMesh::Point& /*X*/,
                    Elem* const /*elem*/,
                    const vector<const vector<double>*>& /*var_data*/,
                    const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                    double /*time*/,
                    void* ctx)
{
    const TetherData* const tether_data = reinterpret_cast<TetherData*>(ctx);

    PP = 2.0 * tether_data->c1_s * (FF - tensor_inverse_transpose(FF, NDIM));
    return;
} // PK1_stress_function

// Tether (penalty) force functions.
void
tether_force_function(VectorValue<double>& F,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double /*time*/,
                      void* ctx)
{
    const TetherData* const tether_data = reinterpret_cast<TetherData*>(ctx);

    const std::vector<double>& U = *var_data[0];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = tether_data->kappa_s_body * (X(d) - x(d)) - tether_data->eta_s_body * U[d];
    }
    return;
} // tether_force_function

void
tether_force_function(VectorValue<double>& F,
                      const VectorValue<double>& n,
                      const VectorValue<double>& N,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      const unsigned short /*side*/,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double /*time*/,
                      void* ctx)
{
    const TetherData* const tether_data = reinterpret_cast<TetherData*>(ctx);

    const std::vector<double>& U = *var_data[0];
    double u_bndry_n = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d) u_bndry_n += n(d) * U[d];

    for (unsigned int d = 0; d < NDIM; ++d)
        F(d) = tether_data->kappa_s_surface * (X(d) - x(d)) - tether_data->eta_s_surface * u_bndry_n * n(d);
    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

// Function prototypes
static ofstream drag_F_stream, lift_F_stream, drag_TAU_stream, lift_TAU_stream, U_L1_norm_stream, U_L2_norm_stream,
    U_max_norm_stream;

void postprocess_data(Pointer<Database> input_db,
                      Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
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
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Setup user-defined kernel function.
        LEInteractor::s_kernel_fcn = &kernel;
        LEInteractor::s_kernel_fcn_stencil_size = 8;

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
        const string exodus_filename = app_initializer->getExodusIIFilename();

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

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = input_db->getDouble("R");
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
        solid_mesh.boundary_info->sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        use_boundary_mesh = input_db->getBoolWithDefault("USE_BOUNDARY_MESH", false);

        compute_fluid_traction = input_db->getBoolWithDefault("COMPUTE_FLUID_TRACTION", false);

        Mesh& mesh = use_boundary_mesh ? boundary_mesh : solid_mesh;

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input database
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
        Pointer<IBStrategy> ib_ops;
        if (use_boundary_mesh)
            ib_ops = new IBFESurfaceMethod(
                "IBFESurfaceMethod",
                app_initializer->getComponentDatabase("IBFESurfaceMethod"),
                &mesh,
                app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        else
            ib_ops =
                new IBFEMethod("IBFEMethod",
                               app_initializer->getComponentDatabase("IBFEMethod"),
                               &mesh,
                               app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_ops,
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
        TetherData tether_data(input_db);
        void* const tether_data_ptr = reinterpret_cast<void*>(&tether_data);
        EquationSystems* equation_systems;
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        if (use_boundary_mesh)
        {
            vector<SystemData> sys_data(1, SystemData(IBFESurfaceMethod::VELOCITY_SYSTEM_NAME, vars));
            Pointer<IBFESurfaceMethod> ibfe_ops = ib_ops;
            ibfe_ops->initializeFEEquationSystems();
            equation_systems = ibfe_ops->getFEDataManager()->getEquationSystems();
            IBFESurfaceMethod::LagSurfaceForceFcnData surface_fcn_data(
                tether_force_function, sys_data, tether_data_ptr);
            ibfe_ops->registerLagSurfaceForceFunction(surface_fcn_data);
        }
        else
        {
            vector<SystemData> sys_data(1, SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, vars));
            Pointer<IBFEMethod> ibfe_ops = ib_ops;
            ibfe_ops->initializeFEEquationSystems();
            equation_systems = ibfe_ops->getFEDataManager()->getEquationSystems();
            IBFEMethod::PK1StressFcnData PK1_stress_data(
                PK1_stress_function, std::vector<IBTK::SystemData>(), tether_data_ptr);
            PK1_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_QUAD_ORDER", "THIRD"));
            ibfe_ops->registerPK1StressFunction(PK1_stress_data);

            IBFEMethod::LagBodyForceFcnData body_fcn_data(tether_force_function, sys_data, tether_data_ptr);
            ibfe_ops->registerLagBodyForceFunction(body_fcn_data);

            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data(tether_force_function, sys_data, tether_data_ptr);
            ibfe_ops->registerLagSurfaceForceFunction(surface_fcn_data);

            if (input_db->getBoolWithDefault("ELIMINATE_PRESSURE_JUMPS", false))
            {
                ibfe_ops->registerStressNormalizationPart();
            }
        }

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
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        if (use_boundary_mesh)
        {
            Pointer<IBFESurfaceMethod> ibfe_ops = ib_ops;
            ibfe_ops->initializeFEData();
        }
        else
        {
            Pointer<IBFEMethod> ibfe_ops = ib_ops;
            ibfe_ops->initializeFEData();
        }
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
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save lift and drag coefficients and the norms of the
        // velocity.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_F_stream.open("C_F_D.curve", ios_base::out | ios_base::trunc);
            lift_F_stream.open("C_F_L.curve", ios_base::out | ios_base::trunc);
            drag_TAU_stream.open("C_T_D.curve", ios_base::out | ios_base::trunc);
            lift_TAU_stream.open("C_T_L.curve", ios_base::out | ios_base::trunc);
            U_L1_norm_stream.open("U_L1.curve", ios_base::out | ios_base::trunc);
            U_L2_norm_stream.open("U_L2.curve", ios_base::out | ios_base::trunc);
            U_max_norm_stream.open("U_max.curve", ios_base::out | ios_base::trunc);

            drag_F_stream.precision(10);
            lift_F_stream.precision(10);
            drag_TAU_stream.precision(10);
            lift_TAU_stream.precision(10);
            U_L1_norm_stream.precision(10);
            U_L2_norm_stream.precision(10);
            U_max_norm_stream.precision(10);
        }

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
                    exodus_io->write_timestep(
                        exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
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
                postprocess_data(input_db,
                                 patch_hierarchy,
                                 navier_stokes_integrator,
                                 mesh,
                                 equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_F_stream.close();
            lift_F_stream.close();
            drag_TAU_stream.close();
            lift_TAU_stream.close();
            U_L1_norm_stream.close();
            U_L2_norm_stream.close();
            U_max_norm_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return true;
} // main

void
postprocess_data(Pointer<Database> input_db,
                 Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
                 const int /*iteration_num*/,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{
    TetherData tether_data(input_db);
    void* const tether_data_ptr = reinterpret_cast<void*>(&tether_data);

    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    double T_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F_integral[d] = 0.0;
        T_integral[d] = 0.0;
    }
    System* x_system;
    System* U_system;

    if (use_boundary_mesh)
    {
        x_system = &equation_systems->get_system(IBFESurfaceMethod::COORDS_SYSTEM_NAME);
        U_system = &equation_systems->get_system(IBFESurfaceMethod::VELOCITY_SYSTEM_NAME);
    }
    else
    {
        x_system = &equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
        U_system = &equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    }
    NumericVector<double>* x_vec = x_system->solution.get();
    NumericVector<double>* x_ghost_vec = x_system->current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system->solution.get();
    NumericVector<double>* U_ghost_vec = U_system->current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system->get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    NumericVector<double>& X_vec = x_system->get_vector("INITIAL_COORDINATES");

    std::vector<std::vector<unsigned int> > WSS_o_dof_indices(NDIM);
    System* TAU_system;
    NumericVector<double>* TAU_ghost_vec = NULL;
    if (compute_fluid_traction)
    {
        TAU_system = &equation_systems->get_system(IBFESurfaceMethod::TAU_OUT_SYSTEM_NAME);

        TAU_ghost_vec = TAU_system->current_local_solution.get();
    }

    std::unique_ptr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<libMesh::Point>& q_point = fe->get_xyz();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();

    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule_face = QBase::build(QGAUSS, dim - 1, SEVENTH);
    fe_face->attach_quadrature_rule(qrule_face.get());
    const vector<double>& JxW_face = fe_face->get_JxW();
    const vector<libMesh::Point>& q_point_face = fe_face->get_xyz();
    const vector<libMesh::Point>& normal_face = fe_face->get_normals();
    const vector<vector<double> >& phi_face = fe_face->get_phi();
    const vector<vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;

    TensorValue<double> FF, FF_inv_trans;
    boost::multi_array<double, 2> x_node, X_node, U_node, TAU_node;

    VectorValue<double> F, N, U, n, x, X, TAU;

    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);
        get_values_for_interpolation(X_node, X_vec, dof_indices);
        if (compute_fluid_traction) get_values_for_interpolation(TAU_node, *TAU_ghost_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X, qp, X_node, phi);
            interpolate(x, qp, x_node, phi);
            jacobian(FF, qp, x_node, dphi);
            interpolate(U, qp, U_node, phi);
            if (compute_fluid_traction) interpolate(TAU, qp, TAU_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U(d);
            }
            if (use_boundary_mesh)
                tether_force_function(F, n, N, FF, x, X, elem, 0, var_data, grad_var_data, loop_time, tether_data_ptr);
            else
                tether_force_function(F, n, x, q_point[qp], elem, var_data, grad_var_data, loop_time, tether_data_ptr);

            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F(d) * JxW[qp];
                if (compute_fluid_traction) T_integral[d] += TAU(d) * JxW[qp];
            }
        }
        if (!use_boundary_mesh)
        {
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor_ptr(side)) continue;
                fe_face->reinit(elem, side);
                const unsigned int n_qp_face = qrule_face->n_points();
                for (unsigned int qp = 0; qp < n_qp_face; ++qp)
                {
                    interpolate(x, qp, x_node, phi_face);
                    jacobian(FF, qp, x_node, dphi_face);
                    interpolate(U, qp, U_node, phi_face);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        U_qp_vec[d] = U(d);
                    }
                    N = normal_face[qp];
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    n = (FF_inv_trans * N).unit();

                    tether_force_function(F,
                                          n,
                                          N,
                                          FF,
                                          x,
                                          q_point_face[qp],
                                          elem,
                                          side,
                                          var_data,
                                          grad_var_data,
                                          loop_time,
                                          tether_data_ptr);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_integral[d] += F(d) * JxW_face[qp];
                    }
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    SAMRAI_MPI::sumReduction(T_integral, NDIM);
    static const double rho = 1.0;
    static const double U_max = 1.0;
    static const double D = 1.0;
    if (SAMRAI_MPI::getRank() == 0)
    {
        drag_F_stream << loop_time << " " << -F_integral[0] / (0.5 * rho * U_max * U_max * D) << endl;
        lift_F_stream << loop_time << " " << -F_integral[1] / (0.5 * rho * U_max * U_max * D) << endl;
        if (compute_fluid_traction)
        {
            drag_TAU_stream << loop_time << " " << T_integral[0] / (0.5 * rho * U_max * U_max * D) << endl;
            lift_TAU_stream << loop_time << " " << T_integral[1] / (0.5 * rho * U_max * U_max * D) << endl;
        }
    }
    return;
} // postprocess_data
