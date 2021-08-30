// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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
#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_subdivision_support.h"
#include "libmesh/mesh_tools.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include "ibamr/CFINSForcing.h"
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBFESurfaceMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

struct DiskData
{
public:
    DiskData(std::vector<double> r, double w, double dx) : r_c(r), omega_s(w), d_dx(dx)
    {
        return;
    }
    std::vector<double> r_c;
    double omega_s, d_dx;
};

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double omega_s_0 = 0.0;
static double omega_s_1 = 0.0;
static double R_s = 1.0;
void
tether_force_function(VectorValue<double>& F,
                      const libMesh::VectorValue<double>& /*n*/,
                      const libMesh::VectorValue<double>& /*N*/,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      unsigned short int /*side*/,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double time,
                      void* ctx)
{
    DiskData data = *(static_cast<DiskData*>(ctx));
    const double& w = data.omega_s;
    const std::vector<double>& r_c = data.r_c;
    const double dx = data.d_dx;

    const std::vector<double>& U = *var_data[0];
    double rx = X(0) - r_c[0];
    double ry = X(1) - r_c[1];
    const double th = std::atan2(ry, rx);
    F(0) =
        kappa_s * (R_s * std::cos(th + w * time) + r_c[0] - x(0)) - eta_s * (-R_s * w * std::sin(th + w * time) - U[0]);
    F(1) =
        kappa_s * (R_s * std::sin(th + w * time) + r_c[1] - x(1)) - eta_s * (R_s * w * std::cos(th + w * time) - U[1]);
    std::vector<double> d(NDIM);
    d[0] = std::abs(R_s * std::cos(th + w * time) + r_c[0] - x(0));
    d[1] = std::abs(R_s * std::sin(th + w * time) + r_c[1] - x(1));
    if ((d[0] > 0.2 * dx) || (d[1] > 0.2 * dx))
    {
        TBOX_ERROR("Structure has moved too much.\n");
    }
    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

// Function prototypes
static ofstream drag_force_stream, sxx_component_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                      Pointer<CFINSForcing> polymericStressForcing,
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
        const string exodus_filename_0 = app_initializer->getExodusIIFilename("TopLeft");
        const string exodus_filename_1 = app_initializer->getExodusIIFilename("TopRight");
        const string exodus_filename_2 = app_initializer->getExodusIIFilename("BotRight");
        const string exodus_filename_3 = app_initializer->getExodusIIFilename("BotLeft");

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
        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
        omega_s_0 = input_db->getDouble("OMEGA_S");
        double R = input_db->getDouble("RADIUS");
        double B = input_db->getDouble("B");

        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double mfac = input_db->getDouble("MFAC");
        const double ds = mfac * dx;
        std::cout << dx << " " << mfac << "\n";
        std::cout << ds << "\n";
        string elem_type = input_db->getString("ELEM_TYPE");
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
                    Node& n = *elem->node_ptr(k);
                    n = R * n.unit() /*+libMesh::Point(0.25,0.25)*/;
                }
            }
        }
        solid_mesh.prepare_for_use();

        BoundaryMesh roller_tl(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(roller_tl);

        BoundaryMesh roller_tr = roller_tl;
        BoundaryMesh roller_br = roller_tl;
        BoundaryMesh roller_bl = roller_tl;

        MeshTools::Modification::translate(roller_tl, -B, B);
        MeshTools::Modification::translate(roller_tr, B, B);
        MeshTools::Modification::translate(roller_br, B, -B);
        MeshTools::Modification::translate(roller_bl, -B, -B);

        roller_tl.prepare_for_use();
        roller_tr.prepare_for_use();
        roller_br.prepare_for_use();
        roller_bl.prepare_for_use();

        omega_s_1 = -1.0 * omega_s_0;
        std::vector<double> tl_r(NDIM), tr_r(NDIM), br_r(NDIM), bl_r(NDIM);
        tl_r[0] = -B;
        tl_r[1] = B;
        tr_r[0] = B;
        tr_r[1] = B;
        br_r[0] = B;
        br_r[1] = -B;
        bl_r[0] = -B;
        bl_r[1] = -B;
        DiskData tl(tl_r, omega_s_0, dx), tr(tr_r, omega_s_1, dx), br(br_r, omega_s_0, dx), bl(bl_r, omega_s_1, dx);

        std::vector<MeshBase*> meshes;
        meshes.push_back(&roller_tl);
        meshes.push_back(&roller_tr);
        meshes.push_back(&roller_br);
        meshes.push_back(&roller_bl);

        // MeshTools::Subdivision::prepare_subdivision_mesh (mesh, false);

        // Print information about the subdivision mesh to the screen.
        // mesh.print_info();

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
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        Pointer<IBFESurfaceMethod> ib_method_ops =
            new IBFESurfaceMethod("IBFEMethod",
                                  app_initializer->getComponentDatabase("IBFEMethod"),
                                  meshes,
                                  app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
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
        vector<SystemData> sys_data(1, SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, vars));

        IBFESurfaceMethod::LagSurfaceForceFcnData body_fcn_data_tl(tether_force_function, sys_data, &tl);
        IBFESurfaceMethod::LagSurfaceForceFcnData body_fcn_data_tr(tether_force_function, sys_data, &tr);
        IBFESurfaceMethod::LagSurfaceForceFcnData body_fcn_data_br(tether_force_function, sys_data, &br);
        IBFESurfaceMethod::LagSurfaceForceFcnData body_fcn_data_bl(tether_force_function, sys_data, &bl);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data_tl, 0);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data_tr, 1);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data_br, 2);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data_bl, 3);
        EquationSystems* equation_systems_0 = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* equation_systems_1 = ib_method_ops->getFEDataManager(1)->getEquationSystems();
        EquationSystems* equation_systems_2 = ib_method_ops->getFEDataManager(2)->getEquationSystems();
        EquationSystems* equation_systems_3 = ib_method_ops->getFEDataManager(3)->getEquationSystems();

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

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create Eulerian body force function specification objects.
        Pointer<CFINSForcing> polymericStressForcing;
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }
        if (input_db->keyExists("ComplexFluid"))
        {
            polymericStressForcing = new CFINSForcing("PolymericStressForcing",
                                                      app_initializer->getComponentDatabase("ComplexFluid"),
                                                      navier_stokes_integrator,
                                                      grid_geometry,
                                                      adv_diff_integrator,
                                                      visit_data_writer);
            time_integrator->registerBodyForceFunction(polymericStressForcing);
        }

        std::unique_ptr<ExodusII_IO> exodus_io_0(uses_exodus ? new ExodusII_IO(*(meshes[0])) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_io_1(uses_exodus ? new ExodusII_IO(*(meshes[1])) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_io_2(uses_exodus ? new ExodusII_IO(*(meshes[2])) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_io_3(uses_exodus ? new ExodusII_IO(*(meshes[3])) : NULL);

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
                exodus_io_0->write_timestep(
                    exodus_filename_0, *equation_systems_0, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_io_1->write_timestep(
                    exodus_filename_1, *equation_systems_1, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_io_2->write_timestep(
                    exodus_filename_2, *equation_systems_2, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_io_3->write_timestep(
                    exodus_filename_3, *equation_systems_3, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save lift and drag coefficients and the norms of the
        // velocity.
        if (IBTK_MPI::getRank() == 0)
        {
            drag_force_stream.open("CD_Force_integral.curve", ios_base::out | ios_base::trunc);
            drag_force_stream.precision(10);
            sxx_component_stream.open("Sxx_component.curve", ios_base::out | ios_base::trunc);
            sxx_component_stream.precision(10);
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
                    exodus_io_0->write_timestep(
                        exodus_filename_0, *equation_systems_0, iteration_num / viz_dump_interval + 1, loop_time);
                    exodus_io_1->write_timestep(
                        exodus_filename_1, *equation_systems_1, iteration_num / viz_dump_interval + 1, loop_time);
                    exodus_io_2->write_timestep(
                        exodus_filename_2, *equation_systems_2, iteration_num / viz_dump_interval + 1, loop_time);
                    exodus_io_3->write_timestep(
                        exodus_filename_3, *equation_systems_3, iteration_num / viz_dump_interval + 1, loop_time);
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
        }

        // Close the logging streams.
        if (IBTK_MPI::getRank() == 0)
        {
            drag_force_stream.close();
            sxx_component_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main
