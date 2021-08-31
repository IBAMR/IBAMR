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

#include "InterpolationUtilities.h"

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double rho = 1.0;
static double dx = -1.0;
static bool ERROR_ON_MOVE = false;
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
                      double /*time*/,
                      void* /*ctx*/)
{
    const std::vector<double>& U = *var_data[0];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s * (X(d) - x(d)) - eta_s * U[d];
    }
    std::vector<double> d(NDIM);
    d[0] = std::abs(x(0) - X(0));
    d[1] = std::abs(X(1) - x(1));
    if (ERROR_ON_MOVE && ((d[0] > 0.25 * dx) || (d[1] > 0.25 * dx)))
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
void
setInsideOfCylinder(const int d_idx, Pointer<PatchHierarchy<NDIM> > patch_hierarchy, const double time, const double R);

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
        dx = input_db->getDouble("DX");
        ERROR_ON_MOVE = input_db->getBool("ERROR_ON_MOVE");
        const double mfac = input_db->getDouble("MFAC");
        const double ds = mfac * dx;
        std::cout << dx << " " << mfac << "\n";
        std::cout << ds << "\n";
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 1.0;
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

        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        Mesh& mesh = boundary_mesh;

        // MeshRefinement mesh_refinement (mesh);
        // mesh_refinement.uniformly_refine (3);
        // MeshTools::Modification::flatten (mesh);

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
        rho = input_db->getDouble("RHO");

        mesh.print_info();

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
                                  &mesh,
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
        IBFESurfaceMethod::LagSurfaceForceFcnData body_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(body_fcn_data);
        EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

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

        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

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
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
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
            if (dump_viz_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 adv_diff_integrator,
                                 polymericStressForcing,
                                 mesh,
                                 equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
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

void
postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                 Pointer<CFINSForcing> polymericStressForcing,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname)
{
    pout << "Outputting data at time: " << loop_time << "\n And iteration num : " << iteration_num << "\n";
    {
        // Output files
        string file_name = data_dump_dirname + "/hier_data.";
        char temp_buf[128];
        sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
        file_name += temp_buf;
        Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
        hier_db->create(file_name);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        ComponentSelector hier_data;
        if (polymericStressForcing)
        {
            hier_data.setFlag(var_db->mapVariableAndContextToIndex(polymericStressForcing->getVariable(),
                                                                   adv_diff_integrator->getCurrentContext()));
            pout << "current ctx " << adv_diff_integrator->getCurrentContext()->getName() << "\n";
        }
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                               navier_stokes_integrator->getCurrentContext()));
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                               navier_stokes_integrator->getCurrentContext()));
        patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
        hier_db->putDouble("loop_time", loop_time);
        hier_db->putInteger("iteration_num", iteration_num);
        hier_db->close();
    }
    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    double F_int = 0.0;

    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;

    System& x_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);

    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system.solution.get();
    NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<libMesh::Point>& q_point = fe->get_xyz();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();

    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule_face(QBase::build(QGAUSS, dim - 1, SEVENTH));
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
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF, FF_inv_trans;
    boost::multi_array<double, 2> x_node, U_node;
    VectorValue<double> F, N, U, n, x;

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

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x, qp, x_node, phi);
            jacobian(FF, qp, x_node, dphi);
            interpolate(U, qp, U_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U(d);
            }
            tether_force_function(
                F, n, N, FF, x, q_point[qp], elem, 0, var_data, grad_var_data, loop_time, force_fcn_ctx);
            F_int += F(0) * JxW[qp];
            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F(d) * JxW[qp];
            }
        }

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
                tether_force_function(
                    F, n, N, FF, x, q_point_face[qp], elem, side, var_data, grad_var_data, loop_time, force_fcn_ctx);
                for (int d = 0; d < NDIM; ++d)
                {
                    F_integral[d] += F(d) * JxW_face[qp];
                }
            }
        }
    }
    IBTK_MPI::sumReduction(F_integral, NDIM);

    if (IBTK_MPI::getRank() == 0)
    {
        drag_force_stream << loop_time << " " << -F_integral[0] << endl;
    }

    // Interpolate sxx value along cylinder surface
    if (polymericStressForcing)
    {
        if (IBTK_MPI::getRank() == 0) sxx_component_stream << loop_time << " ";
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
        Pointer<Patch<NDIM> > patch = level->getPatch(PatchLevel<NDIM>::Iterator(level)());
        const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
        const double* dx = p_geom->getDx();
        double xp, yp, sxx;
        std::vector<double> X(2);
        int num_pts = static_cast<int>(3.0 / dx[0]);
        for (int i = 0; i < num_pts; ++i)
        {
            xp = -4.0 + i * dx[0];
            yp = 0.0;
            X[0] = xp;
            X[1] = yp;
            sxx = polymericStressForcing->getViscosity() / polymericStressForcing->getRelaxationTime() *
                  (IBTK::InterpolationUtilities::interpolate(
                       X,
                       polymericStressForcing->getVariableIdx(),
                       polymericStressForcing->getVariable(),
                       patch_hierarchy,
                       adv_diff_integrator->getPhysicalBcCoefs(polymericStressForcing->getVariable()),
                       1) -
                   1.0);
            if (IBTK_MPI::getRank() == 0) sxx_component_stream << xp << " " << yp << " " << sxx << " ";
        }
        double dr = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1]);
        num_pts = static_cast<int>(M_PI / dr);
        for (int i = 0; i < num_pts; ++i)
        {
            xp = -std::cos(M_PI * i / num_pts);
            yp = std::sin(M_PI * i / num_pts);
            X[0] = xp;
            X[1] = yp;
            sxx = polymericStressForcing->getViscosity() / polymericStressForcing->getRelaxationTime() *
                  (IBTK::InterpolationUtilities::interpolate(
                       X,
                       polymericStressForcing->getVariableIdx(),
                       polymericStressForcing->getVariable(),
                       patch_hierarchy,
                       adv_diff_integrator->getPhysicalBcCoefs(polymericStressForcing->getVariable()),
                       1) -
                   1.0);
            if (IBTK_MPI::getRank() == 0) sxx_component_stream << xp << " " << yp << " " << sxx << " ";
        }
        num_pts = static_cast<int>(3.0 / dx[0]);
        for (int i = 0; i < num_pts; ++i)
        {
            xp = 1.0 + i * dx[0];
            yp = 0.0;
            X[0] = xp;
            X[1] = yp;
            sxx = polymericStressForcing->getViscosity() / polymericStressForcing->getRelaxationTime() *
                  (IBTK::InterpolationUtilities::interpolate(
                       X,
                       polymericStressForcing->getVariableIdx(),
                       polymericStressForcing->getVariable(),
                       patch_hierarchy,
                       adv_diff_integrator->getPhysicalBcCoefs(polymericStressForcing->getVariable()),
                       1) -
                   1.0);
            if (IBTK_MPI::getRank() == 0) sxx_component_stream << xp << " " << yp << " " << sxx << " ";
        }
        if (IBTK_MPI::getRank() == 0) sxx_component_stream << endl;
    }
    return;
} // postprocess_data
