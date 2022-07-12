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
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IIMethod.h>
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

// Elasticity model data.
namespace ModelData
{
// The tether penalty functions each require some data that is set in the
// input file. This data is passed to each object through the void *ctx
// context data pointer. Here we collect all relevant tether data in a struct:
struct TetherData
{
    const double kappa_s_surface;
    const double eta_s_surface;

    TetherData(tbox::Pointer<tbox::Database> input_db)
        : kappa_s_surface(input_db->getDouble("KAPPA_S_SURFACE")), eta_s_surface(input_db->getDouble("ETA_S_SURFACE"))
    {
    }
};
static double ds = 0.0;
static double dx = 0.0;

bool compute_fluid_traction = false;

void
tether_force_function(VectorValue<double>& F,
                      const VectorValue<double>& n,
                      const VectorValue<double>& /*N*/,
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

void postprocess_data(tbox::Pointer<tbox::Database> input_db,
                      tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
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

    // Set a tight tolerance usign with IIM for the ksp sovler
    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-10");
    PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-10");

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        tbox::Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        tbox::Pointer<tbox::Database> input_db = app_initializer->getInputDatabase();

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
        const double mfac = input_db->getDouble("MFAC");
        ds = mfac * dx;
        const double R_D = input_db->getDouble("R_D");
        string elem_type = input_db->getString("ELEM_TYPE");
        if (NDIM == 2)
        {
            MeshTools::Generation::build_square(solid_mesh,
                                                static_cast<int>(ceil(R_D / ds)),
                                                static_cast<int>(ceil(R_D / ds)),
                                                -R_D / 2.0,
                                                R_D / 2.0,
                                                -R_D / 2.0,
                                                R_D / 2.0,
                                                Utility::string_to_enum<ElemType>(elem_type));
        }
        else
        {
            MeshTools::Generation::build_cube(solid_mesh,
                                              static_cast<int>(ceil(R_D / ds)),
                                              static_cast<int>(ceil(R_D / ds)),
                                              static_cast<int>(ceil(R_D / ds)),
                                              R_D / 2.0,
                                              -R_D / 2.0,
                                              -R_D / 2.0,
                                              R_D / 2.0,
                                              -R_D / 2.0,
                                              R_D / 2.0,
                                              Utility::string_to_enum<ElemType>(elem_type));
        }
        solid_mesh.prepare_for_use();

        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        compute_fluid_traction = input_db->getBoolWithDefault("COMPUTE_FLUID_TRACTION", false);

        Mesh& mesh = boundary_mesh;

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
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
        tbox::Pointer<IBStrategy> ib_ops;
        ib_ops = new IIMethod("IIMethod",
                              app_initializer->getComponentDatabase("IIMethod"),
                              &mesh,
                              app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        tbox::Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_ops,
                                              navier_stokes_integrator);
        tbox::Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        tbox::Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        tbox::Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        tbox::Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        tbox::Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
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

        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));
        tbox::Pointer<IIMethod> ibfe_ops = ib_ops;

        // Whether to use discontinuous basis functions with element-local support
        // We ask this before initializing the FE equation system
        const bool USE_DISCON_ELEMS = input_db->getBool("USE_DISCON_ELEMS");
        if (USE_DISCON_ELEMS) ibfe_ops->registerDisconElemFamilyForJumps();

        ibfe_ops->initializeFEEquationSystems();
        equation_systems = ibfe_ops->getFEDataManager()->getEquationSystems();
        IIMethod::LagSurfaceForceFcnData surface_fcn_data(tether_force_function, sys_data, tether_data_ptr);
        ibfe_ops->registerLagSurfaceForceFunction(surface_fcn_data);

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            tbox::Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            tbox::Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
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
            tbox::Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        tbox::Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ibfe_ops->initializeFEData();
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
} // main

void
postprocess_data(tbox::Pointer<tbox::Database> input_db,
                 tbox::Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 tbox::Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
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

    x_system = &equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    U_system = &equation_systems->get_system(IIMethod::VELOCITY_SYSTEM_NAME);
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
        TAU_system = &equation_systems->get_system(IIMethod::TAU_OUT_SYSTEM_NAME);

        TAU_ghost_vec = TAU_system->current_local_solution.get();
    }

    std::unique_ptr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();

    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule_face = QBase::build(QGAUSS, dim - 1, SEVENTH);
    fe_face->attach_quadrature_rule(qrule_face.get());

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;

    TensorValue<double> FF, FF_inv_trans;
    boost::multi_array<double, 2> x_node, X_node, U_node, TAU_node;

    VectorValue<double> F, N, U, n, x, X, TAU;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const auto elem = *el_it;
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
            tether_force_function(F, n, N, FF, x, X, elem, 0, var_data, grad_var_data, loop_time, tether_data_ptr);

            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F(d) * JxW[qp];
                if (compute_fluid_traction) T_integral[d] += TAU(d) * JxW[qp];
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
