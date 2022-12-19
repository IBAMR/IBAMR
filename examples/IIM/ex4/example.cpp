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
// --------------------------------------------------------------------

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/centroid_partitioner.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/parallel.h>
#include <libmesh/quadrature.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Application-specific includes.

// Elasticity model data.
namespace ModelData
{
static const double TOL = sqrt(std::numeric_limits<double>::epsilon());
static double dx = 0.0;

System *Y_solid_system, *U_b_solid_system;
System* dY_solid_system;

// Tether (penalty) stress function.

// Tether (penalty) force functions.
static VectorValue<double> F_mean, T_mean;
static double R;

static double kappa_s = 0.0;
static double eta_s = 0.0;
static double DT = 0.0;

void
calculateGeomQuantitiesOfStructure(double& vol,                             // volume of the body
                                   VectorValue<double>& x_com,              // new center of the mass
                                   EquationSystems* solid_equation_systems) // mass density of the body (assumed to be
                                                                            // uniform)
{
    // Get the structure mesh for codim-0 solid.
    // For now the eqs are setup only for one part but this will be extended
    // to multiple parts

    MeshBase& mesh = solid_equation_systems->get_mesh();
    System& X_system = solid_equation_systems->get_system("RIGID_BODY_COORDS_SYSTEM");

    const unsigned int dim = mesh.mesh_dimension();

    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);

    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType fe_type = X_dof_map.variable_type(0);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    PetscVector<double>& X_petsc = dynamic_cast<PetscVector<double>&>(*X_system.current_local_solution.get());
    X_petsc.close();
    Vec X_global_vec = X_petsc.vec();
    Vec X_local_ghost_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
    double* X_local_ghost_soln;
    VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);

    // 3D identity tensor.

    vol = 0.0;
    x_com.zero();

    boost::multi_array<double, 2> X_node;

    // double X_qp_new[NDIM], X_qp_current[NDIM], R_qp_current[NDIM], R_qp_new[NDIM];
    VectorValue<double> X_qp, R_qp;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi);

            x_com += X_qp * JxW[qp];

            vol += JxW[qp];
        }
    }
    SAMRAI_MPI::sumReduction(&x_com(0), NDIM);
    SAMRAI_MPI::sumReduction(&vol);

    x_com /= vol;

    VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);
    X_system.solution->close();

    return;

} // calculateGeomQuantitiesOfStructure

void
tether_force_function(VectorValue<double>& F,
                      const VectorValue<double>& /*n*/,
                      const VectorValue<double>& /*N*/,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the boundary mesh
                      const libMesh::Point& X_bndry, // X_bndry gives reference coordinates on the boundary mesh
                      Elem* const elem,
                      const unsigned short /*side*/,
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
        x_solid[d] = Y_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = U_b_solid_system->point_value(d, X_bndry, interior_parent);
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

void
calculateFluidForceAndTorque(VectorValue<double>& F, // net force  acting on the body
                             VectorValue<double>& T, // net torque acting on the body
                             VectorValue<double> x_com_current,
                             const double /*loop_time*/,
                             EquationSystems* equation_systems)
{
    MeshBase& mesh = equation_systems->get_mesh();

    const unsigned int dim = mesh.mesh_dimension();
    F.zero();
    T.zero();

    System& x_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    System& TAU_system = equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);

    NumericVector<double>* TAU_vec = TAU_system.solution.get();
    NumericVector<double>* TAU_ghost_vec = TAU_system.current_local_solution.get();
    TAU_vec->localize(*TAU_ghost_vec);
    DofMap& TAU_dof_map = TAU_system.get_dof_map();
    std::vector<std::vector<unsigned int> > TAU_dof_indices(NDIM);
    std::unique_ptr<FEBase> fe(FEBase::build(dim, TAU_dof_map.variable_type(0)));

    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double> >& phi = fe->get_phi();

    boost::multi_array<double, 2> x_node, TAU_node;
    VectorValue<double> F_qp, x_qp, W_qp, TAU_qp, R_qp;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
            TAU_dof_map.dof_indices(elem, TAU_dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(TAU_node, *TAU_ghost_vec, TAU_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            interpolate(TAU_qp, qp, TAU_node, phi);

            R_qp = x_qp - x_com_current;

            F += TAU_qp * JxW[qp];
            T += R_qp.cross(TAU_qp) * JxW[qp];
        }
    }
    SAMRAI_MPI::sumReduction(&F(0), NDIM);
    SAMRAI_MPI::sumReduction(&T(0), NDIM);

    x_ghost_vec->close();
    TAU_ghost_vec->close();

    return;

} // calculateFluidForceAndTorque

void
update_solid_vel_pos(const double* params,
                     VectorValue<double>& Y,
                     VectorValue<double>& U,
                     VectorValue<double> x_com,
                     EquationSystems* equation_bndry_systems,
                     EquationSystems& solid_equation_systems,
                     const double current_time,
                     const double dts)
{
    // Ensure ghost values are communicated.

    const double& m = params[0];
    const double& k = params[1];
    const double& c = params[2];
    const double new_time = current_time + dts;

    MeshBase& mesh = solid_equation_systems.get_mesh();
    System& Y_system = solid_equation_systems.get_system("RIGID_BODY_COORDS_SYSTEM");
    const unsigned int Y_sys_num = Y_system.number();
    System& U_b_system = solid_equation_systems.get_system("RIGID_BODY_VELOCITY_SYSTEM");
    const unsigned int U_b_sys_num = U_b_system.number();
    System& dY_system = solid_equation_systems.get_system("RIGID_BODY_COORDS_MAPPING_SYSTEM");
    const unsigned int dY_sys_num = dY_system.number();

    NumericVector<double>& Y_coords = *Y_system.solution;
    NumericVector<double>& U_b_coords = *U_b_system.solution;
    NumericVector<double>& dY_coords = *dY_system.solution;

    VectorValue<double> F_s, Torque, Y_n, U_n, U_iter, Y_iter;
    F_s.zero();
    Torque.zero();
    Y_iter = { 100, 100 };
    Y_n = Y;
    U_n = U;
    U_iter = U;

    calculateFluidForceAndTorque(F_s, Torque, x_com, new_time, equation_bndry_systems);

    int iter = 0;
    while ((Y - Y_iter).norm() > TOL)
    {
        // m d2y/dt2 + c*dy/dt + k*y = F_y.
        U = U_iter;
        Y = Y_iter;

        Y_iter = 0.5 * dts * (U_n + U_iter) + Y_n;

        U_iter = ((F_s / m - k * Y_iter / m) * dts + U_n) / (1.0 + c * dts / m);

        iter += 1;
    }

    pout << " Number of time integration iterations = " << iter << "\n\n";

    const auto node_end = mesh.local_nodes_end();
    for (auto it = mesh.local_nodes_begin(); it != node_end; ++it)
    {
        Node* n = *it;
        if (n->n_vars(dY_sys_num))
        {
            const libMesh::Point& X = *n;
            libMesh::Point S = X;
            S = X + Y;
            VectorValue<double> U_b;
            U_b.zero();
            U_b = U;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int Y_dof_index = n->dof_number(Y_sys_num, d, 0);
                Y_coords.set(Y_dof_index, S(d));
                const int U_b_dof_index = n->dof_number(U_b_sys_num, d, 0);
                U_b_coords.set(U_b_dof_index, U_b(d));
                const int dY_dof_index = n->dof_number(dY_sys_num, d, 0);
                dY_coords.set(dY_dof_index, S(d) - X(d));
            }
        }
    }
    Y_coords.close();
    Y_system.get_dof_map().enforce_constraints_exactly(Y_system, &Y_coords);
    Y_system.solution->localize(*Y_system.current_local_solution);

    dY_coords.close();
    dY_system.get_dof_map().enforce_constraints_exactly(dY_system, &dY_coords);
    dY_system.solution->localize(*dY_system.current_local_solution);

    U_b_coords.close();
    U_b_system.get_dof_map().enforce_constraints_exactly(U_b_system, &U_b_coords);
    U_b_system.solution->localize(*U_b_system.current_local_solution);

    return;
} // update_solid_vel_pos

} // namespace ModelData
using namespace ModelData;

// Function prototypes
static ofstream position_stream, l_max_stream;

void postprocess_Convergence(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
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

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-10");
    PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-10");

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        dx = input_db->getDouble("DX");
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string viz_dump_dirname = app_initializer->getVizDumpDirectory();
        const string exodus_solid_filename = viz_dump_dirname + "/disk.ex2";
        const string exodus_bndry_filename = viz_dump_dirname + "/disk_bndry.ex2";

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

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);

        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        R = 0.5;
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
            //  triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = false;
            triangle.smooth_after_generating() = false;
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

        BoundaryMesh bndry_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        solid_mesh.get_boundary_info().sync(bndry_mesh);
        bndry_mesh.prepare_for_use();

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
        DT = input_db->getDouble("DT");

        const bool second_order_mesh = (elem_type != "TRI3" && elem_type != "QUAD4");
        const Order fe_order = second_order_mesh ? SECOND : FIRST;
        const FEFamily fe_family = LAGRANGE;

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
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid option is: STAGGERED");
        }
        Pointer<IIMethod> ib_method_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         &bndry_mesh,
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

        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        ib_method_ops->initializeFEEquationSystems();

        EquationSystems* bndry_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));

        libMesh::EquationSystems* solid_equation_systems(new EquationSystems(solid_mesh));
        Y_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("RIGID_BODY_COORDS_SYSTEM");
        dY_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("RIGID_BODY_COORDS_MAPPING_SYSTEM");
        U_b_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("RIGID_BODY_VELOCITY_SYSTEM");

        //~ Order fe_order = FIRST;
        //~ FEFamily fe_family = LAGRANGE;
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "Y_" << d;
            Y_solid_system->add_variable(os.str(), fe_order, fe_family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "dY_" << d;
            dY_solid_system->add_variable(os.str(), fe_order, fe_family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_b_" << d;
            U_b_solid_system->add_variable(os.str(), fe_order, fe_family);
        }

        solid_equation_systems->init();

        IIMethod::LagSurfaceForceFcnData tether_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(tether_fcn_data);

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
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_solid_io(uses_exodus ? new ExodusII_IO(solid_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_io(uses_exodus ? new ExodusII_IO(bndry_mesh) : NULL);

        // Initialize FE data.
        ib_method_ops->initializeFEData();

        MeshBase& mesh = solid_equation_systems->get_mesh();
        System& Y_init_solid_system = solid_equation_systems->get_system("RIGID_BODY_COORDS_SYSTEM");
        const unsigned int Y_init_solid_sys_num = Y_init_solid_system.number();
        NumericVector<double>& Y_init_solid_coords = *Y_init_solid_system.solution;
        const auto node_end = mesh.local_nodes_end();
        for (auto it = mesh.local_nodes_begin(); it != node_end; ++it)
        {
            Node* n = *it;
            if (n->n_vars(Y_init_solid_sys_num))
            {
                TBOX_ASSERT(n->n_vars(Y_init_solid_sys_num) == NDIM);
                const libMesh::Point& Y = *n;
                libMesh::Point y = Y;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const int dof_index = n->dof_number(Y_init_solid_sys_num, d, 0);
                    Y_init_solid_coords.set(dof_index, y(d));
                }
            }
        }
        Y_init_solid_coords.close();
        Y_init_solid_system.get_dof_map().enforce_constraints_exactly(Y_init_solid_system, &Y_init_solid_coords);
        Y_init_solid_system.solution->localize(*Y_init_solid_system.current_local_solution);

        System& dY_init_solid_system = solid_equation_systems->get_system("RIGID_BODY_COORDS_MAPPING_SYSTEM");
        const unsigned int dY_init_solid_sys_num = dY_init_solid_system.number();
        NumericVector<double>& dY_init_solid_coords = *dY_init_solid_system.solution;
        for (auto it = mesh.local_nodes_begin(); it != node_end; ++it)
        {
            Node* n = *it;
            if (n->n_vars(dY_init_solid_sys_num))
            {
                TBOX_ASSERT(n->n_vars(dY_init_solid_sys_num) == NDIM);
                const libMesh::Point& Y = *n;

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const int dY_dof_index = n->dof_number(dY_init_solid_sys_num, d, 0);
                    const int Y_dof_index = n->dof_number(Y_init_solid_sys_num, d, 0);
                    dY_init_solid_coords.set(dY_dof_index, Y_init_solid_coords(Y_dof_index) - Y(d));
                }
            }
        }
        dY_init_solid_coords.close();
        dY_init_solid_system.get_dof_map().enforce_constraints_exactly(dY_init_solid_system, &dY_init_solid_coords);
        dY_init_solid_system.solution->localize(*dY_init_solid_system.current_local_solution);

        System& U_b_init_solid_system = solid_equation_systems->get_system("RIGID_BODY_VELOCITY_SYSTEM");
        NumericVector<double>& U_b_vec = *U_b_init_solid_system.solution;
        U_b_vec.zero();
        U_b_vec.close();
        U_b_vec.localize(*U_b_init_solid_system.current_local_solution);

        Y_solid_system->assemble_before_solve = false;
        Y_solid_system->assemble();

        dY_solid_system->assemble_before_solve = false;
        dY_solid_system->assemble();

        U_b_solid_system->assemble_before_solve = false;
        U_b_solid_system->assemble();

        pout << "\nRigid body variables setup done.\n";

        // Initialize hierarchy configuration and data on all patches.
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

        // Open streams to save lift and drag coefficients and the norms of the
        // velocity.

        // State variables for the rigid body.
        const double m = input_db->getDouble("BODY_MASS");
        const double k = input_db->getDouble("BODY_SPRING_CONSTANT");
        const double c = input_db->getDouble("BODY_DAMPING_CONSTANT");
        double params[3] = { m, k, c };

        if (SAMRAI_MPI::getRank() == 0)
        {
            position_stream.open("Displacement_2DOF_dt_" + std::to_string(double(DT)) + "_kappa_" +
                                     std::to_string(int(kappa_s)) + ".curve",
                                 ios_base::out | ios_base::trunc);
            l_max_stream.open("Lmax_2DOF_dt_" + std::to_string(double(DT)) + "_kappa_" + std::to_string(int(kappa_s)) +
                                  ".curve",
                              ios_base::out | ios_base::trunc);

            position_stream.precision(10);

            l_max_stream.precision(10);
        }

        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        double vol = 0.0;
        VectorValue<double> x_com_current, U, Y;
        x_com_current.zero();
        Y.zero();
        U.zero();
        double current_time;

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

            calculateGeomQuantitiesOfStructure(vol, x_com_current, solid_equation_systems);
            // Predict the body position.
            update_solid_vel_pos(
                params, Y, U, x_com_current, bndry_equation_systems, *solid_equation_systems, current_time, 0.5 * dt);

            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // Update the body configuration.
            update_solid_vel_pos(params,
                                 Y,
                                 U,
                                 x_com_current,
                                 bndry_equation_systems,
                                 *solid_equation_systems,
                                 current_time + 0.5 * dt,
                                 0.5 * dt);

            pout << "cylinder position = " << Y << "\n";
            pout << "cylinder velocity = " << U << "\n";

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
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            postprocess_Convergence(patch_hierarchy,
                                    navier_stokes_integrator,
                                    bndry_mesh,
                                    bndry_equation_systems,
                                    iteration_num,
                                    loop_time,
                                    postproc_data_dump_dirname);

            if (SAMRAI_MPI::getRank() == 0 && (iteration_num % 20 == 0))
            {
                position_stream << loop_time << "\t" << Y(0) << "\t" << Y(1) << "\t" << U(0) << "\t" << U(1) << endl;
            }
        }
        Y_solid_system->clear();
        dY_solid_system->clear();
        U_b_solid_system->clear();

        solid_equation_systems->clear();
        bndry_equation_systems->clear();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main

void
postprocess_Convergence(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                        Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                        Mesh& mesh,
                        EquationSystems* equation_systems,
                        const int /*iteration_num*/,
                        const double loop_time,
                        const string& /*data_dump_dirname*/)
{
    const unsigned int dim = mesh.mesh_dimension();
    System& x_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IIMethod::VELOCITY_SYSTEM_NAME);

    NumericVector<double>& X0_vec = x_system.get_vector("INITIAL_COORDINATES");

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
    const vector<vector<double> >& phi = fe->get_phi();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;

    TensorValue<double> FF_qp;
    boost::multi_array<double, 2> x_node, U_node, TAU_node, X_node;
    VectorValue<double> F_qp, U_qp, x_qp, X_qp, W_qp, TAU_qp, N, n, X;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();

    DofMap& U_dof_map = U_system.get_dof_map();
    std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);

    VectorValue<double> tau1, tau2;
    double Lmax_diff = 0.0;

    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* elem = *el_it;
        fe->reinit(elem);
        const int n_qp = qrule->n_points();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
            U_dof_map.dof_indices(elem, U_dof_indices[d], d);
        }

        get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(X_node, X0_vec, dof_indices);

        const Elem* const interior_par = elem->interior_parent();

        for (int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            interpolate(X_qp, qp, X_node, phi);

            std::vector<double> x_solid(NDIM, 0.0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_solid[d] = Y_solid_system->point_value(d, X_qp, interior_par);
            }

            double disx = 0.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                disx += (x_solid[d] - x_qp(d)) * (x_solid[d] - x_qp(d));
            }

            disx = sqrt(disx);

            Lmax_diff = std::max(Lmax_diff, disx);
        }
    }

    SAMRAI_MPI::maxReduction(&Lmax_diff, 1);

    if (SAMRAI_MPI::getRank() == 0)
    {
        l_max_stream << loop_time << " " << Lmax_diff << endl;
    }

    return;
} // postprocess_Convergence
