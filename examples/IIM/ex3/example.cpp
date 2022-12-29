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

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStrategySet.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
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

VectorValue<double> x_com_initial;

System *Y_solid_system, *U_b_solid_system;
System* dY_solid_system;
System *x_new_airfoil_system, *u_new_airfoil_system;

// Tether (penalty) stress function.

// Tether (penalty) force functions.
static VectorValue<double> F_mean, T_mean;
static double R;

struct AirfoilParams
{
    double rho_s;
    double theta_m;
    double fr;
};

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
UpdateSurfaceMeshCoordinates(Mesh& boundary_mesh,
                             EquationSystems* airfoil_equation_systems,
                             EquationSystems* solid_systems)
{
    ReplicatedMesh* mesh = &boundary_mesh;
    EquationSystems* equation_systems = airfoil_equation_systems;

    MeshBase& mesh_solid = solid_systems->get_mesh();

    System& X_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    auto it = mesh->local_nodes_begin();
    const auto end_it = mesh->local_nodes_end();
    for (; it != end_it; ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            const libMesh::Point& X = *n;
            libMesh::Point x;

            const auto el_begin = mesh_solid.active_local_elements_begin();
            const auto el_end = mesh_solid.active_local_elements_end();
            for (auto el_it = el_begin; el_it != el_end; ++el_it)
            {
                const Elem* elem_c = *el_it;
                if (elem_c->contains_point(X))
                {
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x(d) = x_new_airfoil_system->point_value(d, X, elem_c);
                    }
                }
            }

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num, d, 0);
                X_coords.set(dof_index, x(d));
            }
        }
    }
    X_coords.close();
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    copy_and_synch(X_coords, *X_system.current_local_solution);
    copy_and_synch(X_coords, X_system.get_vector("INITIAL_COORDINATES"));
    return;
} // UpdateSurfaceMeshCoordinates

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
                      double /*current_time*/,
                      void* /*ctx*/)
{
    // tether_force_function() is called on elements of the boundary mesh.  Here
    // we look up the element in the solid mesh that the current boundary
    // element was extracted from.
    const Elem* const interior_parent = elem->interior_parent();

    VectorValue<double> X_solid;
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

    double ds = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ds += (x_solid[d] - x_bndry(d)) * (x_solid[d] - x_bndry(d));
    }
    ds = sqrt(ds);
    TBOX_ASSERT(ds < 1.5 * dx);
    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) =
            kappa_s * (x_solid[d] - x_bndry(d)) + eta_s * (u_solid[d] - u_bndry[d]); // (u_solid_n - u_bndry_n) * n(d);
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
        // get_values_for_interpolation(X0_node, X0_vec, dof_indices);

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
    SAMRAI_MPI::sumReduction(&F(0), NDIM + 1);
    SAMRAI_MPI::sumReduction(&T(0), NDIM + 1);

    x_ghost_vec->close();
    TAU_ghost_vec->close();

    return;

} // calculateFluidForceAndTorque

void
update_solid_vel_pos(VectorValue<double>& Y,
                     VectorValue<double>& U,
                     VectorValue<double> Y_n,
                     VectorValue<double> U_n,
                     VectorValue<double> sum_Y,
                     VectorValue<double> x_com,
                     const double vol,
                     EquationSystems* equation_bndry_systems,
                     EquationSystems& solid_equation_systems,
                     const double current_time,
                     const double dts,
                     void* ctx)
{
    // Ensure ghost values are communicated.

    AirfoilParams* params = static_cast<AirfoilParams*>(ctx);
    double rho_s = params->rho_s;
    double freq = params->fr;
    double theta_m = params->theta_m;

    const double m = vol * rho_s;

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

    VectorValue<double> F_s, Torque, Y_iter, U_iter, Omega;
    F_s.zero();
    Torque.zero();

    Y_iter = { 100, 100 };
    U_iter = U_n;

    calculateFluidForceAndTorque(F_s, Torque, x_com, new_time, equation_bndry_systems);

    U_iter = (F_s / m) * dts + U_n;
    Y_iter = dts * U_iter + Y_n;

    Y = Y_iter;
    U = U_iter;

    double T_iter = theta_m * sin(2 * M_PI * freq * current_time);
    double T_M = -theta_m * sin(2 * M_PI * freq * (current_time - dts));
    double W_iter = (2 * M_PI * theta_m) * cos(2 * M_PI * freq * current_time);
    const TensorValue<double> Q(cos(T_iter), -sin(T_iter), 0.0, sin(T_iter), cos(T_iter), 0.0, 0.0, 0.0, 1.0);
    const TensorValue<double> QM(cos(T_M), -sin(T_M), 0.0, sin(T_M), cos(T_M), 0.0, 0.0, 0.0, 1.0);
    const TensorValue<double> P(1, 0, sum_Y(0), 0, 1, sum_Y(1), 0, 0, 1);
    const TensorValue<double> PM(1, 0, -sum_Y(0), 0, 1, -sum_Y(1), 0, 0, 1);

    Omega.zero();
    Omega(2) = W_iter;

    const auto node_end = mesh.local_nodes_end();
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != node_end; ++it)
    {
        Node* n = *it;
        if (n->n_vars(dY_sys_num))
        {
            const libMesh::Point& X = *n;
            libMesh::Point S = X;
            S = P * (Q * (PM * X));
            VectorValue<double> U_b;
            U_b.zero();
            U_b = U_iter + Omega.cross(X - Y_n);
            S += Y_iter;

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
static ofstream position_stream, velocity_stream, CL_stream, CD_stream;

void postprocess_CL_CD(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
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

        Pointer<Database> airfoil_params_db = app_initializer->getComponentDatabase("AirfoilParams");

        AirfoilParams airfoil_params;
        airfoil_params.rho_s = airfoil_params_db->getDoubleWithDefault("RHOS", 1.0);
        airfoil_params.theta_m = airfoil_params_db->getDoubleWithDefault("THETA_M", 0.0);
        airfoil_params.fr = airfoil_params_db->getDoubleWithDefault("FREQUENCY", 0.0);

        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string viz_dump_dirname = app_initializer->getVizDumpDirectory();
        const string exodus_solid_filename = viz_dump_dirname + "/airfoil.ex2";
        const string exodus_bndry_filename = viz_dump_dirname + "/airfoil_bndry.ex2";

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
        solid_mesh.read(input_db->getString("MESH_FILENAME"));
        solid_mesh.get_boundary_info().clear_boundary_node_ids();

        string elem_type = input_db->getString("ELEM_TYPE");
        R = 0.05;

        using MeshTools::Modification::scale;
        scale(solid_mesh, 0.1, 0.1, 0.1);

        double theta = input_db->getDouble("INITIAL_THETA") * M_PI / 180.0;
        TensorValue<double> Rz(
            std::cos(theta), -std::sin(theta), 0.0, std::sin(theta), std::cos(theta), 0.0, 0.0, 0.0, 1.0);
        const auto node_end = solid_mesh.nodes_end();
        for (MeshBase::node_iterator it = solid_mesh.nodes_begin(); it != node_end; ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            const libMesh::Point s = X;
            X = Rz * s;
        }
        MeshRefinement mesh_refinement_solid(solid_mesh);
        mesh_refinement_solid.uniformly_refine(1);

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

        // Whether to use discontinuous basis functions with element-local support

        const bool USE_DISCON_ELEMS = input_db->getBool("USE_DISCON_ELEMS");

        if (USE_DISCON_ELEMS) ib_method_ops->registerDisconElemFamilyForJumps();

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
        const auto local_node_end = mesh.local_nodes_end();
        for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != local_node_end; ++it)
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
        for (auto it = mesh.local_nodes_begin(); it != local_node_end; ++it)
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

        // Open streams to save force coefficients

        if (SAMRAI_MPI::getRank() == 0)
        {
            position_stream.open("position_dt_" + std::to_string(double(DT)) + "_kappa_" +
                                     std::to_string(int(kappa_s)) + ".curve",
                                 ios_base::out | ios_base::trunc);
            velocity_stream.open("velocity_dt_" + std::to_string(double(DT)) + "_kappa_" +
                                     std::to_string(int(kappa_s)) + ".curve",
                                 ios_base::out | ios_base::trunc);
            CD_stream.open("CD_dt_" + std::to_string(double(DT)) + "_kappa_" + std::to_string(int(kappa_s)) + ".curve",
                           ios_base::out | ios_base::trunc);
            CL_stream.open("CL_dt_" + std::to_string(double(DT)) + "_kappa_" + std::to_string(int(kappa_s)) + ".curve",
                           ios_base::out | ios_base::trunc);

            position_stream.precision(10);
            velocity_stream.precision(10);
            CD_stream.precision(10);
            CL_stream.precision(10);
        }

        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        double vol = 0.0;
        VectorValue<double> x_com_current, U, Y, U_n, Y_n, sum_Y, Omega;
        x_com_current.zero();
        x_com_initial.zero();
        Y.zero();
        U.zero();
        Y_n.zero();
        U_n.zero();
        sum_Y.zero();
        double current_time;
        calculateGeomQuantitiesOfStructure(vol, x_com_initial, solid_equation_systems);

        Omega.zero();

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
            update_solid_vel_pos(Y,
                                 U,
                                 Y_n,
                                 U_n,
                                 sum_Y,
                                 x_com_current,
                                 vol,
                                 bndry_equation_systems,
                                 *solid_equation_systems,
                                 current_time,
                                 0.5 * dt,
                                 static_cast<void*>(&airfoil_params));
            Y_n = Y;
            U_n = U;
            sum_Y += Y;
            x_new_airfoil_system = &solid_equation_systems->get_system("RIGID_BODY_COORDS_SYSTEM");
            u_new_airfoil_system = &solid_equation_systems->get_system("RIGID_BODY_VELOCITY_SYSTEM");
            //~ UpdateSurfaceMeshCoordinates(bndry_mesh, bndry_equation_systems, solid_equation_systems);

            time_integrator->advanceHierarchy(dt);
            //~ UpdateSurfaceMeshCoordinates(bndry_mesh, bndry_equation_systems, solid_equation_systems);
            loop_time += dt;

            // Update the body configuration.
            update_solid_vel_pos(Y,
                                 U,
                                 Y_n,
                                 U_n,
                                 sum_Y,
                                 x_com_current,
                                 vol,
                                 bndry_equation_systems,
                                 *solid_equation_systems,
                                 current_time + 0.5 * dt,
                                 0.5 * dt,
                                 static_cast<void*>(&airfoil_params));
            Y_n = Y;
            U_n = U;
            sum_Y += Y;
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

            postprocess_CL_CD(patch_hierarchy,
                              navier_stokes_integrator,
                              bndry_mesh,
                              bndry_equation_systems,
                              iteration_num,
                              loop_time,
                              postproc_data_dump_dirname);

            if (SAMRAI_MPI::getRank() == 0 && (iteration_num % 20 == 0))
            {
                position_stream << loop_time << "\t" << Y(0) << "\t" << Y(1) << endl;
                velocity_stream << loop_time << "\t" << U(0) << "\t" << U(1) << endl;
            }
        }
        Y_solid_system->clear();
        dY_solid_system->clear();
        U_b_solid_system->clear();

        solid_equation_systems->clear();
        bndry_equation_systems->clear();

        if (SAMRAI_MPI::getRank() == 0)
        {
            position_stream.close();
            velocity_stream.close();
            CL_stream.close();
            CD_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main

void
postprocess_CL_CD(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                  Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                  Mesh& mesh,
                  EquationSystems* equation_systems,
                  const int /*iteration_num*/,
                  const double loop_time,
                  const string& /*data_dump_dirname*/)
{
    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    double T_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d) T_integral[d] = 0.0;

    System& x_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IIMethod::VELOCITY_SYSTEM_NAME);
    System& TAU_system = equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);

    NumericVector<double>* TAU_vec = TAU_system.solution.get();
    NumericVector<double>* TAU_ghost_vec = TAU_system.current_local_solution.get();
    TAU_vec->localize(*TAU_ghost_vec);
    DofMap& TAU_dof_map = TAU_system.get_dof_map();
    std::vector<std::vector<unsigned int> > TAU_dof_indices(NDIM);
    std::unique_ptr<FEBase> fe(FEBase::build(dim, TAU_dof_map.variable_type(0)));

    NumericVector<double>& X0_vec = x_system.get_vector("INITIAL_COORDINATES");

    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system.solution.get();
    NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF_qp;
    boost::multi_array<double, 2> x_node, U_node, TAU_node, X0_node;
    VectorValue<double> F_qp, U_qp, x_qp, W_qp, TAU_qp, N, n, X;

    const MeshBase::element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
            TAU_dof_map.dof_indices(elem, TAU_dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(TAU_node, *TAU_ghost_vec, TAU_dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);
        get_values_for_interpolation(X0_node, X0_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X, qp, X0_node, phi);
            interpolate(x_qp, qp, x_node, phi);
            jacobian(FF_qp, qp, x_node, dphi);
            interpolate(U_qp, qp, U_node, phi);
            interpolate(TAU_qp, qp, TAU_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U_qp(d);
            }
            tether_force_function(
                F_qp, n, N, FF_qp, x_qp, X, elem, 0, var_data, grad_var_data, loop_time, force_fcn_ctx);
            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F_qp(d) * JxW[qp];
                T_integral[d] += TAU_qp(d) * JxW[qp];
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    SAMRAI_MPI::sumReduction(T_integral, NDIM);

    if (SAMRAI_MPI::getRank() == 0)
    {
        CD_stream << loop_time << " " << -F_integral[0] << endl;
        CL_stream << loop_time << " " << -F_integral[1] << endl;
    }

    return;
} // postprocess_CL_CD
