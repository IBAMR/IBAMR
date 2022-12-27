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

#include <ibamr/FEMechanicsExplicitIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStrategySet.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/MathUtilities.h>
#include <tbox/Utilities.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_submatrix.h>
#include <libmesh/dense_subvector.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_solver_package.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/fe.h>
#include <libmesh/getpot.h>
#include <libmesh/gnuplot_io.h>
#include <libmesh/libmesh.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/perf_log.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/petsc_macro.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/solver_configuration.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/zero_function.h>

#include <boost/multi_array.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
static BoundaryInfo* lag_bdry_info;

// Tether (penalty) force function for the solid blocks.
static double kappa_s = 1.0e6;
static double eta_s = 1.0e6;
static double DX = 0.01;
static double DT = 0.001;

static bool use_volumetric_term;
static std::string stress_funtion;

static ofstream l_max_stream;

System* x_new_solid_system;
System* x_new_surface_system;
System* Tau_new_surface_system;
System* u_new_solid_system;

EquationSystems* boundary_systems;

void
solid_surface_force_function(VectorValue<double>& F,
                             const VectorValue<double>& /*n*/,
                             const VectorValue<double>& /*N*/,
                             const TensorValue<double>& /*FF*/,
                             const libMesh::Point& /*x*/,
                             const libMesh::Point& X,
                             Elem* const /*elem*/,
                             const unsigned short int /*side*/,
                             const vector<const vector<double>*>& /*var_data*/,
                             const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    MeshBase& mesh_bndry = boundary_systems->get_mesh();
    std::vector<double> x_surface(NDIM, 0.0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const auto el_begin = mesh_bndry.active_local_elements_begin();
        const auto el_end = mesh_bndry.active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* elem_bndry = *el_it;
            if (elem_bndry->contains_point(X))
                F(d) = Tau_new_surface_system->point_value(d, X, elem_bndry); //&side_elem);
        }
    }

    return;
}

// Stress tensor functions.
static double c1_s = 0.05;
static double p0_s = 0.0;
static double beta_s = 0.0;
static double shear_mod = 0.0;
static double pr = 0.0;
static double bulk_mod = 0.0;

void
calculateGeomQuantitiesOfStructure(double& vol,                // mass of the body
                                   VectorValue<double>& x_com, // new center of the mass
                                   const FEMechanicsExplicitIntegrator* const fem_solver,
                                   EquationSystems* equation_systems)
{
    // Get the structure mesh for codim-0 solid.
    // For now the eqs are setup only for one part but this will be extended
    // to multiple parts

    // Extract the FE system and DOF map, and setup the FE object.

    System& X_system = equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
    X_system.solution->localize(*X_system.current_local_solution);
    MeshBase& mesh = equation_systems->get_mesh();

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
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    x_com.zero();
    vol = 0.0; // volume of the body
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
    SAMRAI_MPI::sumReduction(&x_com(0), NDIM + 1);
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

    const std::vector<double>& u_bndry = *var_data[0];

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.
    std::vector<double> x_solid(NDIM, 0.0), u_solid(NDIM, 0.0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_solid_system->point_value(d, X_bndry, interior_parent);
    }
    double dx_length = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx_length += (x_solid[d] - x_bndry(d)) * (x_solid[d] - x_bndry(d));
    }
    dx_length = sqrt(dx_length);
    TBOX_ASSERT(dx_length < 0.5 * DX);

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s * (x_solid[d] - x_bndry(d)) +
               eta_s * (u_solid[d] - u_bndry[d]); //(u_cylinder_n - u_bndry_n) * n(d);
    }
    return;
} // tether_force_function

void
PK1_dev_stress_function_disk(TensorValue<double>& PP,
                             const TensorValue<double>& FF,
                             const libMesh::Point& /*X*/,
                             const libMesh::Point& /*s*/,
                             Elem* const /*elem*/,
                             const vector<const vector<double>*>& /*var_data*/,
                             const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    static const TensorValue<double> IO(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    const TensorValue<double> CC = FF.transpose() * FF;
    const double J = FF.det();
    const double I1 = (FF.transpose() * FF).tr();
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    PP = shear_mod * pow(J, -2.0 / 3.0) * (FF - (I1 / 3.0) * FF_inv_trans);
    return;
} // PK1_dev_stress_function

void
PK1_dil_stress_function_disk(TensorValue<double>& PP,
                             const TensorValue<double>& FF,
                             const libMesh::Point& /*X*/,
                             const libMesh::Point& /*s*/,
                             Elem* const /*elem*/,
                             const vector<const vector<double>*>& /*var_data*/,
                             const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    const double J = FF.det();
    PP = bulk_mod * J * log(J) * FF_inv_trans;
    return;
} // PK1_dil_stress_function_disk

Real
initial_jacobian(const libMesh::Point& /*p*/,
                 const Parameters& /*parameters*/,
                 const string& /*system*/,
                 const string& var_name)
{
    // could include something more complicated
    if (var_name == "Avg J")
        return 1.0;
    else
        return 0.0;
}

void
apply_initial_jacobian(EquationSystems& es, const string& system_name)
{
    libmesh_assert_equal_to(system_name, "JacobianDeterminant");
    ExplicitSystem& system = es.get_system<ExplicitSystem>("JacobianDeterminant");
    es.parameters.set<Real>("time") = system.time = 0;
    system.project_solution(initial_jacobian, NULL, es.parameters);
}

} // namespace ModelData
using namespace ModelData;

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
    SAMRAIManager::startup();

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
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string viz_dump_dirname = app_initializer->getVizDumpDirectory();
        const string exodus_filename = viz_dump_dirname + "/disk.ex2";
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
        ReplicatedMesh mesh(init.comm(), NDIM);
        DX = input_db->getDouble("DX");
        double MAX_LEVELS = input_db->getDouble("MAX_LEVELS");
        DT = input_db->getDouble("DT");
        const double ds = input_db->getDouble("MFAC") * DX;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 0.2;
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(mesh);
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
            MeshTools::Generation::build_sphere(mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        // Ensure nodes on the surface are on the analytic boundary.
        MeshBase::element_iterator el_end = mesh.elements_end();
        for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
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
        const auto node_end = mesh.nodes_end();
        for (MeshBase::node_iterator it = mesh.nodes_begin(); it != node_end; ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(0) += 0.6;
            X(1) += 0.5;
        }

        mesh.prepare_for_use();

        BoundaryMesh boundary_mesh(mesh.comm(), mesh.mesh_dimension() - 1);
        mesh.get_boundary_info().sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        c1_s = input_db->getDouble("C1_S");
        pr = input_db->getDouble("POISSON_RATIO");
        shear_mod = input_db->getDouble("SHEAR_MOD");
        bulk_mod = input_db->getDouble("BULK_MOD");
        p0_s = input_db->getDouble("P0_S");
        beta_s = input_db->getDouble("BETA_S");
        eta_s = input_db->getDouble("ETA_S");
        lag_bdry_info = &mesh.get_boundary_info();

        // Setup the model parameters.
        kappa_s = input_db->getDouble("KAPPA_S");

        use_volumetric_term = input_db->getBool("USE_VOLUMETRIC_TERM");
        stress_funtion = input_db->getString("STRESS_FUNCTION");

        // Setup the time stepping parameters.
        const double loop_time_end = input_db->getDouble("END_TIME");
        double dt = input_db->getDouble("DT");

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

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<FEMechanicsExplicitIntegrator> fem_solver = new FEMechanicsExplicitIntegrator(
            "FEMechanicsExplicitIntegrator",
            app_initializer->getComponentDatabase("FEMechanicsExplicitIntegrator"),
            &mesh,
            app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        Pointer<IIMethod> ibfe_bndry_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         &boundary_mesh,
                         app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        vector<Pointer<IBStrategy> > ib_method_ops(2);
        ib_method_ops[0] = fem_solver;
        ib_method_ops[1] = ibfe_bndry_ops;

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ibfe_bndry_ops,
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

        VectorValue<double> x_com;
        x_com.zero();
        double vol = 0.0;

        // attach velocity

        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> velocity_data(1);
        velocity_data[0] = SystemData(fem_solver->getVelocitySystemName(), vars);

        ibfe_bndry_ops->initializeFEEquationSystems();

        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));

        // Configure the FE solver.
        FEMechanicsBase::LagSurfaceForceFcnData solid_surface_force_data(solid_surface_force_function, velocity_data);
        fem_solver->registerLagSurfaceForceFunction(solid_surface_force_data);

        FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_disk, velocity_data);
        PK1_dev_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        fem_solver->registerPK1StressFunction(PK1_dev_stress_data);

        IIMethod::LagSurfaceForceFcnData surface_fcn_data(tether_force_function, sys_data);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_fcn_data);

        EquationSystems* bndry_equation_systems = ibfe_bndry_ops->getFEDataManager()->getEquationSystems();

        FEMechanicsBase::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function_disk);
        PK1_dil_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
        fem_solver->registerPK1StressFunction(PK1_dil_stress_data);

        fem_solver->initializeFEEquationSystems();
        EquationSystems* equation_systems = fem_solver->getEquationSystems();
        ExplicitSystem& jac_system = equation_systems->add_system<ExplicitSystem>("JacobianDeterminant");
        jac_system.attach_init_function(apply_initial_jacobian);

        // **************Create Eulerian initial condition specification objects.**************** //
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

        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_io(uses_exodus ? new ExodusII_IO(boundary_mesh) : NULL);

        ibfe_bndry_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Initialize solver data.
        fem_solver->initializeFEData();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = 0;
        double loop_time = 0.0;
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

                exodus_bndry_io->write_timestep(
                    exodus_bndry_filename, *bndry_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save volume of structure.
        ofstream volume_stream, x_com_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            volume_stream.open("volume_MAX_LEVELS_" + std::to_string(double(MAX_LEVELS)) + "_kappa_" +
                                   std::to_string(int(kappa_s)) + "_eta_" + std::to_string(double(eta_s)) + "_DT_" +
                                   std::to_string(double(DT)) + ".curve",
                               ios_base::out | ios_base::trunc);
            x_com_stream.open("x_com_MAX_LEVELS_" + std::to_string(double(MAX_LEVELS)) + "_kappa_" +
                                  std::to_string(int(kappa_s)) + "_eta_" + std::to_string(double(eta_s)) + "_DT_" +
                                  std::to_string(double(DT)) + ".curve",
                              ios_base::out | ios_base::trunc);
            l_max_stream.open("Lmax_MAX_LEVELS_" + std::to_string(double(MAX_LEVELS)) + "_kappa_" +
                                  std::to_string(int(kappa_s)) + "_eta_" + std::to_string(double(eta_s)) + "_DT_" +
                                  std::to_string(double(DT)) + ".curve",
                              ios_base::out | ios_base::trunc);

            x_com_stream.precision(10);
            volume_stream.precision(10);
            l_max_stream.precision(10);
        }

        calculateGeomQuantitiesOfStructure(vol, x_com, fem_solver, equation_systems);
        const double n_cycles = input_db->getDouble("NCYCLE");
        // Main time step loop.
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end))
        {
            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            boundary_systems = bndry_equation_systems;
            // to compute average J for each element

            System& X_system = equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
            x_new_solid_system = &X_system;
            u_new_solid_system = &equation_systems->get_system<System>(fem_solver->getVelocitySystemName());
            Tau_new_surface_system = &bndry_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_surface_system = &bndry_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            dt = time_integrator->getMaximumTimeStepSize();

            for (int ii = 0; ii < static_cast<int>(n_cycles); ii++)
            {
                fem_solver->preprocessIntegrateData(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                    loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles,
                                                    /*num_cycles*/ 1);
                fem_solver->modifiedTrapezoidalStep(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                    loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles);
                fem_solver->postprocessIntegrateData(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                     loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles,
                                                     /*num_cycles*/ 1);
            }

            x_new_solid_system = &equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
            u_new_solid_system = &equation_systems->get_system<System>(fem_solver->getVelocitySystemName());

            time_integrator->advanceHierarchy(dt);
            boundary_systems = bndry_equation_systems;
            Tau_new_surface_system = &bndry_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_surface_system = &bndry_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            for (int ii = 0; ii < static_cast<int>(n_cycles); ii++)
            {
                fem_solver->preprocessIntegrateData(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                    loop_time +
                                                        (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles,
                                                    /*num_cycles*/ 1);
                fem_solver->modifiedTrapezoidalStep(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                    loop_time +
                                                        (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles);
                fem_solver->postprocessIntegrateData(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                     loop_time +
                                                         (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles,
                                                     /*num_cycles*/ 1);
            }

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
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0))
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

                    exodus_bndry_io->write_timestep(exodus_bndry_filename,
                                                    *bndry_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            double J_integral = 0.0;

            NumericVector<double>* X_vec = X_system.solution.get();
            NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
            copy_and_synch(*X_vec, *X_ghost_vec);
            DofMap& X_dof_map = X_system.get_dof_map();
            vector<vector<unsigned int> > X_dof_indices(NDIM);
            std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
            std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
            fe->attach_quadrature_rule(qrule.get());
            const vector<double>& JxW = fe->get_JxW();
            const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();
            TensorValue<double> FF;
            boost::multi_array<double, 2> X_node;
            const auto el_begin = mesh.active_local_elements_begin();
            const auto el_end = mesh.active_local_elements_end();
            for (auto el_it = el_begin; el_it != el_end; ++el_it)
            {
                const Elem* elem = *el_it;
                fe->reinit(elem);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_dof_map.dof_indices(elem, X_dof_indices[d], d);
                }
                const int n_qp = qrule->n_points();
                get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    jacobian(FF, qp, X_node, dphi);
                    J_integral += abs(FF.det()) * JxW[qp];
                }
            }
            J_integral = SAMRAI_MPI::sumReduction(J_integral);

            calculateGeomQuantitiesOfStructure(vol, x_com, fem_solver, equation_systems);

            postprocess_Convergence(patch_hierarchy,
                                    navier_stokes_integrator,
                                    boundary_mesh,
                                    bndry_equation_systems,
                                    iteration_num,
                                    loop_time,
                                    postproc_data_dump_dirname);

            if (SAMRAI_MPI::getRank() == 0)
            {
                volume_stream.precision(12);
                volume_stream.setf(ios::fixed, ios::floatfield);
                volume_stream << loop_time << " " << J_integral << endl;
                x_com_stream << loop_time << "\t " << x_com(0) << "\t " << x_com(1) << "\t " << vol << endl;
            }

        } // end of time loop

        if (SAMRAI_MPI::getRank() == 0)
        {
            x_com_stream.close();
            volume_stream.close();
        }
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
                x_solid[d] = x_new_solid_system->point_value(d, X_qp, interior_par);
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
