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
#include <libmesh/dense_vector.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/fe.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>

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
static double kappa_s_FSI_block = 1.0e6;
static double kappa_s_FSI_beam = 1.0e6;
static double eta_FSI_beam = 1.0e6;
static double kappa_s_block = 1.0e6;
static double kappa_s_surface_block = 1e6;
static double DX = 0.01;
static double D = 1.0;
static double mu_s, lambda_s, bulk_mod;
static std::string stress_funtion;

System* x_new_solid_system;
System* u_new_solid_system;
System* x_new_surface_system;
System* Tau_new_surface_system;

EquationSystems* boundary_systems;

static BoundaryInfo* beam_boundary_info;

// Stress tensor functions.
static double c1_s = 0.05;
static double p0_s = 0.0;

void
block_tether_force_function(VectorValue<double>& F,
                            const TensorValue<double>& /*FF*/,
                            const libMesh::Point& x,
                            const libMesh::Point& X,
                            Elem* const elem,
                            const vector<const vector<double>*>& /*var_data*/,
                            const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double /*time*/,
                            void* /*ctx*/)
{
    const libMesh::Point cp_elem = elem->centroid();

    if (cp_elem(0) < 0.25)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F(d) = kappa_s_block * (X(d) - x(d));
        }
    }
    else
    {
        F.zero();
    }

    return;
} // block_tether_force_function

void
FSI_tether_force_function(VectorValue<double>& F,
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
    const libMesh::Point cp_elem = elem->centroid();

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.
    std::vector<double> x_solid(NDIM, 0.0);
    std::vector<double> u_solid(NDIM, 0.0);

    const std::vector<double>& U = *var_data[0];

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_solid_system->point_value(d, X_bndry, interior_parent);
    }

    // The tether force is proportional to the mismatch between the positions
    // and velocities.

    double disp = 0.0;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        disp += (x_solid[d] - x_bndry(d)) * (x_solid[d] - x_bndry(d));
    }
    disp = sqrt(disp);
    TBOX_ASSERT(disp < DX);

    if (cp_elem(0) < 0.25)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F(d) = kappa_s_FSI_block * (x_solid[d] - x_bndry(d));
        }
    }
    else
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F(d) = kappa_s_FSI_beam * (x_solid[d] - x_bndry(d)) + eta_FSI_beam * (u_solid[d] - U[d]);
        }
    }
    return;
} // FSI_tether_force_function

void
solid_surface_force_function(VectorValue<double>& F,
                             const VectorValue<double>& /*n*/,
                             const VectorValue<double>& /*N*/,
                             const TensorValue<double>& /*FF*/,
                             const libMesh::Point& x,
                             const libMesh::Point& X,
                             Elem* const elem,
                             const unsigned short int side,
                             const vector<const vector<double>*>& /*var_data*/,
                             const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    MeshBase& mesh_bndry = boundary_systems->get_mesh();
    std::unique_ptr<Elem> side_elem = elem->side_ptr(side);

    const libMesh::Point cp_elem = elem->centroid();

    if (cp_elem(0) < 0.25)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F(d) = kappa_s_surface_block * (X(d) - x(d));
        }
    }
    else
    {
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
    }
    return;
} // solid_surface_force_function

void
PK1_dev_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const elem,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    static const TensorValue<double> IO(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);

    std::vector<double> x_surface(NDIM, 0.0);

    //  Unmodified St. Venant-Kirchhoff Model
#if 0
    const TensorValue<double> EE = 0.5 * (CC - II);
    PP = FF * (2 * mu_s * EE + lambda_s * EE.tr() * II);
#endif

//  Modified St. Venant-Kirchhoff Model
#if 0
    const TensorValue<double> CC_sq = CC * CC;
    const double J_c = CC.det();
    const double J_c_1_3 = pow(J_c, -1.0 / 3.0);
    const double I1 = CC.tr();
    const double I2 = 0.5 * (I1*I1 - CC_sq.tr());
    const double I1_bar = J_c_1_3 * I1;
    const double I2_bar = J_c_1_3 * J_c_1_3 * I2;
    const TensorValue<double> CC_inv_trans = tensor_inverse_transpose(CC, NDIM);
    const TensorValue<double> dW_dCC =
        (lambda_s / 4.0) * (I1_bar - 3.0) * (J_c_1_3 * II - (I1_bar / 3.0) * CC_inv_trans) +
        (mu_s / 2.0) *
            (J_c_1_3 * (J_c_1_3 * CC - II) + (1.0 / 3.0) * (I1_bar * (1.0 - I1_bar) + 2.0 * I2_bar) * CC_inv_trans);
    PP = 2.0 * FF * dW_dCC;
#endif

#if 0
     PP = 2.0 * c1_s * FF;

#endif

    const libMesh::Point cp_elem = elem->centroid();

    if (cp_elem(0) < 0.25) // if ( r < 0.05 )
    {
        const TensorValue<double> EE = 0.5 * (CC - II);
        PP = 2.0 * c1_s * (FF - tensor_inverse_transpose(FF, NDIM)); // FF * (2 * mu_s * EE + lambda_s * EE.tr() * II);
    }
    else
    {
        //  Unmodified St. Venant-Kirchhoff Model

        //  Modified Neo-Hookean Model
#if 1
        const double J = FF.det();
        const double I1 = (FF.transpose() * FF).tr();
        PP = mu_s * pow(J, -2.0 / 3.0) * (FF - (I1 / 3.0) * FF_inv_trans);

#endif
    }

    return;
} // PK1_dev_stress_function

void
beam_PK1_dil_stress_function(TensorValue<double>& PP,
                             const TensorValue<double>& FF,
                             const libMesh::Point& /*X*/,
                             const libMesh::Point& /*s*/,
                             Elem* const elem,
                             const std::vector<const std::vector<double>*>& /*var_data*/,
                             const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    double J = FF.det();
    TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);

    const libMesh::Point cp_elem = elem->centroid();

    if (cp_elem(0) < 0.25)
    {
        PP = 0.0;
    }
    else
    {
        PP = bulk_mod * J * log(J) * FF_inv_trans;
    }

    return;

} // beam_PK1_dil_stress_function

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

static ofstream drag_stream, lift_stream, A_x_posn_stream, A_y_posn_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      const FEMechanicsExplicitIntegrator* fem_solver,
                      ReplicatedMesh& beam_mesh,
                      EquationSystems* beam_equation_systems,
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
        ReplicatedMesh beam_mesh(init.comm(), NDIM);
        DX = input_db->getDouble("DX");
        D = input_db->getDouble("D");
        const double n_cycles = input_db->getDouble("NCYCLE");
        string beam_elem_type = input_db->getString("BEAM_ELEM_TYPE");
        beam_mesh.read(input_db->getString("BEAM_MESH_FILENAME"), NULL);

        const auto node_end = beam_mesh.nodes_end();
        for (MeshBase::node_iterator n_it = beam_mesh.nodes_begin(); n_it != node_end; ++n_it)
        {
            Node& n = **n_it;

            n(1) += 2 * D;
            n(0) += 2 * D;
        }

        MeshRefinement mesh_refinement_beam1(beam_mesh);
        mesh_refinement_beam1.uniformly_refine(1);

        beam_boundary_info = &beam_mesh.get_boundary_info();
        beam_mesh.prepare_for_use();
        BoundaryMesh boundary_mesh(beam_mesh.comm(), beam_mesh.mesh_dimension() - 1);
        beam_mesh.get_boundary_info().sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        mu_s = input_db->getDouble("MU_S");
        bulk_mod = input_db->getDouble("BULK_MOD");
        lambda_s = input_db->getDouble("LAMBDA_S");
        c1_s = input_db->getDouble("C1_S");
        p0_s = input_db->getDouble("P0_S");
        double pr = input_db->getDouble("POISSON_RATIO");

        lag_bdry_info = &beam_mesh.get_boundary_info();

        // Setup the model parameters.
        kappa_s_FSI_block = input_db->getDouble("KAPPA_S_FSI_BLOCK");
        kappa_s_FSI_beam = input_db->getDouble("KAPPA_S_FSI_BEAM");
        eta_FSI_beam = input_db->getDouble("ETA_FSI_BEAM");
        kappa_s_block = input_db->getDouble("KAPPA_S_BLOCK");
        kappa_s_surface_block = input_db->getDouble("KAPPA_S_SURFACE_BLOCK");

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
            &beam_mesh,
            app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        Pointer<IIMethod> ibfe_bndry_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         &boundary_mesh,
                         app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        vector<Pointer<IBStrategy> > ib_method_ops(2);
        ib_method_ops[0] = fem_solver;
        ib_method_ops[1] = ibfe_bndry_ops;
        Pointer<IBStrategySet> ib_method_set = new IBStrategySet(ib_method_ops.begin(), ib_method_ops.end());

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

        // attach velocity

        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> velocity_data(1);
        velocity_data[0] = SystemData(fem_solver->getVelocitySystemName(), vars);

        const bool USE_DISCON_ELEMS = input_db->getBool("USE_DISCON_ELEMS");
        const bool USE_NORMALIZED_PRESSURE_JUMP = input_db->getBool("USE_NORMALIZED_PRESSURE_JUMP");

        if (USE_DISCON_ELEMS) ibfe_bndry_ops->registerDisconElemFamilyForJumps();
        if (USE_NORMALIZED_PRESSURE_JUMP) ibfe_bndry_ops->registerPressureJumpNormalization();

        ibfe_bndry_ops->initializeFEEquationSystems();

        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));

        // Configure the FE solver.

        FEMechanicsBase::LagBodyForceFcnData block_tether_force_data(block_tether_force_function, velocity_data);
        fem_solver->registerLagBodyForceFunction(block_tether_force_data);

        FEMechanicsBase::LagSurfaceForceFcnData solid_surface_force_data(solid_surface_force_function, velocity_data);
        fem_solver->registerLagSurfaceForceFunction(solid_surface_force_data);

        FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function, velocity_data);
        PK1_dev_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        fem_solver->registerPK1StressFunction(PK1_dev_stress_data);

        FEMechanicsBase::PK1StressFcnData beam_PK1_dil_stress_data(beam_PK1_dil_stress_function, velocity_data);
        beam_PK1_dil_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "CONSTANT"));
        fem_solver->registerPK1StressFunction(beam_PK1_dil_stress_data);

        IIMethod::LagSurfaceForceFcnData surface_FSI_fcn_data(FSI_tether_force_function, sys_data);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_fcn_data);
        EquationSystems* bndry_equation_systems = ibfe_bndry_ops->getFEDataManager()->getEquationSystems();

        fem_solver->initializeFEEquationSystems();
        EquationSystems* equation_systems = fem_solver->getEquationSystems();

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

        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(beam_mesh) : NULL);
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

        // Open streams to save lift and drag coefficients.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.open("C_D_N3_MFAC6_" + std::to_string(int(10000 * pr)) + "_NCYCLE.curve",
                             ios_base::out | ios_base::trunc);
            lift_stream.open("C_L_N3_MFAC6_" + std::to_string(int(10000 * pr)) + "_NCYCLE.curve",
                             ios_base::out | ios_base::trunc);
            A_x_posn_stream.open("A_x_N3_MFAC6_" + std::to_string(int(10000 * pr)) + "_NCYCLE.curve",
                                 ios_base::out | ios_base::trunc);
            A_y_posn_stream.open("A_y_N3_MFAC6_" + std::to_string(int(10000 * pr)) + "_NCYCLE.curve",
                                 ios_base::out | ios_base::trunc);
        }

        // Open streams to save volume of structure.

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

            System& U_system = equation_systems->get_system<System>(fem_solver->getVelocitySystemName());
            u_new_solid_system = &U_system;

            Tau_new_surface_system = &bndry_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);

            x_new_surface_system = &bndry_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            dt = time_integrator->getMaximumTimeStepSize();

            for (int ii = 0; ii < static_cast<int>(n_cycles); ii++)
            {
                fem_solver->preprocessIntegrateData(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                    loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles,
                                                    /*num_cycles*/ 1);
                fem_solver->backwardEulerStep(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
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
                fem_solver->backwardEulerStep(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                              loop_time + (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles);
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
            const bool last_step = !time_integrator->stepsRemaining();
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

            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting state data...\n\n";
                postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 fem_solver,
                                 beam_mesh,
                                 equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }

        } // end of time loop

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
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

    SAMRAIManager::shutdown();
    return 0;
} // main

void
postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                 const FEMechanicsExplicitIntegrator* const fem_solver,
                 ReplicatedMesh& beam_mesh,
                 EquationSystems* beam_equation_systems,
                 const int /*iteration_num*/,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{
    double F_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
    ReplicatedMesh* mesh_solid = &beam_mesh;
    EquationSystems* equation_systems = beam_equation_systems;

    System& F_system = equation_systems->get_system<System>(fem_solver->getForceSystemName());
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
    const auto el_begin = mesh_solid->active_local_elements_begin();
    const auto el_end = mesh_solid->active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* elem = *el_it;
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

    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    if (SAMRAI_MPI::getRank() == 0)
    {
        drag_stream.precision(12);
        drag_stream.setf(ios::fixed, ios::floatfield);
        drag_stream << loop_time << " " << -F_integral[0] << endl;
        lift_stream.precision(12);
        lift_stream.setf(ios::fixed, ios::floatfield);
        lift_stream << loop_time << " " << -F_integral[1] << endl;
    }

    System& X_system = beam_equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
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
    if (SAMRAI_MPI::getRank() == 0)
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
