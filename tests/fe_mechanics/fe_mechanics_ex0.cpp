// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include <ibtk/AppInitializer.h>

#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include <SAMRAI_config.h>

#include <ibamr/app_namespaces.h>

using namespace IBAMR;
using namespace IBTK;
using namespace libMesh;
using namespace std;

// Elasticity model data.
namespace ModelData
{
static BoundaryInfo* lag_bdry_info;

// Tether (penalty) force function for the solid blocks.
static double kappa = 1.0e6;
static double damping = 0.0;
static double traction_force = 6.25e3;
static double load_time;

static bool use_static_pressure;
static bool use_dynamic_pressure;
static bool use_volumetric_term;
static std::string stress_funtion;

double time_ramp(double time);

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
                             double time,
                             void* /*ctx*/)
{
    const auto dim = elem->dim();
    if (lag_bdry_info->has_boundary_id(elem, side, dim == 2 ? 3 : 4))
    {
        F = kappa * (X - x);
    }
    else if (lag_bdry_info->has_boundary_id(elem, side, dim == 2 ? 1 : 2))
    {
        F = libMesh::Point(0.0, time_ramp(time) * traction_force, 0.0);
    }
    else
    {
        F.zero();
    }
}

static double shear_mod, bulk_mod;
void
PK1_dev_stress_function_mod(TensorValue<double>& PP,
                            const TensorValue<double>& FF,
                            const libMesh::Point& /*x*/,
                            const libMesh::Point& /*X*/,
                            Elem* const /*elem*/,
                            const vector<const vector<double>*>& /*var_data*/,
                            const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double /*time*/,
                            void* /*ctx*/)
{
    const RealTensor CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    const RealTensor dI1_bar_dFF = pow(FF.det(), -2.0 / 3.0) * (FF - I1 / 3.0 * tensor_inverse_transpose(FF, NDIM));
    PP = shear_mod * dI1_bar_dFF;
}

void
PK1_dev_stress_function_unmod(TensorValue<double>& PP,
                              const TensorValue<double>& FF,
                              const libMesh::Point& /*x*/,
                              const libMesh::Point& /*X*/,
                              Elem* const /*elem*/,
                              const vector<const vector<double>*>& /*var_data*/,
                              const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                              double /*time*/,
                              void* /*ctx*/)
{
    PP = shear_mod * FF;
}

void
PK1_dev_stress_function_dev(TensorValue<double>& PP,
                            const TensorValue<double>& FF,
                            const libMesh::Point& /*x*/,
                            const libMesh::Point& /*X*/,
                            Elem* const /*elem*/,
                            const vector<const vector<double>*>& /*var_data*/,
                            const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double /*time*/,
                            void* /*ctx*/)
{
    const RealTensor CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    PP = shear_mod * (FF - I1 / 3.0 * tensor_inverse_transpose(FF, NDIM));
}

void
PK1_dil_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*x*/,
                        const libMesh::Point& /*X*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = bulk_mod * log(FF.det()) * tensor_inverse_transpose(FF, NDIM);
}

void
solid_body_force_function(VectorValue<double>& F,
                          const TensorValue<double>& /*FF*/,
                          const libMesh::Point& /*x*/,
                          const libMesh::Point& /*X*/,
                          Elem* const /*elem*/,
                          const vector<const vector<double>*>& var_data,
                          const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                          double /*time*/,
                          void* /*ctx*/)
{
    VectorValue<double> U;
    std::copy(var_data[0]->begin(), var_data[0]->end(), &U(0));
    F = -damping * U;
}

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

double
time_ramp(double time)
{
    double displacement_factor;
    if (time < load_time)
        displacement_factor = -2.0 * (time / load_time) * (time / load_time) * (time / load_time) +
                              3.0 * (time / load_time) * (time / load_time);
    else
        displacement_factor = 1.0;

    return displacement_factor;
}

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
    // Initialize libMesh and SAMRAI.
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

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        Mesh mesh(init.comm(), NDIM);
        if (input_db->keyExists("MESH_NAME"))
        {
            mesh.read(input_db->getString("MESH_NAME"));
        }
        else
        {
            const string elem_type = input_db->getStringWithDefault("ELEM_TYPE", "TRI3");
            int ps_x = input_db->getDouble("N_X");
            int ps_y = (NDIM > 1 ? input_db->getDouble("N_Y") : 0);
            int ps_z = (NDIM > 2 ? input_db->getDouble("N_Z") : 0);

            // Build grid mesh
            MeshTools::Generation::build_cube(mesh,
                                              ps_x,
                                              (NDIM > 1) ? ps_y : 0,
                                              (NDIM > 2) ? ps_z : 0,
                                              0.0,
                                              1.0,
                                              0.0,
                                              1.0,
                                              0.0,
                                              1.0,
                                              Utility::string_to_enum<ElemType>(elem_type));

            // Map unit square onto cook's membrane
            MeshBase::const_node_iterator nd = mesh.nodes_begin();
            const MeshBase::const_node_iterator end_nd = mesh.nodes_end();
            for (; nd != end_nd; ++nd)
            {
                libMesh::Point s = **nd;
                (**nd)(0) = 4.8 * s(0) + 2.6;
                (**nd)(1) = -2.8 * s(0) * s(1) + 4.4 * s(0) + 4.4 * s(1) + 2.0;
            }
        }

        mesh.prepare_for_use();
        lag_bdry_info = &mesh.get_boundary_info();

        // Setup the model parameters.
        shear_mod = input_db->getDouble("SHEAR_MOD");
        bulk_mod = input_db->getDouble("BULK_MOD");
        damping = input_db->getDouble("ETA");
        kappa = input_db->getDouble("KAPPA");
        traction_force = input_db->getDouble("TRACTION_FORCE");
        load_time = input_db->getDouble("LOAD_TIME");

        use_static_pressure = input_db->getBool("USE_STATIC_PRESSURE");
        use_dynamic_pressure = input_db->getBool("USE_DYNAMIC_PRESSURE");
        PressureProjectionType pressure_proj_type = UNKNOWN_PRESSURE_TYPE;
        if (use_static_pressure || use_dynamic_pressure)
        {
            pressure_proj_type =
                IBAMR::string_to_enum<PressureProjectionType>(input_db->getString("PRESSURE_PROJECTION_TYPE"));
        }
        stress_funtion = input_db->getString("STRESS_FUNCTION");
        use_volumetric_term = input_db->getBool("USE_VOLUMETRIC_TERM") && !use_static_pressure && !use_dynamic_pressure;

        // Setup the time stepping parameters.
        const double loop_time_end = input_db->getDouble("END_TIME");
        const double dt = input_db->getDouble("DT");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<FEMechanicsExplicitIntegrator> fem_solver =
            new FEMechanicsExplicitIntegrator("FEMechanicsExplicitIntegrator",
                                              app_initializer->getComponentDatabase("FEMechanicsExplicitIntegrator"),
                                              &mesh);

        // attach velocity
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> velocity_data(1);
        velocity_data[0] = SystemData(FEMechanicsBase::VELOCITY_SYSTEM_NAME, vars);

        // Configure the FE solver.
        FEMechanicsBase::LagSurfaceForceFcnData solid_surface_force_data(solid_surface_force_function);
        FEMechanicsBase::LagBodyForceFcnData solid_body_force_data(solid_body_force_function, velocity_data);
        fem_solver->registerLagSurfaceForceFunction(solid_surface_force_data);
        fem_solver->registerLagBodyForceFunction(solid_body_force_data);

        if (stress_funtion == "MODIFIED")
        {
            FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_mod);
            PK1_dev_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
            fem_solver->registerPK1StressFunction(PK1_dev_stress_data);
        }
        else if (stress_funtion == "UNMODIFIED")
        {
            FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_unmod);
            PK1_dev_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
            fem_solver->registerPK1StressFunction(PK1_dev_stress_data);
        }
        else if (stress_funtion == "DEVIATORIC")
        {
            FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_dev);
            PK1_dev_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
            fem_solver->registerPK1StressFunction(PK1_dev_stress_data);
        }
        else
        {
            FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_mod);
            PK1_dev_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
            fem_solver->registerPK1StressFunction(PK1_dev_stress_data);
        }

        if (use_volumetric_term == true)
        {
            FEMechanicsBase::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
            PK1_dil_stress_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
            fem_solver->registerPK1StressFunction(PK1_dil_stress_data);
        }

        fem_solver->initializeFEEquationSystems();
        if (use_static_pressure)
        {
            fem_solver->registerStaticPressurePart(pressure_proj_type);
        }
        if (use_dynamic_pressure)
        {
            fem_solver->registerDynamicPressurePart(pressure_proj_type);
        }
        EquationSystems* equation_systems = fem_solver->getEquationSystems();
        ExplicitSystem& jac_system = equation_systems->add_system<ExplicitSystem>("JacobianDeterminant");
        jac_system.add_variable("Avg J", CONSTANT, MONOMIAL);
        jac_system.attach_init_function(apply_initial_jacobian);

        const std::string time_stepping_scheme =
            input_db->getStringWithDefault("TIME_STEPPING_SCHEME", "MODIFIED_TRAPEZOIDAL_RULE");

        // Set up visualization plot file writers.
        unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

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
            if (uses_exodus)
            {
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save volume of structure.
        ofstream volume_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            volume_stream.open("volume.curve", ios_base::out | ios_base::trunc);
        }

        // Main time step loop.
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end))
        {
            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            // to compute average J for each element
            double J_integral = 0.0;
            double L_1_error_J = 0.0;
            double L_2_error_J = 0.0;
            double L_oo_temp = 0.0;
            double L_oo_error_J = 0.0;
            double elem_area = 0.0;
            System& X_system = equation_systems->get_system<System>(FEMechanicsBase::COORDS_SYSTEM_NAME);
            NumericVector<double>* X_vec = X_system.solution.get();
            NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
            X_vec->localize(*X_ghost_vec);
            DofMap& X_dof_map = X_system.get_dof_map();
            vector<vector<unsigned int> > X_dof_indices(NDIM);
            unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
            unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
            fe->attach_quadrature_rule(qrule.get());
            const vector<double>& JxW = fe->get_JxW();
            const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();
            TensorValue<double> FF;
            boost::multi_array<double, 2> X_node;
            const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
            const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
            for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
            {
                elem_area = 0.0;
                L_oo_temp = 0.0;
                J_integral = 0.0;

                Elem* const elem = *el_it;
                fe->reinit(elem);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_dof_map.dof_indices(elem, X_dof_indices[d], d);
                }
                // std::cout << "\n Current element is " << elem;
                const int n_qp = qrule->n_points();
                get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    jacobian(FF, qp, X_node, dphi);
                    J_integral += abs(FF.det()) * JxW[qp];
                    L_1_error_J += abs(1.0 - FF.det()) * JxW[qp];
                    L_2_error_J += (1.0 - FF.det()) * (1.0 - FF.det()) * JxW[qp];
                    L_oo_temp += abs(1.0 - FF.det()) * JxW[qp];
                    // if (abs(1.0 - FF.det()) > L_oo_error_J ) L_oo_error_J = abs(1.0 - FF.det());
                    elem_area += JxW[qp];

                    // std::cout << "\n At qp, abs error in J is " << abs(1.0 - FF.det());
                }
                // std::cout << "\n";
                L_oo_temp = L_oo_temp / elem_area;
                if (L_oo_temp > L_oo_error_J) L_oo_error_J = L_oo_temp;
                J_integral = J_integral / elem_area;
                ExplicitSystem& jac_system = equation_systems->get_system<ExplicitSystem>("JacobianDeterminant");
                const unsigned int jac_system_num = jac_system.number(); // identifier number for the system
                const dof_id_type dof_id_J = elem->dof_number(jac_system_num, 0, /*comp*/ 0);
                jac_system.solution->set(dof_id_J, J_integral);
                // jac_system.solution->localize(*jac_system.current_local_solution);
            }
            jac_system.solution->close();
            L_1_error_J = SAMRAI_MPI::sumReduction(L_1_error_J);
            L_2_error_J = SAMRAI_MPI::sumReduction(L_2_error_J);
            L_oo_error_J = SAMRAI_MPI::maxReduction(L_oo_error_J);
            if (SAMRAI_MPI::getRank() == 0) L_2_error_J = sqrt(L_2_error_J);

            pout << "\nThe L_1 error of J = " << L_1_error_J << "\n";
            pout << "\nThe L_2 error of J = " << L_2_error_J << "\n";
            pout << "\nThe L_oo error of J = " << L_oo_error_J << "\n";
            pout << "\n";

            fem_solver->preprocessIntegrateData(loop_time, loop_time + dt, /*num_cycles*/ 1);
            if (time_stepping_scheme == "FORWARD_EULER")
            {
                fem_solver->forwardEulerStep(loop_time, loop_time + dt);
            }
            else if (time_stepping_scheme == "MODIFIED_EULER")
            {
                fem_solver->modifiedEulerStep(loop_time, loop_time + dt);
            }
            else if (time_stepping_scheme == "MIDPOINT")
            {
                fem_solver->midpointStep(loop_time, loop_time + dt);
            }
            else if (time_stepping_scheme == "TRAPEZOIDAL_RULE" || time_stepping_scheme == "SSPRK2")
            {
                fem_solver->trapezoidalStep(loop_time, loop_time + dt);
            }
            else if (time_stepping_scheme == "MODIFIED_TRAPEZOIDAL_RULE")
            {
                fem_solver->modifiedTrapezoidalStep(loop_time, loop_time + dt);
            }
            else
            {
                TBOX_ERROR("unknown time stepping scheme: " << time_stepping_scheme << "\n");
            }
            fem_solver->postprocessIntegrateData(loop_time, loop_time + dt, /*num_cycles*/ 1);

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
                if (uses_exodus)
                {
                    exodus_io->write_timestep(
                        exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
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

        } // end of time loop

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main
