// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
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
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>

#include "IBFEVolumetricStabilization.h"

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Model data for the compressed block benchmark (Vadala-Roth et al. 2020, Sec. 4.1).
namespace ModelData
{
using namespace IBAMR::FEVolumetricStabilization;

static StabilizationFormulation stress_formulation = StabilizationFormulation::MODIFIED_INVARIANTS;
static double G_s = 80.194;
static double kappa_stab = 374.239;
static double eta_s = 4.0097;
static double kappa_d = 1.0e5;
static double traction_load = 200.0;
static double load_time = 1.0;

static constexpr double block_x_min = 10.0;
static constexpr double block_x_max = 30.0;
static constexpr double block_y_min = 15.0;
static constexpr double block_y_max = 25.0;
static constexpr double load_x_min = 15.0;
static constexpr double load_x_max = 25.0;

static double J_min_global = 1.0e8;
static double J_max_global = -1.0e8;

double
time_ramp(const double time)
{
    if (time < load_time)
    {
        const double tau = time / load_time;
        return -2.0 * tau * tau * tau + 3.0 * tau * tau;
    }
    return 1.0;
} // time_ramp

double
static_pressure_fcn(const double J)
{
    return computeLogarithmicStaticPressure(J, kappa_stab);
} // static_pressure_fcn

void
PK1_dev_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*x*/,
                        const libMesh::Point& /*X*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    const double J = FF.det();
    J_min_global = std::min(J_min_global, J);
    J_max_global = std::max(J_max_global, J);
    PP = computeNeoHookeanPK1(FF, G_s, stress_formulation, NDIM);
    return;
} // PK1_dev_stress_function

void
PK1_dil_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*x*/,
                        const libMesh::Point& /*X*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = computeLogarithmicPenaltyPK1(FF, kappa_stab, NDIM);
    return;
} // PK1_dil_stress_function

void
solid_body_force_function(VectorValue<double>& F,
                          const TensorValue<double>& /*FF*/,
                          const libMesh::Point& /*x*/,
                          const libMesh::Point& /*X*/,
                          Elem* const /*elem*/,
                          const vector<const vector<double>*>& var_data,
                          const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                          double /*time*/,
                          void* /*ctx*/)
{
    VectorValue<double> U;
    std::copy(var_data[0]->begin(), var_data[0]->end(), &U(0));
    F = -eta_s * U;
    return;
} // solid_body_force_function

libMesh::Point
compute_side_midpoint(const Elem* const elem, const unsigned short int side)
{
    libMesh::Point mid;
    unsigned int n_side_nodes = 0;
    for (unsigned int n = 0; n < elem->n_nodes(); ++n)
    {
        if (elem->is_node_on_side(n, side))
        {
            mid += elem->point(n);
            ++n_side_nodes;
        }
    }
    return mid / static_cast<double>(n_side_nodes);
} // compute_side_midpoint

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
                             const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                             double time,
                             void* /*ctx*/)
{
    F.zero();
    const double tol = 1.0e-8;
    const libMesh::Point side_mid = compute_side_midpoint(elem, side);
    if (std::abs(side_mid(1) - block_y_min) < tol)
    {
        // Zero vertical displacement along the bottom face.
        F(1) = kappa_d * (X(1) - x(1));
    }
    else if (std::abs(side_mid(1) - block_y_max) < tol)
    {
        // Zero horizontal displacement along top face and centered traction.
        F(0) = kappa_d * (X(0) - x(0));
        if (X(0) >= load_x_min - tol && X(0) <= load_x_max + tol)
        {
            F(1) = -time_ramp(time) * traction_load;
        }
    }
    return;
} // solid_surface_force_function

double
compute_block_area(MeshBase& mesh, const IBFEMethod* const ib_method_ops, EquationSystems* const equation_systems)
{
    System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
    NumericVector<double>* X_vec = X_system.solution.get();
    NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
    copy_and_synch(*X_vec, *X_ghost_vec);
    DofMap& X_dof_map = X_system.get_dof_map();
    vector<vector<unsigned int>> X_dof_indices(NDIM);
    std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<VectorValue<double>>>& dphi = fe->get_dphi();
    TensorValue<double> FF;
    boost::multi_array<double, 2> X_node;
    double area = 0.0;
    for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            jacobian(FF, qp, X_node, dphi);
            area += std::abs(FF.det()) * JxW[qp];
        }
    }
    return IBTK_MPI::sumReduction(area);
} // compute_block_area

double
compute_top_center_uy(const MeshBase& mesh,
                      const IBFEMethod* const ib_method_ops,
                      EquationSystems* const equation_systems,
                      const dof_id_type center_node_id)
{
    System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
    NumericVector<double>* X_vec = X_system.solution.get();
    NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
    copy_and_synch(*X_vec, *X_ghost_vec);
    const DofMap& X_dof_map = X_system.get_dof_map();
    vector<dof_id_type> X_idxs;
    double local_uy = -1.0e100;
    const Node& center_node = mesh.node_ref(center_node_id);
    if (center_node.processor_id() == static_cast<processor_id_type>(IBTK_MPI::getRank()))
    {
        X_dof_map.dof_indices(&center_node, X_idxs, 1);
        local_uy = (*X_ghost_vec)(X_idxs[0]) - center_node(1);
    }
    return IBTK_MPI::maxReduction(local_uy);
} // compute_top_center_uy
} // namespace ModelData
using namespace ModelData;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
#endif
        const string exodus_filename = app_initializer->getExodusIIFilename();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        ReplicatedMesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const string elem_type = input_db->getString("ELEM_TYPE");
        const int n_x = std::max(4, static_cast<int>(std::lround((block_x_max - block_x_min) / ds)));
        const int n_y = std::max(2, static_cast<int>(std::lround((block_y_max - block_y_min) / ds)));
        MeshTools::Generation::build_square(mesh,
                                            n_x,
                                            n_y,
                                            block_x_min,
                                            block_x_max,
                                            block_y_min,
                                            block_y_max,
                                            Utility::string_to_enum<ElemType>(elem_type));
        mesh.prepare_for_use();

        IBTK::StableCentroidPartitioner partitioner;
        partitioner.partition(mesh);

        dof_id_type top_center_node_id = 0;
        double min_dist = 1.0e100;
        const libMesh::Point target_center(20.0, block_y_max);
        for (auto nd_it = mesh.nodes_begin(); nd_it != mesh.nodes_end(); ++nd_it)
        {
            const Node* const node = *nd_it;
            const double dist = (*node - target_center).norm();
            if (dist < min_dist)
            {
                min_dist = dist;
                top_center_node_id = node->id();
            }
        }

        G_s = input_db->getDouble("G_S");
        const double nu_stab = input_db->getDouble("NU_STAB");
        kappa_stab = computeNumericalBulkModulus(G_s, nu_stab);
        eta_s = input_db->getDouble("ETA");
        kappa_d = input_db->getDouble("KAPPA_D");
        traction_load = input_db->getDouble("TRACTION");
        load_time = input_db->getDouble("LOAD_TIME");
        stress_formulation = stringToFormulation(input_db->getString("STRESS_FORMULATION"));
        const string volumetric_mode = input_db->getString("VOLUMETRIC_MODE");

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        IBFEMethod::PK1StressFcnData PK1_dev_data(PK1_dev_stress_function);
        PK1_dev_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        ib_method_ops->registerPK1StressFunction(PK1_dev_data);

        if (volumetric_mode == "PK1_PENALTY" && kappa_stab > 0.0)
        {
            IBFEMethod::PK1StressFcnData PK1_dil_data(PK1_dil_stress_function);
            PK1_dil_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
            ib_method_ops->registerPK1StressFunction(PK1_dil_data);
        }
        else if (volumetric_mode != "PK1_PENALTY" && volumetric_mode != "STATIC_PRESSURE" && volumetric_mode != "NONE")
        {
            TBOX_ERROR("Unsupported VOLUMETRIC_MODE: " << volumetric_mode << "\n");
        }

        const std::vector<int> vel_vars = { 0, 1 };
        const std::vector<SystemData> velocity_data = { SystemData(ib_method_ops->getVelocitySystemName(), vel_vars) };
        IBFEMethod::LagBodyForceFcnData body_force_data(solid_body_force_function, velocity_data);
        IBFEMethod::LagSurfaceForceFcnData surface_force_data(solid_surface_force_function);
        ib_method_ops->registerLagBodyForceFunction(body_force_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_force_data);

        ib_method_ops->initializeFEEquationSystems();
        if (volumetric_mode == "STATIC_PRESSURE" && kappa_stab > 0.0)
        {
            const PressureProjectionType proj_type = IBAMR::string_to_enum<PressureProjectionType>(
                input_db->getStringWithDefault("PRESSURE_PROJECTION_TYPE", "STABILIZED_PROJECTION"));
            ib_method_ops->registerStaticPressurePart(proj_type, static_pressure_fcn);
        }
        EquationSystems* const equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);
            const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);
            u_bc_coefs[d] = new muParserRobinBcCoefs(
                bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }
        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io = uses_exodus ? std::make_unique<ExodusII_IO>(mesh) : nullptr;
        if (uses_exodus)
        {
            exodus_io->append(RestartManager::getManager()->isFromRestart());
        }

        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        app_initializer.setNull();

        const double area_0 = compute_block_area(mesh, ib_method_ops, equation_systems);
        pout << "formulation=" << input_db->getString("STRESS_FORMULATION") << " nu_stab=" << nu_stab
             << " kappa_stab=" << kappa_stab << " volumetric_mode=" << volumetric_mode << " nodes=" << mesh.n_nodes()
             << " elems=" << mesh.n_elem() << " area_0=" << std::setprecision(10) << area_0 << "\n";

        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_io->write_timestep(exodus_filename, *equation_systems, 1, loop_time);
            }
        }

        const double loop_time_end = time_integrator->getEndTime();
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
            const double dt = time_integrator->getMaximumTimeStepSize();

            J_min_global = 1.0e8;
            J_max_global = -1.0e8;
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
            iteration_num += 1;

            J_min_global = IBTK_MPI::minReduction(J_min_global);
            J_max_global = IBTK_MPI::maxReduction(J_max_global);

            const bool last_step = !time_integrator->stepsRemaining() || IBTK::rel_equal_eps(loop_time, loop_time_end);
            if (iteration_num % 50 == 0 || last_step)
            {
                const double area = compute_block_area(mesh, ib_method_ops, equation_systems);
                const double rel_area_err = std::abs(area - area_0) / area_0;
                const double top_center_uy =
                    compute_top_center_uy(mesh, ib_method_ops, equation_systems, top_center_node_id);
                const double max_abs_J_minus_1 = std::max(std::abs(J_min_global - 1.0), std::abs(J_max_global - 1.0));
                pout << "step=" << iteration_num << " t=" << loop_time << " top_center_uy=" << top_center_uy
                     << " rel_area_err=" << rel_area_err << " J_min=" << J_min_global << " J_max=" << J_max_global
                     << " max_abs_J_minus_1=" << max_abs_J_minus_1 << "\n";
            }

            if (dump_viz_data && viz_dump_interval > 0 && (iteration_num % viz_dump_interval == 0 || last_step))
            {
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
            if (dump_restart_data && restart_dump_interval > 0 &&
                (iteration_num % restart_dump_interval == 0 || last_step))
            {
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && timer_dump_interval > 0 && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                TimerManager::getManager()->print(plog);
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
    return 0;
} // main
