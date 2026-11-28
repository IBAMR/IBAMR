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
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity and tethering model data for the 2D pressure-loaded elastic band.
namespace ModelData
{
static double G_s = 200.0;
static double kappa_s = 1.0e6;
static double J_band_min = 1.0e8;
static double J_band_max = -1.0e8;

void
tether_force_function(VectorValue<double>& F,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& X,
                      const libMesh::Point& s,
                      Elem* const elem,
                      const vector<const vector<double>*>& /*var_data*/,
                      const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                      double /*time*/,
                      void* /*ctx*/)
{
    if (elem->subdomain_id() == 0)
    {
        F = kappa_s * (s - X);
    }
    else
    {
        F.zero();
    }
    return;
} // tether_force_function

void
PK1_stress_function(TensorValue<double>& PP,
                    const TensorValue<double>& FF,
                    const libMesh::Point& /*X*/,
                    const libMesh::Point& /*s*/,
                    Elem* const elem,
                    const vector<const vector<double>*>& /*var_data*/,
                    const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                    double /*time*/,
                    void* /*ctx*/)
{
    if (elem->subdomain_id() == 1)
    {
        const double J = FF.det();
        J_band_min = std::min(J_band_min, J);
        J_band_max = std::max(J_band_max, J);
        // Unmodified neo-Hookean stress without volumetric stabilization.
        PP = G_s * FF;
    }
    else
    {
        PP.zero();
    }
    return;
} // PK1_stress_function

double
compute_band_area(MeshBase& mesh, const IBFEMethod* const ib_method_ops, EquationSystems* const equation_systems)
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
        if (elem->subdomain_id() != 1)
        {
            continue;
        }
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
} // compute_band_area

void
compute_displacements(const MeshBase& mesh,
                      const IBFEMethod* const ib_method_ops,
                      EquationSystems* const equation_systems,
                      const dof_id_type center_node_id,
                      double& center_x_disp,
                      double& max_block_disp)
{
    System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
    NumericVector<double>* X_vec = X_system.solution.get();
    NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
    copy_and_synch(*X_vec, *X_ghost_vec);
    const DofMap& X_dof_map = X_system.get_dof_map();
    vector<dof_id_type> X_idxs;

    double local_center_disp = -1.0e100;
    const Node& center_node = mesh.node_ref(center_node_id);
    if (center_node.processor_id() == static_cast<processor_id_type>(IBTK_MPI::getRank()))
    {
        X_dof_map.dof_indices(&center_node, X_idxs, 0);
        local_center_disp = (*X_ghost_vec)(X_idxs[0]) - center_node(0);
    }
    center_x_disp = IBTK_MPI::maxReduction(local_center_disp);

    double local_max_block_disp = 0.0;
    for (auto nd_it = mesh.local_nodes_begin(); nd_it != mesh.local_nodes_end(); ++nd_it)
    {
        const Node* const node = *nd_it;
        if ((*node)(1) < 0.1 - 1.0e-12 || (*node)(1) > 0.9 + 1.0e-12)
        {
            double disp_sq = 0.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(node, X_idxs, d);
                const double diff = (*X_ghost_vec)(X_idxs[0]) - (*node)(d);
                disp_sq += diff * diff;
            }
            local_max_block_disp = std::max(local_max_block_disp, std::sqrt(disp_sq));
        }
    }
    max_block_disp = IBTK_MPI::maxReduction(local_max_block_disp);
    return;
} // compute_displacements
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

        // Build a conformal mesh spanning [0.95, 1.05] x [0.0, 1.0] with aligned clamp interfaces at y = 0.1 and 0.9.
        ReplicatedMesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const string elem_type = input_db->getString("ELEM_TYPE");
        const int n_x = std::max(2, static_cast<int>(std::lround(0.1 / ds)));
        const int n_y_seg = std::max(2, static_cast<int>(std::lround(0.1 / ds)));
        const int n_y = 10 * n_y_seg;
        MeshTools::Generation::build_square(
            mesh, n_x, n_y, 0.95, 1.05, 0.0, 1.0, Utility::string_to_enum<ElemType>(elem_type));
        mesh.prepare_for_use();

        IBTK::StableCentroidPartitioner partitioner;
        partitioner.partition(mesh);

        for (auto el_it = mesh.active_elements_begin(); el_it != mesh.active_elements_end(); ++el_it)
        {
            Elem* const elem = *el_it;
            const double yc = elem->vertex_average()(1);
            if (yc < 0.1 || yc > 0.9)
            {
                elem->subdomain_id() = 0;
            }
            else
            {
                elem->subdomain_id() = 1;
            }
        }

        dof_id_type center_node_id = 0;
        double min_center_dist = 1.0e100;
        const libMesh::Point target_center(1.0, 0.5);
        for (auto nd_it = mesh.nodes_begin(); nd_it != mesh.nodes_end(); ++nd_it)
        {
            const Node* const node = *nd_it;
            const double dist = (*node - target_center).norm();
            if (dist < min_center_dist)
            {
                min_center_dist = dist;
                center_node_id = node->id();
            }
        }

        G_s = input_db->getDouble("G_S");
        kappa_s = input_db->getDouble("KAPPA_S");

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
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n");
        }

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

        IBFEMethod::PK1StressFcnData PK1_stress_data(PK1_stress_function);
        IBFEMethod::LagBodyForceFcnData tether_force_data(tether_force_function);
        ib_method_ops->registerPK1StressFunction(PK1_stress_data);
        ib_method_ops->registerLagBodyForceFunction(tether_force_data);
        ib_method_ops->initializeFEEquationSystems();
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

        const double band_area_0 = compute_band_area(mesh, ib_method_ops, equation_systems);
        pout << "kernel=" << input_db->getString("IB_DELTA_FUNCTION") << " elems=" << mesh.n_elem()
             << " band_area_0=" << std::setprecision(10) << band_area_0 << "\n";

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

            J_band_min = 1.0e8;
            J_band_max = -1.0e8;
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
            iteration_num += 1;

            J_band_min = IBTK_MPI::minReduction(J_band_min);
            J_band_max = IBTK_MPI::maxReduction(J_band_max);

            const bool last_step = !time_integrator->stepsRemaining() || IBTK::rel_equal_eps(loop_time, loop_time_end);
            if (iteration_num % 100 == 0 || last_step)
            {
                const double band_area = compute_band_area(mesh, ib_method_ops, equation_systems);
                const double rel_area_err = std::abs(band_area - band_area_0) / band_area_0;
                double center_x_disp = 0.0;
                double max_block_disp = 0.0;
                compute_displacements(
                    mesh, ib_method_ops, equation_systems, center_node_id, center_x_disp, max_block_disp);
                pout << "step=" << iteration_num << " t=" << loop_time << " center_x_disp=" << center_x_disp
                     << " rel_area_err=" << rel_area_err << " J_min=" << J_band_min << " J_max=" << J_band_max
                     << " max_block_disp=" << max_block_disp << "\n";
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
