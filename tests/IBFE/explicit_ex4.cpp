// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/xdr_io.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
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

/*
 * Test based on IBFE/explicit/ex4 that can verify whether or not a bunch of
 * things related to IBFEMethod work. This test program is run in various
 * different ways to check, e.g.,
 *
 * 1. Basic reproducibility of the IBFE stack (IBExplicitHierarchyIntegrator,
 *    INSStaggeredHierarchyIntegrator, and IBFEMethod)
 * 2. IBFE restart code
 * 3. IBFE part inactivation
 * 4. IBFECentroidPostProcessor
 * 5. IBTK::MergingLoadBalancer
 * 6. IBTK::MarkerPatchHierarchy
 */

// Elasticity model data.
namespace ModelData
{
// Coordinate mapping function.
void
coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    X(0) = s(0) + 0.6;
    X(1) = s(1) + 0.5;
#if (NDIM == 3)
    X(2) = s(2) + 0.5;
#endif
    return;
} // coordinate_mapping_function

// Stress tensor functions.
static double c1_s = 0.05;
static double p0_s = 0.0;
static double beta_s = 0.0;
void
PK1_dev_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = 2.0 * c1_s * FF;
    return;
} // PK1_dev_stress_function

void
PK1_dev_inactive_stress_function(TensorValue<double>& PP,
                                 const TensorValue<double>& FF,
                                 const libMesh::Point& /*X*/,
                                 const libMesh::Point& /*s*/,
                                 Elem* const /*elem*/,
                                 const vector<const vector<double>*>& /*var_data*/,
                                 const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                 double /*time*/,
                                 void* /*ctx*/)
{
    PP = 1e3 * c1_s * FF;
    return;
} // PK1_dev_inactive_stress_function

void
PK1_dil_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const vector<const vector<double>*>& /*var_data*/,
                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = 2.0 * (-p0_s + beta_s * log(FF.det())) * tensor_inverse_transpose(FF, NDIM);
    return;
} // PK1_dil_stress_function

static double uniform_body_source_strength = 0.0;
void
lag_body_source_function(double& Q,
                         const TensorValue<double>& /*FF*/,
                         const libMesh::Point& /*X*/,
                         const libMesh::Point& /*s*/,
                         Elem* /*elem*/,
                         const std::vector<const std::vector<double>*>& /*system_var_data*/,
                         const std::vector<const std::vector<VectorValue<double> >*>& /*system_grad_var_data*/,
                         double /*data_time*/,
                         void* /*ctx*/)
{
    Q = uniform_body_source_strength;
    return;
} // lag_body_source_function
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
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // suppress warnings caused by using a refinement ratio of 4 and not
    // setting up coarsening correctly
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-16");
    // use lower tolerances in 2D
    if (NDIM == 2)
    {
        PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-14");
        PetscOptionsSetValue(nullptr, "-stokes_ksp_rtol", "1e-14");
    }
    else
    {
        PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-12");
        PetscOptionsSetValue(nullptr, "-stokes_ksp_rtol", "1e-12");
    }

    { // cleanup dynamically allocated objects prior to shutdown
        // prevent a warning about timer initializations
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        // Create a simple FE mesh.
        ReplicatedMesh mesh(init.comm(), NDIM);
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 0.2;
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
            XdrIO io_in(mesh);
            io_in.read(std::string(SOURCE_DIR) + "/" + input_db->getString("mesh_file"));
            if (elem_type == "TRI6") mesh.all_second_order();
        }
        else
        {
            MeshTools::Generation::build_sphere(mesh, 0.01, 1, Utility::string_to_enum<ElemType>(elem_type));
        }
        mesh.prepare_for_use();

        // Possibly assign boundary elements to a finer patch level (depending
        // on what is in the input file)
        {
            libMesh::subdomain_id_type inner_id = 1;
            libMesh::subdomain_id_type outer_id = 2;
            if (input_db->getBoolWithDefault("use_huge_subdomains", false))
            {
                inner_id = std::numeric_limits<libMesh::subdomain_id_type>::max() - 3;
                outer_id = std::numeric_limits<libMesh::subdomain_id_type>::max() - 1;
            }
            const MeshBase::element_iterator el_end = mesh.active_elements_end();
            for (MeshBase::element_iterator el = mesh.active_elements_begin(); el != el_end; ++el)
            {
#if LIBMESH_VERSION_LESS_THAN(1, 7, 0)
                const libMesh::Point centroid = (*el)->centroid();
#else
                const libMesh::Point centroid = (*el)->vertex_average();
#endif
                if (centroid.norm() > 0.75 * R)
                {
                    (*el)->subdomain_id() = outer_id;
                }
                else
                {
                    (*el)->subdomain_id() = inner_id;
                }
            }
        }

        const bool use_amr = input_db->getBoolWithDefault("use_amr", false);
        if (use_amr)
        {
            const MeshBase::element_iterator el_end = mesh.elements_end();
            for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
            {
                Elem* elem = *el;
#if LIBMESH_VERSION_LESS_THAN(1, 7, 0)
                const libMesh::Point centroid = elem->centroid();
#else
                const libMesh::Point centroid = elem->vertex_average();
#endif

                if (centroid(1) > 0.0) elem->set_refinement_flag(Elem::REFINE);
            }
            MeshRefinement mesh_refinement(mesh);
            mesh_refinement.refine_and_coarsen_elements();
        }

        // Ensure nodes on the surface are on the analytic boundary.
#if NDIM == 2
        {
            const MeshBase::element_iterator el_end = mesh.elements_end();
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
        }
#else
        NULL_USE(R);
#endif

        ReplicatedMesh mesh_2(init.comm(), NDIM);
        // extra parameter controlling whether or not to add an inactive mesh
        // to test our ability to turn structures off and on
        const bool use_inactive_mesh = input_db->getBoolWithDefault("use_inactive_mesh", false);
        if (use_inactive_mesh)
        {
            // Add a substantial mesh near the top of the domain
#if NDIM == 2
            MeshTools::Generation::build_square(mesh_2, 10, 12, 0.0, 1.0, 0.15, 0.95, libMesh::QUAD9);
#else
            MeshTools::Generation::build_cube(mesh_2, 5, 6, 7, 0.0, 1.0, 0.0, 1.0, 0.15, 0.95, libMesh::HEX27);
#endif
        }

        // metis does a good job partitioning, but the partitioning relies on
        // random numbers: the seed changed in libMesh commit
        // 98cede90ca8837688ee13aac5e299a3765f083da (between 1.3.1 and
        // 1.4.0). Hence, to achieve consistent partitioning, use a simpler partitioning scheme:
        IBTK::StableCentroidPartitioner partitioner;
        partitioner.partition(mesh);
        if (use_inactive_mesh) partitioner.partition(mesh_2);

        c1_s = input_db->getDouble("C1_S");
        p0_s = input_db->getDouble("P0_S");
        beta_s = input_db->getDouble("BETA_S");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

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
        std::vector<libMesh::MeshBase*> meshes;
        meshes.push_back(&mesh);
        if (use_inactive_mesh) meshes.push_back(&mesh_2);
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           meshes,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
        Pointer<IBExplicitHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        time_integrator->registerLoadBalancer(load_balancer);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IBFE solver.
        ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
        PK1_dev_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        PK1_dil_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data);

        IBFEMethod::PK1StressFcnData PK1_dev_inactive_stress_data(PK1_dev_inactive_stress_function);
        PK1_dev_inactive_stress_data.quad_order = PK1_dev_stress_data.quad_order;
        if (use_inactive_mesh)
        {
            ib_method_ops->registerPK1StressFunction(PK1_dev_inactive_stress_data, 1);
        }
        if (input_db->getBoolWithDefault("ELIMINATE_PRESSURE_JUMPS", false))
        {
            ib_method_ops->registerStressNormalizationPart();
        }
        ib_method_ops->initializeFEEquationSystems();
        if (input_db->isDouble("BODY_SOURCE_STRENGTH"))
        {
            uniform_body_source_strength = input_db->getDouble("BODY_SOURCE_STRENGTH");
            ib_method_ops->registerLagBodySourceFunction(lag_body_source_function);
        }
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

        // Set up post processor to recover computed stresses.
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();

        FEDataManager* other_manager = nullptr;
        const bool use_separate_fe_data_manager = input_db->getBoolWithDefault("use_separate_fe_data_manager", false);
        if (use_separate_fe_data_manager)
        {
            FEDataManager::WorkloadSpec workload_spec;
            Pointer<Database> fe_data_manager_db(new InputDatabase("fe_data_manager_db"));
            other_manager = FEDataManager::getManager(
                fe_data_manager->getFEData(),
                "cloned_fe_data_manager",
                fe_data_manager_db,
                app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                fe_data_manager->getDefaultInterpSpec(),
                fe_data_manager->getDefaultSpreadSpec(),
                workload_spec);
            // Check that we have the same Lagrangian data.
            TBOX_ASSERT(fe_data_manager->getFEData() == other_manager->getFEData());
            TBOX_ASSERT(fe_data_manager->getEquationSystems() == other_manager->getEquationSystems());
        }
        else
        {
            other_manager = fe_data_manager;
        }

        const bool log_postprocessor = input_db->getBoolWithDefault("log_postprocessor", false);
        const int postprocessor_sampling_rate = input_db->getIntegerWithDefault("postprocessor_sampling_rate", 25);

        Pointer<IBFEPostProcessor> ib_post_processor =
            new IBFECentroidPostProcessor("IBFEPostProcessor", other_manager);
        ib_post_processor->registerTensorVariable("FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);

        std::pair<IBTK::TensorMeshFcnPtr, void*> PK1_dev_stress_fcn_data(PK1_dev_stress_function,
                                                                         static_cast<void*>(NULL));
        ib_post_processor->registerTensorVariable("sigma_dev",
                                                  MONOMIAL,
                                                  CONSTANT,
                                                  IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                                                  vector<SystemData>(),
                                                  &PK1_dev_stress_fcn_data);

        std::pair<IBTK::TensorMeshFcnPtr, void*> PK1_dil_stress_fcn_data(PK1_dil_stress_function,
                                                                         static_cast<void*>(NULL));
        ib_post_processor->registerTensorVariable("sigma_dil",
                                                  MONOMIAL,
                                                  CONSTANT,
                                                  IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                                                  vector<SystemData>(),
                                                  &PK1_dil_stress_fcn_data);

        Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
        HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill(
            /*data_idx*/ -1, "LINEAR_REFINE", /*use_cf_bdry_interpolation*/ false, "CONSERVATIVE_COARSEN", "LINEAR");
        FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR",
                                                QGAUSS,
                                                FIFTH,
                                                /*use_adaptive_quadrature*/ false,
                                                /*point_density*/ 2.0,
                                                /*use_consistent_mass_matrix*/ true,
                                                /*use_nodal_quadrature*/ false,
                                                /*allow_rules_with_negative_weights*/ false);
        ib_post_processor->registerInterpolatedScalarEulerianVariable(
            "p_f", LAGRANGE, FIRST, p_var, p_current_ctx, p_ghostfill, p_interp_spec);

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        if (ib_post_processor)
        {
            ib_post_processor->initializeFEData();
        }
        if (use_separate_fe_data_manager)
        {
            other_manager->setPatchHierarchy(patch_hierarchy);
        }

        auto add_markers = [&]()
        {
            System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
            NumericVector<double>& X_vec = *X_system.solution.get();
            std::vector<double> X_vec_global(X_vec.size());
            X_vec.localize(X_vec_global);

            EigenAlignedVector<IBTK::Point> positions;
            const auto n_nodes = mesh.parallel_n_nodes();
            std::vector<dof_id_type> X_idxs;
            for (dof_id_type i = 0; i < n_nodes; ++i)
            {
                const auto& node = mesh.node_ref(i);
                IBTK::Point X_node;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    IBTK::get_nodal_dof_indices(X_system.get_dof_map(), &node, d, X_idxs);
                    X_node[d] = X_vec_global[X_idxs[0]];
                }

                positions.push_back(X_node);
            }


            time_integrator->setMarkers(positions);
        };

        // First test for markers: add them at the start
        if (input_db->getBoolWithDefault("test_markers", false))
        {
            add_markers();
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        plog << std::setprecision(20);
        plog << "Total number of elems: " << mesh.n_elem() << std::endl;
        plog << "Total number of DoFs: " << equation_systems->n_dofs() << std::endl;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            if (use_inactive_mesh)
                if (input_db->getIntegerWithDefault("inactivate_timestep", 0) == iteration_num)
                    ib_method_ops->inactivateLagrangianStructure(1);

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

            iteration_num += 1;
            if (log_postprocessor && iteration_num % 10 == 0)
            {
                // The external FEDataManager is not reinitialized after
                // regrids, so recompute its Lagrangian-Eulerian data every
                // time we use it.
                if (use_separate_fe_data_manager) other_manager->reinitElementMappings();
                ib_post_processor->postProcessData(loop_time);
                // This hard-codes in an internal detail that is not presently documented
                auto& system = equation_systems->get_system("FF reconstruction system");
                plog << system.get_info() << std::endl;

                NumericVector<double>& solution = *system.solution;
                plog << "some values from \"FF reconstruction system\" + 1:" << std::endl;
                for (std::size_t i = solution.first_local_index(); i < solution.last_local_index(); ++i)
                {
                    if (i % (solution.local_size() / postprocessor_sampling_rate) == 0)
                    {
                        // improve the relative error (since we are taking
                        // derivatives it ends up being nearly 1e-5 for some
                        // values) somewhat by adding 1
                        plog << i << ", " << 1.0 + std::abs(solution(i)) << std::endl;
                    }
                }
            }

            if (dump_restart_data && (iteration_num % restart_dump_interval == 0))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }

            // Compute the volume of the structure.
            double J_integral = 0.0;
            System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
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
                const auto elem = *el_it;
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
            J_integral = IBTK_MPI::sumReduction(J_integral);
            if (IBTK_MPI::getRank() == 0)
            {
                plog << std::setprecision(12) << std::fixed << loop_time << " " << J_integral << std::endl;
            }

            // Second test for markers: add them after a restart (some users
            // want to add markers later in a simulation)
            if (iteration_num == 90 && input_db->getBoolWithDefault("test_markers_90", false))
            {
                add_markers();
            }
        }

        // Markers should still be in the same positions as nodes
        if (input_db->getBoolWithDefault("test_markers", false) || (iteration_num > 90 && input_db->getBoolWithDefault("test_markers_90", false)))
        {
            System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
            NumericVector<double>& X_vec = *X_system.solution.get();
            std::vector<double> X_vec_global(X_vec.size());
            X_vec.localize(X_vec_global);
            System& U_system = equation_systems->get_system<System>(ib_method_ops->getVelocitySystemName());
            NumericVector<double>& U_vec = *U_system.solution.get();
            std::vector<double> U_vec_global(U_vec.size());
            U_vec.localize(U_vec_global);

            const auto n_nodes = mesh.parallel_n_nodes();
            TBOX_ASSERT(X_vec_global.size() == NDIM * n_nodes);

            const auto pair = time_integrator->collectAllMarkers();
            TBOX_ASSERT(pair.first.size() == n_nodes);

            double max_distance = 0.0;
            double max_velocity_difference = 0.0;
            std::vector<dof_id_type> X_idxs;
            std::vector<dof_id_type> U_idxs;
            for (dof_id_type i = 0; i < n_nodes; ++i)
            {
                const auto& X_marker = pair.first[i];
                const auto& U_marker = pair.second[i];
                const auto& node = mesh.node_ref(i);
                IBTK::Point X_node;
                IBTK::Point U_node;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    IBTK::get_nodal_dof_indices(X_system.get_dof_map(), &node, d, X_idxs);
                    IBTK::get_nodal_dof_indices(U_system.get_dof_map(), &node, d, U_idxs);
                    X_node[d] = X_vec_global[X_idxs[0]];
                    U_node[d] = U_vec_global[U_idxs[0]];
                }

                // print the last node to make sure we aren't doing something
                // silly like 0 - 0.
                if (i == n_nodes - 1)
                {
                    plog << "Last node X = " << X_node[0] << ", " << X_node[1]
#if NDIM == 3
                         << ", " << X_node[2]
#endif
                         << '\n'
                         << "Last marker X = " << X_marker[0] << ", " << X_marker[1]
#if NDIM == 3
                         << ", " << X_marker[2]
#endif
                         << '\n'
                         << "Last marker U = " << U_marker[0] << ", " << U_marker[1]
#if NDIM == 3
                         << ", " << U_marker[2]
#endif
                         << '\n';
                }
                X_node -= X_marker;
                U_node -= U_marker;
                max_distance = std::max(max_distance, X_node.norm());
                max_velocity_difference = std::max(max_velocity_difference, U_node.norm());
            }

            plog << "Maximum distance between markers and nodes = " << max_distance << '\n';
            plog << "Maximum velocity difference between markers and nodes = " << max_velocity_difference << '\n';
        }

        if (input_db->getBoolWithDefault("log_scratch_partitioning", false))
        {
            Pointer<PatchHierarchy<NDIM> > scratch_hier = ib_method_ops->getScratchHierarchy();
            TBOX_ASSERT(scratch_hier);
            Pointer<PatchLevel<NDIM> > patch_level = scratch_hier->getPatchLevel(scratch_hier->getFinestLevelNumber());
            const BoxArray<NDIM> boxes = patch_level->getBoxes();
            plog << "Scratch hierarchy boxes:\n";
            for (int i = 0; i < boxes.size(); ++i) plog << boxes[i] << '\n';
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main
