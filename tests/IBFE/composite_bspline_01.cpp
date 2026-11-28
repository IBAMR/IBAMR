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
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/FEProjector.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

#include <fstream>
#include <iomanip>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Test driver for composite B-spline (CBS) kernels in IBFE with nodal quadrature.
namespace ModelData
{
static double G_s = 200.0;
static double kappa_s = 1.0e6;
static bool use_half_z = false;
static double J_min_global = 1.0e8;
static double J_max_global = -1.0e8;

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
    if (use_half_z)
    {
        if (s(2) > 0.25 - 1.0e-6)
        {
            F = kappa_s * (s - X);
        }
        else
        {
            F.zero();
        }
    }
    else if (elem->subdomain_id() == 0)
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
    if (use_half_z || elem->subdomain_id() == 1)
    {
        const double J = FF.det();
        J_min_global = std::min(J_min_global, J);
        J_max_global = std::max(J_max_global, J);
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
compute_elastic_volume(MeshBase& mesh, const IBFEMethod* const ib_method_ops, EquationSystems* const equation_systems)
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
    double vol = 0.0;
    for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
    {
        const Elem* const elem = *el_it;
        if (!use_half_z && elem->subdomain_id() != 1)
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
            vol += std::abs(FF.det()) * JxW[qp];
        }
    }
    return IBTK_MPI::sumReduction(vol);
} // compute_elastic_volume
} // namespace ModelData
using namespace ModelData;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-12");
    PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-14");

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        use_half_z = input_db->getBoolWithDefault("USE_HALF_Z", false);
        G_s = input_db->getDouble("G_S");
        kappa_s = input_db->getDouble("KAPPA_S");
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const string elem_type = input_db->getString("ELEM_TYPE");

        ReplicatedMesh mesh(init.comm(), NDIM);
#if (NDIM == 2)
        const int n_x = std::max(2, static_cast<int>(std::lround(0.1 / ds)));
        const int n_y_seg = std::max(2, static_cast<int>(std::lround(0.1 / ds)));
        MeshTools::Generation::build_square(
            mesh, n_x, 10 * n_y_seg, 0.95, 1.05, 0.0, 1.0, Utility::string_to_enum<ElemType>(elem_type));
#else
        const int n_c = std::max(2, static_cast<int>(std::lround(0.25 / ds)));
        MeshTools::Generation::build_cube(
            mesh, n_c, n_c, n_c, 0.35, 0.65, 0.35, 0.65, 0.0, 0.25, Utility::string_to_enum<ElemType>(elem_type));
#endif
        mesh.prepare_for_use();

        IBTK::StableCentroidPartitioner partitioner;
        partitioner.partition(mesh);

#if (NDIM == 2)
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
#else
        for (auto el_it = mesh.active_elements_begin(); el_it != mesh.active_elements_end(); ++el_it)
        {
            (*el_it)->subdomain_id() = 1;
        }
#endif

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"),
            false);
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           false);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator,
                                              false);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"), false);
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry, false);
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
                                        load_balancer,
                                        false);

        IBFEMethod::PK1StressFcnData PK1_stress_data(PK1_stress_function);
        IBFEMethod::LagBodyForceFcnData tether_force_data(tether_force_function);
        ib_method_ops->registerPK1StressFunction(PK1_stress_data);
        ib_method_ops->registerLagBodyForceFunction(tether_force_data);
        ib_method_ops->initializeFEEquationSystems();

        FEDataManager* const fe_data_manager = ib_method_ops->getFEDataManager();
        EquationSystems* const equation_systems = fe_data_manager->getEquationSystems();

        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);
            const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);
            u_bc_coefs[d] = new muParserRobinBcCoefs(
                bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }
        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        app_initializer.setNull();

        // Verify strictly positive diagonal L2 mass matrix entries under nodal quadrature.
        double min_diag_mass = 1.0e100;
        double max_diag_mass = -1.0e100;
        double diag_mass_sum = 0.0;
        {
            Pointer<Database> proj_db(new InputDatabase("proj_db"));
            proj_db->putBool("enable_logging", false);
            FEProjector fe_projector(equation_systems, proj_db);
            System& velocity_system = equation_systems->get_system(ib_method_ops->getVelocitySystemName());
            PetscVector<double>& diagonal_mass = *fe_projector.buildDiagonalL2MassMatrix(velocity_system.name());
            min_diag_mass = diagonal_mass.min();
            max_diag_mass = diagonal_mass.max();
            diag_mass_sum = diagonal_mass.sum();
        }
        TBOX_ASSERT(min_diag_mass > 0.0);

        const double vol_0 = compute_elastic_volume(mesh, ib_method_ops, equation_systems);
        double loop_time = time_integrator->getIntegratorTime();
        const double loop_time_end = time_integrator->getEndTime();
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            const double dt = time_integrator->getMaximumTimeStepSize();
            J_min_global = 1.0e8;
            J_max_global = -1.0e8;
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
            J_min_global = IBTK_MPI::minReduction(J_min_global);
            J_max_global = IBTK_MPI::maxReduction(J_max_global);
        }
        TBOX_ASSERT(J_min_global > 0.0);

        const double vol_final = compute_elastic_volume(mesh, ib_method_ops, equation_systems);
        const double rel_vol_err = std::abs(vol_final - vol_0) / vol_0;

        System& X_system = equation_systems->get_system<System>(ib_method_ops->getCurrentCoordinatesSystemName());
        NumericVector<double>* X_vec = X_system.solution.get();
        NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
        copy_and_synch(*X_vec, *X_ghost_vec);
        const DofMap& X_dof_map = X_system.get_dof_map();
        vector<dof_id_type> X_idxs;
        double local_max_disp = 0.0;
        double local_max_sym_z_disp = 0.0;
        for (auto nd_it = mesh.local_nodes_begin(); nd_it != mesh.local_nodes_end(); ++nd_it)
        {
            const Node* const node = *nd_it;
            double disp_sq = 0.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(node, X_idxs, d);
                const double diff = (*X_ghost_vec)(X_idxs[0]) - (*node)(d);
                disp_sq += diff * diff;
#if (NDIM == 3)
                if (use_half_z && d == 2 && std::abs((*node)(2)) < 1.0e-12)
                {
                    local_max_sym_z_disp = std::max(local_max_sym_z_disp, std::abs(diff));
                }
#endif
            }
            local_max_disp = std::max(local_max_disp, std::sqrt(disp_sq));
        }
        const double max_disp = IBTK_MPI::maxReduction(local_max_disp);
        const double max_sym_z_disp = IBTK_MPI::maxReduction(local_max_sym_z_disp);

        const int min_ghost_width = ib_method_ops->getMinimumGhostCellWidth().max();
        const int fe_ghost_width = fe_data_manager->getGhostCellWidth().max();

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "ndim = " << NDIM << "\n"
                << "kernel = " << input_db->getString("IB_DELTA_FUNCTION") << "\n"
                << "elem_type = " << elem_type << "\n"
                << "use_half_z = " << (use_half_z ? 1 : 0) << "\n"
                << "min_ghost_width = " << min_ghost_width << "\n"
                << "fe_ghost_width = " << fe_ghost_width << "\n"
                << "n_nodes = " << mesh.n_nodes() << "\n"
                << "n_elems = " << mesh.n_elem() << "\n"
                << std::scientific << std::setprecision(8) << "min_diag_mass = " << min_diag_mass << "\n"
                << "max_diag_mass = " << max_diag_mass << "\n"
                << "diag_mass_sum = " << diag_mass_sum << "\n"
                << "vol_0 = " << vol_0 << "\n"
                << "vol_final = " << vol_final << "\n"
                << "rel_vol_err = " << rel_vol_err << "\n"
                << "min_J = " << J_min_global << "\n"
                << "max_J = " << J_max_global << "\n"
                << "max_disp = " << max_disp << "\n";
            if (use_half_z)
            {
                out << "max_sym_z_disp = " << max_sym_z_disp << "\n";
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
    return 0;
} // main
