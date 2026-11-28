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
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>

#include "IBFEVolumetricStabilization.h"

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Test driver for Vadala-Roth (2020) volumetric stabilization in IBFE.
namespace ModelData
{
using namespace IBAMR::FEVolumetricStabilization;

static string material_model = "NEO_HOOKEAN";
static StabilizationFormulation stress_formulation = StabilizationFormulation::MODIFIED_INVARIANTS;
static double G_s = 80.194;
static double c1_mr = 30.0;
static double c2_mr = 10.097;
static double G_T = 80.194;
static double G_L = 120.0;
static double E_L = 400.0;
static libMesh::VectorValue<double> fiber_dir(1.0, 0.5, 0.25);
static double kappa_stab = 374.239;
static double eta_s = 4.0097;
static double kappa_d = 1.0e5;
static double traction_load = 200.0;
static double load_time = 0.05;

static constexpr double block_x_min = 10.0;
static constexpr double block_x_max = 30.0;
static constexpr double block_y_min = 15.0;
static constexpr double block_y_max = 25.0;
static constexpr double block_z_min = 15.0;
static constexpr double block_z_max = 25.0;
static constexpr double load_x_min = 15.0;
static constexpr double load_x_max = 25.0;

static double J_min_global = 1.0e8;
static double J_max_global = -1.0e8;

double
tensor_inf_norm(const TensorValue<double>& A, const unsigned int dim = 3)
{
    double max_abs = 0.0;
    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = 0; j < dim; ++j)
        {
            max_abs = std::max(max_abs, std::abs(A(i, j)));
        }
    }
    return max_abs;
} // tensor_inf_norm

TensorValue<double>
finite_difference_pk1(const TensorValue<double>& FF,
                      const std::function<double(const TensorValue<double>&)>& energy_fcn)
{
    const double eps = 1.0e-6;
    TensorValue<double> PP_fd;
    PP_fd.zero();
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int j = 0; j < NDIM; ++j)
        {
            TensorValue<double> FF_plus = FF;
            TensorValue<double> FF_minus = FF;
            FF_plus(i, j) += eps;
            FF_minus(i, j) -= eps;
            PP_fd(i, j) = (energy_fcn(FF_plus) - energy_fcn(FF_minus)) / (2.0 * eps);
        }
    }
    return PP_fd;
} // finite_difference_pk1

bool
run_unit_checks()
{
    const double bulk_unmod = computeNumericalBulkModulus(80.194, -1.0);
    TBOX_ASSERT(std::abs(bulk_unmod) < 1.0e-14);
    const double bulk_stab = computeNumericalBulkModulus(80.194, 0.4);
    TBOX_ASSERT(std::abs(bulk_stab - 374.2386666666667) < 1.0e-11);

    TensorValue<double> II;
    II.zero();
    II(0, 0) = 1.0;
    II(1, 1) = 1.0;
    II(2, 2) = 1.0;

    TensorValue<double> FF;
    FF.zero();
    FF(0, 0) = 1.15;
    FF(0, 1) = 0.12;
    FF(1, 0) = -0.08;
    FF(1, 1) = 0.92;
#if (NDIM == 3)
    FF(0, 2) = -0.05;
    FF(1, 2) = 0.09;
    FF(2, 0) = 0.04;
    FF(2, 1) = -0.06;
    FF(2, 2) = 1.08;
#else
    FF(2, 2) = 1.0;
#endif

    const libMesh::VectorValue<double> A_test(1.0, 2.0, NDIM == 3 ? -1.0 : 0.0);
    const double mr_G = computeMooneyRivlinShearModulus(c1_mr, c2_mr);
    TBOX_ASSERT(std::abs(mr_G - 80.194) < 1.0e-13);

    // Check zero stress in reference configuration I.
    for (const auto form :
         { StabilizationFormulation::MODIFIED_INVARIANTS, StabilizationFormulation::DEVIATORIC_PROJECTION })
    {
        TBOX_ASSERT(tensor_inf_norm(computeNeoHookeanPK1(II, G_s, form, NDIM)) < 1.0e-14);
        TBOX_ASSERT(tensor_inf_norm(computeMooneyRivlinPK1(II, c1_mr, c2_mr, form, NDIM)) < 1.0e-14);
        TBOX_ASSERT(tensor_inf_norm(computeModifiedReinforcingPK1(II, G_T, G_L, E_L, A_test, form, NDIM)) < 1.0e-13);
    }
    TBOX_ASSERT(tensor_inf_norm(computeLogarithmicPenaltyPK1(II, bulk_stab, NDIM)) < 1.0e-14);
    TBOX_ASSERT(tensor_inf_norm(computeLiuHofstetterMangPK1(II, bulk_stab, NDIM)) < 1.0e-14);

    // Check traceless Cauchy stress for modified invariants and deviatoric projection.
    const double J = FF.det();
    for (const auto form :
         { StabilizationFormulation::MODIFIED_INVARIANTS, StabilizationFormulation::DEVIATORIC_PROJECTION })
    {
        const TensorValue<double> PP_nh = computeNeoHookeanPK1(FF, G_s, form, NDIM);
        const TensorValue<double> PP_mr = computeMooneyRivlinPK1(FF, c1_mr, c2_mr, form, NDIM);
        TBOX_ASSERT(std::abs((PP_nh * FF.transpose()).tr() / J) < 1.0e-12);
        TBOX_ASSERT(std::abs((PP_mr * FF.transpose()).tr() / J) < 1.0e-12);
    }

    // Check analytic PK1 derivatives against central finite differences.
    for (const auto form : { StabilizationFormulation::UNMODIFIED, StabilizationFormulation::MODIFIED_INVARIANTS })
    {
        const TensorValue<double> P_nh = computeNeoHookeanPK1(FF, G_s, form, NDIM);
        const TensorValue<double> P_nh_fd = finite_difference_pk1(
            FF, [&](const TensorValue<double>& F_in) { return computeNeoHookeanEnergy(F_in, G_s, form); });
        TBOX_ASSERT(tensor_inf_norm(P_nh - P_nh_fd, NDIM) < 1.0e-8 * std::max(1.0, tensor_inf_norm(P_nh, NDIM)));

        const TensorValue<double> P_mr = computeMooneyRivlinPK1(FF, c1_mr, c2_mr, form, NDIM);
        const TensorValue<double> P_mr_fd = finite_difference_pk1(
            FF, [&](const TensorValue<double>& F_in) { return computeMooneyRivlinEnergy(F_in, c1_mr, c2_mr, form); });
        TBOX_ASSERT(tensor_inf_norm(P_mr - P_mr_fd, NDIM) < 1.0e-8 * std::max(1.0, tensor_inf_norm(P_mr, NDIM)));

        const TensorValue<double> P_rf = computeModifiedReinforcingPK1(FF, G_T, G_L, E_L, A_test, form, NDIM);
        const TensorValue<double> P_rf_fd =
            finite_difference_pk1(FF,
                                  [&](const TensorValue<double>& F_in)
                                  { return computeModifiedReinforcingEnergy(F_in, G_T, G_L, E_L, A_test, form); });
        TBOX_ASSERT(tensor_inf_norm(P_rf - P_rf_fd, NDIM) < 1.0e-8 * std::max(1.0, tensor_inf_norm(P_rf, NDIM)));
    }

    const TensorValue<double> P_log = computeLogarithmicPenaltyPK1(FF, bulk_stab, NDIM);
    const TensorValue<double> P_log_fd = finite_difference_pk1(
        FF, [&](const TensorValue<double>& F_in) { return computeLogarithmicPenaltyEnergy(F_in.det(), bulk_stab); });
    TBOX_ASSERT(tensor_inf_norm(P_log - P_log_fd, NDIM) < 1.0e-8 * std::max(1.0, tensor_inf_norm(P_log, NDIM)));

    const TensorValue<double> P_lhm = computeLiuHofstetterMangPK1(FF, bulk_stab, NDIM);
    const TensorValue<double> P_lhm_fd = finite_difference_pk1(
        FF, [&](const TensorValue<double>& F_in) { return computeLiuHofstetterMangEnergy(F_in.det(), bulk_stab); });
    TBOX_ASSERT(tensor_inf_norm(P_lhm - P_lhm_fd, NDIM) < 1.0e-8 * std::max(1.0, tensor_inf_norm(P_lhm, NDIM)));

    // Check FEMechanicsBase static pressure sign convention PP = -J * P * F^{-T}.
    const TensorValue<double> FF_inv_trans = IBTK::tensor_inverse_transpose(FF, NDIM);
    const TensorValue<double> PP_from_log_P = -J * computeLogarithmicStaticPressure(J, bulk_stab) * FF_inv_trans;
    const TensorValue<double> PP_from_lhm_P = -J * computeLiuHofstetterMangStaticPressure(J, bulk_stab) * FF_inv_trans;
    TBOX_ASSERT(tensor_inf_norm(PP_from_log_P - P_log) < 1.0e-12);
    TBOX_ASSERT(tensor_inf_norm(PP_from_lhm_P - P_lhm) < 1.0e-12);
    return true;
} // run_unit_checks

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
    if (material_model == "NEO_HOOKEAN")
    {
        PP = computeNeoHookeanPK1(FF, G_s, stress_formulation, NDIM);
    }
    else if (material_model == "MOONEY_RIVLIN")
    {
        PP = computeMooneyRivlinPK1(FF, c1_mr, c2_mr, stress_formulation, NDIM);
    }
    else if (material_model == "REINFORCING")
    {
        PP = computeModifiedReinforcingPK1(FF, G_T, G_L, E_L, fiber_dir, stress_formulation, NDIM);
    }
    else
    {
        TBOX_ERROR("Unsupported MATERIAL_MODEL: " << material_model << "\n");
    }
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
        F(1) = kappa_d * (X(1) - x(1));
    }
    else if (std::abs(side_mid(1) - block_y_max) < tol)
    {
        F(0) = kappa_d * (X(0) - x(0));
        if (X(0) >= load_x_min - tol && X(0) <= load_x_max + tol)
        {
            F(1) = -time_ramp(time) * traction_load;
        }
    }
    return;
} // solid_surface_force_function

double
compute_structure_volume(MeshBase& mesh, const IBFEMethod* const ib_method_ops, EquationSystems* const equation_systems)
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
    double volume = 0.0;
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
            volume += std::abs(FF.det()) * JxW[qp];
        }
    }
    return IBTK_MPI::sumReduction(volume);
} // compute_structure_volume

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

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-12");
    PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-14");

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const bool run_unit = input_db->getBoolWithDefault("RUN_UNIT_CHECKS", false);
        const bool unit_checks_passed = run_unit ? run_unit_checks() : true;

        material_model = input_db->getString("MATERIAL_MODEL");
        stress_formulation = stringToFormulation(input_db->getString("STRESS_FORMULATION"));
        const string volumetric_mode = input_db->getString("VOLUMETRIC_MODE");
        const string proj_type_str = input_db->getStringWithDefault("PRESSURE_PROJECTION_TYPE", "NONE");
        const bool use_default_static_pressure = input_db->getBoolWithDefault("USE_DEFAULT_STATIC_PRESSURE", false);

        G_s = input_db->getDoubleWithDefault("G_S", 80.194);
        c1_mr = input_db->getDoubleWithDefault("C1_MR", 30.0);
        c2_mr = input_db->getDoubleWithDefault("C2_MR", 10.097);
        G_T = input_db->getDoubleWithDefault("G_T", 80.194);
        G_L = input_db->getDoubleWithDefault("G_L", 120.0);
        E_L = input_db->getDoubleWithDefault("E_L", 400.0);
        eta_s = input_db->getDouble("ETA");
        kappa_d = input_db->getDouble("KAPPA_D");
        traction_load = input_db->getDouble("TRACTION");
        load_time = input_db->getDouble("LOAD_TIME");

        const double nu_stab = input_db->getDouble("NU_STAB");
        double base_shear = G_s;
        if (material_model == "MOONEY_RIVLIN")
        {
            base_shear = computeMooneyRivlinShearModulus(c1_mr, c2_mr);
        }
        else if (material_model == "REINFORCING")
        {
            base_shear = G_T;
        }
        kappa_stab = computeNumericalBulkModulus(base_shear, nu_stab);

        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const string elem_type = input_db->getString("ELEM_TYPE");
        ReplicatedMesh mesh(init.comm(), NDIM);
        const int n_x = std::max(2, static_cast<int>(std::lround((block_x_max - block_x_min) / ds)));
        const int n_y = std::max(2, static_cast<int>(std::lround((block_y_max - block_y_min) / ds)));
#if (NDIM == 2)
        MeshTools::Generation::build_square(mesh,
                                            n_x,
                                            n_y,
                                            block_x_min,
                                            block_x_max,
                                            block_y_min,
                                            block_y_max,
                                            Utility::string_to_enum<ElemType>(elem_type));
#else
        const int n_z = std::max(2, static_cast<int>(std::lround((block_z_max - block_z_min) / ds)));
        MeshTools::Generation::build_cube(mesh,
                                          n_x,
                                          n_y,
                                          n_z,
                                          block_x_min,
                                          block_x_max,
                                          block_y_min,
                                          block_y_max,
                                          block_z_min,
                                          block_z_max,
                                          Utility::string_to_enum<ElemType>(elem_type));
#endif
        mesh.prepare_for_use();

        IBTK::StableCentroidPartitioner partitioner;
        partitioner.partition(mesh);

        dof_id_type top_center_node_id = 0;
        double min_dist = 1.0e100;
#if (NDIM == 2)
        const libMesh::Point target_center(20.0, block_y_max);
#else
        const libMesh::Point target_center(20.0, block_y_max, 20.0);
#endif
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

        IBFEMethod::PK1StressFcnData PK1_dev_data(PK1_dev_stress_function);
        PK1_dev_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        ib_method_ops->registerPK1StressFunction(PK1_dev_data);

        if (volumetric_mode == "PK1_PENALTY")
        {
            IBFEMethod::PK1StressFcnData PK1_dil_data(PK1_dil_stress_function);
            PK1_dil_data.quad_order =
                Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
            ib_method_ops->registerPK1StressFunction(PK1_dil_data);
        }
        else if (volumetric_mode != "STATIC_PRESSURE" && volumetric_mode != "NONE")
        {
            TBOX_ERROR("Unsupported VOLUMETRIC_MODE: " << volumetric_mode << "\n");
        }

        vector<int> vel_vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            vel_vars[d] = static_cast<int>(d);
        }
        const vector<SystemData> velocity_data = { SystemData(ib_method_ops->getVelocitySystemName(), vel_vars) };
        IBFEMethod::LagBodyForceFcnData body_force_data(solid_body_force_function, velocity_data);
        IBFEMethod::LagSurfaceForceFcnData surface_force_data(solid_surface_force_function);
        ib_method_ops->registerLagBodyForceFunction(body_force_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_force_data);

        ib_method_ops->initializeFEEquationSystems();
        if (volumetric_mode == "STATIC_PRESSURE")
        {
            const PressureProjectionType proj_type = IBAMR::string_to_enum<PressureProjectionType>(proj_type_str);
            ib_method_ops->registerStaticPressurePart(proj_type,
                                                      use_default_static_pressure ? nullptr : static_pressure_fcn);
        }

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

        const double vol_0 = compute_structure_volume(mesh, ib_method_ops, equation_systems);
        double loop_time = time_integrator->getIntegratorTime();
        const double loop_time_end = time_integrator->getEndTime();
        int n_steps = 0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            const double dt = time_integrator->getMaximumTimeStepSize();
            J_min_global = 1.0e8;
            J_max_global = -1.0e8;
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
            ++n_steps;
            J_min_global = IBTK_MPI::minReduction(J_min_global);
            J_max_global = IBTK_MPI::maxReduction(J_max_global);
        }
        TBOX_ASSERT(J_min_global > 0.0);

        const double vol_final = compute_structure_volume(mesh, ib_method_ops, equation_systems);
        const double rel_vol_err = std::abs(vol_final - vol_0) / vol_0;
        const double top_center_uy = compute_top_center_uy(mesh, ib_method_ops, equation_systems, top_center_node_id);

        double max_abs_P = 0.0;
        if (volumetric_mode == "STATIC_PRESSURE")
        {
            ExplicitSystem& P_system =
                equation_systems->get_system<ExplicitSystem>(ib_method_ops->getPressureSystemName());
            max_abs_P = P_system.solution->linfty_norm();
            TBOX_ASSERT(max_abs_P > 0.0);
        }

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "ndim = " << NDIM << "\n"
                << "elem_type = " << elem_type << "\n"
                << "material_model = " << material_model << "\n"
                << "formulation = " << input_db->getString("STRESS_FORMULATION") << "\n"
                << "volumetric_mode = " << volumetric_mode << "\n"
                << "proj_type = " << proj_type_str << "\n"
                << "use_default_static_pressure = " << (use_default_static_pressure ? 1 : 0) << "\n"
                << "unit_checks_passed = " << (unit_checks_passed ? 1 : 0) << "\n"
                << "steps = " << n_steps << "\n"
                << std::fixed << std::setprecision(6) << "min_J = " << J_min_global << "\n"
                << "max_J = " << J_max_global << "\n"
                << std::scientific << std::setprecision(6) << "rel_vol_err = " << rel_vol_err << "\n"
                << "top_center_uy = " << top_center_uy << "\n";
            if (volumetric_mode == "STATIC_PRESSURE")
            {
                out << "max_abs_P = " << max_abs_P << "\n";
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
    return 0;
} // main
