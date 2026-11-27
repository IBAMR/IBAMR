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
#include <tbox/Utilities.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_base.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/quadrature.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <boost/multi_array.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Test based on IBFE/explicit/ex12: verifies C^1 barrier smoothness, Nanson
// area pullback identity, 2D cylindrical Hertz contact inversion, and live
// 2D parallel-plate compression against analytic 2D Hertz contact theory.

namespace TestEx12
{
struct NeoHookeanCtx
{
    double G_shear = 0.4666666666666667;
    double kappa_bulk = 23.333333333333332;
    double J_cap = 0.20;
};

struct WallBarrierCtx
{
    double k_rep = 28.0;
    double delta_rep = 0.50;
    double wall_lo_0 = 9.48;
    double wall_hi_0 = 22.52;
    double max_approach = 1.68;
    double t_ramp = 60.0;
    double fillet_frac = 0.12;
};

struct ContactMetrics
{
    double approach = 0.0;
    double F_bot = 0.0;
    double F_top = 0.0;
    double F_avg = 0.0;
    double d_min_bot = 0.0;
    double d_min_top = 0.0;
    double span_y = 0.0;
    double area = 0.0;
};

inline void
compute2DFInvTranspose(const TensorValue<double>& FF, double J_safe, TensorValue<double>& FF_inv_T)
{
    FF_inv_T.zero();
    FF_inv_T(0, 0) = FF(1, 1) / J_safe;
    FF_inv_T(0, 1) = -FF(1, 0) / J_safe;
    FF_inv_T(1, 0) = -FF(0, 1) / J_safe;
    FF_inv_T(1, 1) = FF(0, 0) / J_safe;
    FF_inv_T(2, 2) = 1.0;
}

void
PK1_dev_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const std::vector<const std::vector<double>*>& /*var_data*/,
                        const std::vector<const std::vector<VectorValue<double>>*>& /*grad_var_data*/,
                        double /*time*/,
                        void* ctx)
{
    const auto* c = static_cast<const NeoHookeanCtx*>(ctx);
    PP.zero();
    const double J_raw = FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0);
    const double J_safe = std::max(J_raw, c->J_cap);
    TensorValue<double> FF_inv_T;
    compute2DFInvTranspose(FF, J_safe, FF_inv_T);
    const double I1_3d = FF(0, 0) * FF(0, 0) + FF(0, 1) * FF(0, 1) + FF(1, 0) * FF(1, 0) + FF(1, 1) * FF(1, 1) + 1.0;
    const double coeff = c->G_shear * std::pow(J_safe, -2.0 / 3.0);
    const double mean_I1 = I1_3d / 3.0;
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int A = 0; A < NDIM; ++A)
        {
            PP(i, A) = coeff * (FF(i, A) - mean_I1 * FF_inv_T(i, A));
        }
    }
}

void
PK1_dil_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const std::vector<const std::vector<double>*>& /*var_data*/,
                        const std::vector<const std::vector<VectorValue<double>>*>& /*grad_var_data*/,
                        double /*time*/,
                        void* ctx)
{
    const auto* c = static_cast<const NeoHookeanCtx*>(ctx);
    PP.zero();
    const double J_raw = FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0);
    const double J_safe = std::max(J_raw, c->J_cap);
    TensorValue<double> FF_inv_T;
    compute2DFInvTranspose(FF, J_safe, FF_inv_T);
    const double log_J = (J_raw >= c->J_cap) ? std::log(J_raw) : (std::log(c->J_cap) + (J_raw - c->J_cap) / c->J_cap);
    const double p_dil = c->kappa_bulk * log_J;
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int A = 0; A < NDIM; ++A)
        {
            PP(i, A) = p_dil * FF_inv_T(i, A);
        }
    }
}

inline double
wallApproachDisplacement(double time, const WallBarrierCtx& c)
{
    if (c.max_approach <= 0.0 || time <= 0.0) return 0.0;
    if (c.t_ramp <= 0.0 || time >= c.t_ramp) return c.max_approach;
    const double tau = time / c.t_ramp;
    const double a = std::clamp(c.fillet_frac, 1.0e-4, 0.45);
    double integral = 0.0;
    if (tau < a)
    {
        integral = 0.5 * tau - (0.5 * a / M_PI) * std::sin(M_PI * tau / a);
    }
    else if (tau <= 1.0 - a)
    {
        integral = tau - 0.5 * a;
    }
    else
    {
        const double u = tau - (1.0 - a);
        integral = (1.0 - 1.5 * a) + 0.5 * u + (0.5 * a / M_PI) * std::sin(M_PI * u / a);
    }
    return c.max_approach * std::clamp(integral / (1.0 - a), 0.0, 1.0);
}

inline double
evalRepulsionPressure(double d, double k_rep, double delta_rep)
{
    if (k_rep <= 0.0 || delta_rep <= 0.0 || d >= delta_rep) return 0.0;
    const double d_cap = 0.15 * delta_rep;
    if (d >= d_cap)
    {
        const double s = 1.0 - d / delta_rep;
        return k_rep * s * s * (delta_rep / d);
    }
    const double s_c = 1.0 - d_cap / delta_rep;
    const double ratio_c = delta_rep / d_cap;
    const double p_c = k_rep * s_c * s_c * ratio_c;
    const double slope_mag = (k_rep / delta_rep) * (2.0 * s_c * ratio_c + s_c * s_c * ratio_c * ratio_c);
    return p_c + slope_mag * (d_cap - d);
}

void
wallRepulsionSurfaceForceFcn(VectorValue<double>& F_s,
                             const VectorValue<double>& n,
                             const VectorValue<double>& N,
                             const TensorValue<double>& FF,
                             const libMesh::Point& x,
                             const libMesh::Point& /*X*/,
                             Elem* const /*elem*/,
                             unsigned short int /*side*/,
                             const std::vector<const std::vector<double>*>& /*var_data*/,
                             const std::vector<const std::vector<VectorValue<double>>*>& /*grad_var_data*/,
                             double time,
                             void* ctx)
{
    const auto* c = static_cast<const WallBarrierCtx*>(ctx);
    F_s.zero();
    if (!c || c->k_rep <= 0.0) return;

    const double half_dh = 0.5 * wallApproachDisplacement(time, *c);
    const double d_lo = x(1) - (c->wall_lo_0 + half_dh);
    const double d_hi = (c->wall_hi_0 - half_dh) - x(1);
    const double t_cur_y =
        evalRepulsionPressure(d_lo, c->k_rep, c->delta_rep) - evalRepulsionPressure(d_hi, c->k_rep, c->delta_rep);
    if (t_cur_y == 0.0) return;

    const double J = std::abs(FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0));
    const TensorValue<double> FF_inv_T = IBTK::tensor_inverse_transpose(FF, NDIM);
    const double da_dA0 = J * ((FF_inv_T * N) * n);
    F_s(1) = da_dA0 * t_cur_y;
}

inline double
hertz2DParallelPlatesForce(double R, double G, double kappa, double eps)
{
    if (eps <= 0.0) return 0.0;
    const double E_star = 4.0 * G * kappa / (kappa + G);
    const double dh_target = 2.0 * R * eps;
    auto dh_of_f = [&](double f_val) -> double
    {
        if (f_val <= 0.0) return 0.0;
        const double arg = 4.0 * M_PI * E_star * R / f_val;
        if (arg <= std::exp(1.0)) return 1.0e99;
        return (2.0 * f_val / (M_PI * E_star)) * (std::log(arg) - 1.0);
    };
    const double f_peak = 4.0 * M_PI * E_star * R / std::exp(2.0);
    if (dh_target >= dh_of_f(f_peak)) return std::numeric_limits<double>::quiet_NaN();
    double lo = 0.0;
    double hi = f_peak;
    for (int it = 0; it < 80; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (dh_of_f(mid) < dh_target)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

ContactMetrics
computeContactMetrics(MeshBase& mesh,
                      EquationSystems* equation_systems,
                      const std::string& coords_sys_name,
                      double time,
                      const WallBarrierCtx& c)
{
    System& X_system = equation_systems->get_system<System>(coords_sys_name);
    NumericVector<double>* X_vec = X_system.solution.get();
    NumericVector<double>* X_ghost = X_system.current_local_solution.get();
    IBTK::copy_and_synch(*X_vec, *X_ghost);

    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int>> X_dof_indices(NDIM);

    std::unique_ptr<FEBase> fe_vol(FEBase::build(NDIM, X_dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule_vol = QBase::build(QGAUSS, NDIM, FIFTH);
    fe_vol->attach_quadrature_rule(qrule_vol.get());
    const std::vector<double>& JxW_vol = fe_vol->get_JxW();
    const std::vector<std::vector<RealGradient>>& dphi_vol = fe_vol->get_dphi();

    std::unique_ptr<FEBase> fe_face(FEBase::build(NDIM, X_dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule_face = QBase::build(QGAUSS, NDIM - 1, FIFTH);
    fe_face->attach_quadrature_rule(qrule_face.get());
    const std::vector<double>& JxW_face = fe_face->get_JxW();
    const std::vector<std::vector<Real>>& phi_face = fe_face->get_phi();
    const std::vector<std::vector<RealGradient>>& dphi_face = fe_face->get_dphi();
    const std::vector<libMesh::Point>& normal_face = fe_face->get_normals();

    TensorValue<double> FF, FF_inv_T;
    boost::multi_array<double, 2> X_node;

    const double dh = wallApproachDisplacement(time, c);
    const double w_lo = c.wall_lo_0 + 0.5 * dh;
    const double w_hi = c.wall_hi_0 - 0.5 * dh;

    double local_F_bot = 0.0;
    double local_F_top = 0.0;
    double local_area = 0.0;
    double local_d_bot = 1.0e30;
    double local_d_top = 1.0e30;
    double local_ymin = 1.0e30;
    double local_ymax = -1.0e30;

    for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
    {
        Elem* const elem = *el_it;
        for (unsigned int d = 0; d < NDIM; ++d) X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        IBTK::get_values_for_interpolation(X_node, *X_ghost, X_dof_indices);

        fe_vol->reinit(elem);
        for (unsigned int qp = 0; qp < qrule_vol->n_points(); ++qp)
        {
            IBTK::jacobian(FF, qp, X_node, dphi_vol);
            local_area += std::abs(FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0)) * JxW_vol[qp];
        }

        for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
            if (elem->neighbor_ptr(side) != nullptr) continue;
            fe_face->reinit(elem, side);
            const unsigned int n_qp = qrule_face->n_points();
            const unsigned int n_basis = static_cast<unsigned int>(phi_face.size());
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                IBTK::jacobian(FF, qp, X_node, dphi_face);
                const double J = std::abs(FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0));
                FF_inv_T = IBTK::tensor_inverse_transpose(FF, NDIM);
                const double da = J * (FF_inv_T * normal_face[qp]).norm() * JxW_face[qp];

                double y_qp = 0.0;
                for (unsigned int k = 0; k < n_basis; ++k) y_qp += X_node[k][1] * phi_face[k][qp];
                local_ymin = std::min(local_ymin, y_qp);
                local_ymax = std::max(local_ymax, y_qp);

                const double d_lo = y_qp - w_lo;
                const double d_hi = w_hi - y_qp;
                local_d_bot = std::min(local_d_bot, d_lo);
                local_d_top = std::min(local_d_top, d_hi);
                local_F_bot += evalRepulsionPressure(d_lo, c.k_rep, c.delta_rep) * da;
                local_F_top += evalRepulsionPressure(d_hi, c.k_rep, c.delta_rep) * da;
            }
        }
    }

    ContactMetrics out;
    out.approach = dh;
    out.F_bot = IBTK_MPI::sumReduction(local_F_bot);
    out.F_top = IBTK_MPI::sumReduction(local_F_top);
    out.F_avg = 0.5 * (out.F_bot + out.F_top);
    out.area = IBTK_MPI::sumReduction(local_area);
    out.d_min_bot = IBTK_MPI::minReduction(local_d_bot);
    out.d_min_top = IBTK_MPI::minReduction(local_d_top);
    const double ymin = IBTK_MPI::minReduction(local_ymin);
    const double ymax = IBTK_MPI::maxReduction(local_ymax);
    out.span_y = ymax - ymin;
    return out;
}

void
runAnalyticChecks(double R, double G_shear, double kappa_bulk, double k_rep, double delta_rep)
{
    const double d_cap = 0.15 * delta_rep;
    const double eps_fd = 1.0e-7 * delta_rep;

    // 1. Cutoff at d >= delta_rep and C^1 slope at d_cap
    TBOX_ASSERT(evalRepulsionPressure(delta_rep, k_rep, delta_rep) == 0.0);
    TBOX_ASSERT(evalRepulsionPressure(delta_rep + 1.0e-6, k_rep, delta_rep) == 0.0);
    const double p_cap = evalRepulsionPressure(d_cap, k_rep, delta_rep);
    const double p_left = evalRepulsionPressure(d_cap - eps_fd, k_rep, delta_rep);
    const double p_right = evalRepulsionPressure(d_cap + eps_fd, k_rep, delta_rep);
    const double slope_left = (p_cap - p_left) / eps_fd;
    const double slope_right = (p_right - p_cap) / eps_fd;
    TBOX_ASSERT(std::abs(slope_left - slope_right) / std::abs(slope_left) < 1.0e-5);

    // 2. Nanson pullback identity F_s(1) == (da/dA_0) * p_rep(d) for angled unit N
    WallBarrierCtx ctx;
    ctx.k_rep = k_rep;
    ctx.delta_rep = delta_rep;
    ctx.wall_lo_0 = 0.0;
    ctx.wall_hi_0 = 100.0;
    ctx.max_approach = 0.0;

    TensorValue<double> FF;
    FF.zero();
    FF(0, 0) = 1.4;
    FF(1, 1) = 0.75;
    FF(2, 2) = 1.0;
    const double inv_sqrt5 = 1.0 / std::sqrt(5.0);
    const VectorValue<double> N(2.0 * inv_sqrt5, -1.0 * inv_sqrt5, 0.0);
    const TensorValue<double> FF_inv_T = IBTK::tensor_inverse_transpose(FF, NDIM);
    const VectorValue<double> weighted = FF_inv_T * N;
    const VectorValue<double> n = weighted.unit();
    const double expected_da_dA0 = std::abs(FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0)) * weighted.norm();
    const libMesh::Point x_test(0.0, 0.4 * delta_rep, 0.0);
    VectorValue<double> F_s;
    wallRepulsionSurfaceForceFcn(F_s, n, N, FF, x_test, x_test, nullptr, 0, {}, {}, 0.0, &ctx);
    const double expected_Fs_y = expected_da_dA0 * evalRepulsionPressure(0.4 * delta_rep, k_rep, delta_rep);
    TBOX_ASSERT(std::abs(F_s(1) - expected_Fs_y) < 1.0e-12);
    TBOX_ASSERT(std::abs(F_s(0)) < 1.0e-15);

    // 3. 2D cylindrical Hertz equation residual
    const double eps_ref = 0.12;
    const double f_hz = hertz2DParallelPlatesForce(R, G_shear, kappa_bulk, eps_ref);
    const double E_star = 4.0 * G_shear * kappa_bulk / (kappa_bulk + G_shear);
    const double dh_eval = (2.0 * f_hz / (M_PI * E_star)) * (std::log(4.0 * M_PI * E_star * R / f_hz) - 1.0);
    TBOX_ASSERT(std::abs(dh_eval - 2.0 * R * eps_ref) < 1.0e-12);
}
} // namespace TestEx12

using namespace TestEx12;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double R = input_db->getDoubleWithDefault("CELL_R", 6.0);
        const int r_refine = input_db->getIntegerWithDefault("MESH_REFINE_LEVEL", 3);
        const std::string elem_type = input_db->getStringWithDefault("ELEM_TYPE", "QUAD9");
        const double xc = input_db->getDoubleWithDefault("CELL_X0", 16.0);
        const double yc = input_db->getDoubleWithDefault("CELL_Y0", 16.0);

        NeoHookeanCtx nh_ctx;
        nh_ctx.G_shear = input_db->getDouble("G_SHEAR");
        nh_ctx.kappa_bulk = input_db->getDouble("KAPPA_BULK");

        WallBarrierCtx wall_ctx;
        wall_ctx.k_rep = input_db->getDouble("K_REP");
        wall_ctx.delta_rep = input_db->getDouble("DELTA_REP");
        wall_ctx.wall_lo_0 = input_db->getDouble("WALL_LO");
        wall_ctx.wall_hi_0 = input_db->getDouble("WALL_HI_0");
        wall_ctx.max_approach = input_db->getDouble("MAX_APPROACH");
        wall_ctx.t_ramp = input_db->getDouble("T_RAMP");
        wall_ctx.fillet_frac = input_db->getDouble("FILLET_FRAC");

        runAnalyticChecks(R, nh_ctx.G_shear, nh_ctx.kappa_bulk, wall_ctx.k_rep, wall_ctx.delta_rep);

        Mesh solid_mesh(ibtk_init.getLibMeshInit().comm(), NDIM);
        MeshTools::Generation::build_sphere(solid_mesh, R, r_refine, Utility::string_to_enum<ElemType>(elem_type));

        for (auto el = solid_mesh.elements_begin(); el != solid_mesh.elements_end(); ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor_ptr(side) != nullptr) continue;
                for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                {
                    if (!elem->is_node_on_side(k, side)) continue;
                    Node& n = elem->node_ref(k);
                    n = R * n.unit();
                }
            }
        }
        for (Node* node : solid_mesh.node_ptr_range())
        {
            (*node)(0) += xc;
            (*node)(1) += yc;
        }
        solid_mesh.prepare_for_use();

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &solid_mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           false);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        time_integrator->registerLoadBalancer(load_balancer);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        ib_method_ops->initializeFEEquationSystems();
        const auto dev_order =
            Utility::string_to_enum<Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
        const auto dil_order =
            Utility::string_to_enum<Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "THIRD"));
        ib_method_ops->registerPK1StressFunction(
            IBFEMethod::PK1StressFcnData(PK1_dev_stress_function, {}, &nh_ctx, QGAUSS, dev_order));
        ib_method_ops->registerPK1StressFunction(
            IBFEMethod::PK1StressFcnData(PK1_dil_stress_function, {}, &nh_ctx, QGAUSS, dil_order));
        ib_method_ops->registerLagSurfaceForceFunction(
            IBFEMethod::LagSurfaceForceFcnData(wallRepulsionSurfaceForceFcn, {}, &wall_ctx));

        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
        const std::string coords_name = ib_method_ops->getCurrentCoordinatesSystemName();

        double loop_time = time_integrator->getIntegratorTime();
        const double loop_time_end = time_integrator->getEndTime();
        const ContactMetrics m0 = computeContactMetrics(solid_mesh, equation_systems, coords_name, loop_time, wall_ctx);
        const double span0 = m0.span_y;
        const double area0 = m0.area;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            const double dt = std::min(time_integrator->getMaximumTimeStepSize(), loop_time_end - loop_time);
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
        }

        const ContactMetrics m_end =
            computeContactMetrics(solid_mesh, equation_systems, coords_name, loop_time, wall_ctx);
        const double eps_end = std::max(0.0, 1.0 - m_end.span_y / span0);
        const double f_hz_end = hertz2DParallelPlatesForce(R, nh_ctx.G_shear, nh_ctx.kappa_bulk, eps_end);
        const double f_imb_pct = 100.0 * std::abs(m_end.F_top - m_end.F_bot) / std::max(1.0e-12, m_end.F_avg);
        const double hz_err_pct = 100.0 * std::abs(m_end.F_avg - f_hz_end) / std::max(1.0e-12, f_hz_end);
        const double area_drift_pct = 100.0 * (m_end.area - area0) / area0;

        const double max_hz_err_pct = input_db->getDoubleWithDefault("MAX_HERTZ_REL_ERR_PCT", 10.0);
        TBOX_ASSERT(std::max(m0.F_bot, m0.F_top) == 0.0);
        TBOX_ASSERT(m_end.d_min_bot > 0.0 && m_end.d_min_top > 0.0);
        TBOX_ASSERT(eps_end > 0.08);
        TBOX_ASSERT(f_imb_pct < 1.0);
        TBOX_ASSERT(std::abs(area_drift_pct) < 0.5);
        TBOX_ASSERT(hz_err_pct <= max_hz_err_pct);

        if (IBTK_MPI::getRank() == 0)
        {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(4) << "ANALYTIC_CHECKS PASS\n"
                << "t_end=" << loop_time << " F0_init=" << std::max(m0.F_bot, m0.F_top) << " eps_end=" << eps_end
                << "\nF_bot=" << m_end.F_bot << " F_top=" << m_end.F_top << " F_avg=" << m_end.F_avg
                << " F_hertz_2d=" << f_hz_end << "\nhertz_rel_err_pct=" << hz_err_pct
                << " F_imbalance_pct=" << f_imb_pct << "\nd_min_bot=" << m_end.d_min_bot
                << " d_min_top=" << m_end.d_min_top << " area_drift_pct=" << area_drift_pct;
            pout << oss.str() << "\n";
            std::ofstream out_stream("output");
            out_stream << oss.str() << "\n";
            out_stream.flush();
            TBOX_ASSERT(out_stream.good());
        }
    }

    return 0;
}
