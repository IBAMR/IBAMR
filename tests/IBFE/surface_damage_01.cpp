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

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/MemoryDatabase.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/face_tri6.h>
#include <libmesh/fe_base.h>
#include <libmesh/libmesh.h>
#include <libmesh/node.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>
#include <libmesh/replicated_mesh.h>

#include <petscsys.h>

#include <IBFESurfaceDamage.h>
#include <SAMRAI_config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ibamr/app_namespaces.h>

// Fast 3D CMake attest-style verification for IBAMR::FESurfaceDamage:
//   1. Exact thin-sheet strain-energy & FD PK1 consistency (Barthes-Biesel NH & Skalak).
//   2. Irreversibility, damage initiation threshold, uncapped rupture, and sticky fold freeze.
//   3. Peerlings Helmholtz filter on curved 1280- and 5120-TRI6 spheres (constant preservation,
//      surface integral conservation, Y_{2,0} exact eigenvalue convergence, polar spike stability).
//   4. Collective MPI SAMRAI restart (v2 Allgather/Allgatherv + swapped-ownership prune + v1 load).

namespace
{
using namespace IBAMR::FESurfaceDamage;

std::unique_ptr<libMesh::ReplicatedMesh>
buildSphereTri6Mesh(const libMesh::LibMeshInit& init, const double radius, const int n_subdiv)
{
    auto mesh = std::make_unique<libMesh::ReplicatedMesh>(init.comm(), 2);
    mesh->set_spatial_dimension(3);
    mesh->allow_renumbering(false);

    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    std::vector<std::array<double, 3>> V = {
        { -1, phi, 0 },  { 1, phi, 0 },  { -1, -phi, 0 }, { 1, -phi, 0 }, { 0, -1, phi },  { 0, 1, phi },
        { 0, -1, -phi }, { 0, 1, -phi }, { phi, 0, -1 },  { phi, 0, 1 },  { -phi, 0, -1 }, { -phi, 0, 1 },
    };
    std::vector<std::array<unsigned int, 3>> T = {
        { 0, 11, 5 },  { 0, 5, 1 },  { 0, 1, 7 },  { 0, 7, 10 }, { 0, 10, 11 }, { 1, 5, 9 }, { 5, 11, 4 },
        { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },  { 3, 9, 4 },  { 3, 4, 2 },   { 3, 2, 6 }, { 3, 6, 8 },
        { 3, 8, 9 },   { 4, 9, 5 },  { 2, 4, 11 }, { 6, 2, 10 }, { 8, 6, 7 },   { 9, 8, 1 },
    };

    for (int s = 0; s < n_subdiv; ++s)
    {
        std::map<std::pair<unsigned int, unsigned int>, unsigned int> mc;
        auto gm = [&](unsigned int a, unsigned int b) -> unsigned int
        {
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            auto it = mc.find(key);
            if (it != mc.end()) return it->second;
            const unsigned int id = static_cast<unsigned int>(V.size());
            V.push_back({ 0.5 * (V[a][0] + V[b][0]), 0.5 * (V[a][1] + V[b][1]), 0.5 * (V[a][2] + V[b][2]) });
            mc[key] = id;
            return id;
        };
        std::vector<std::array<unsigned int, 3>> T2;
        T2.reserve(4 * T.size());
        for (const auto& t : T)
        {
            const unsigned int ab = gm(t[0], t[1]);
            const unsigned int bc = gm(t[1], t[2]);
            const unsigned int ca = gm(t[2], t[0]);
            T2.push_back({ t[0], ab, ca });
            T2.push_back({ t[1], bc, ab });
            T2.push_back({ t[2], ca, bc });
            T2.push_back({ ab, bc, ca });
        }
        T.swap(T2);
    }

    auto mapPoint = [&](const std::array<double, 3>& p) -> libMesh::Point
    {
        const double rn = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        return libMesh::Point(radius * p[0] / rn, radius * p[1] / rn, radius * p[2] / rn);
    };

    for (std::size_t i = 0; i < V.size(); ++i)
    {
        mesh->add_point(mapPoint(V[i]), static_cast<libMesh::dof_id_type>(i));
    }

    std::map<std::pair<unsigned int, unsigned int>, libMesh::dof_id_type> edge_mid;
    auto getMid = [&](unsigned int a, unsigned int b) -> libMesh::dof_id_type
    {
        auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = edge_mid.find(key);
        if (it != edge_mid.end()) return it->second;
        const libMesh::dof_id_type mid_id = static_cast<libMesh::dof_id_type>(mesh->n_nodes());
        const std::array<double, 3> avg = { 0.5 * (V[a][0] + V[b][0]),
                                            0.5 * (V[a][1] + V[b][1]),
                                            0.5 * (V[a][2] + V[b][2]) };
        mesh->add_point(mapPoint(avg), mid_id);
        edge_mid[key] = mid_id;
        return mid_id;
    };

    for (std::size_t ti = 0; ti < T.size(); ++ti)
    {
        unsigned int a = T[ti][0], b = T[ti][1], c = T[ti][2];
        const libMesh::Point& pa = mesh->point(a);
        const libMesh::Point& pb = mesh->point(b);
        const libMesh::Point& pc = mesh->point(c);
        if ((pb - pa).cross(pc - pa) * ((pa + pb + pc) * (1.0 / 3.0)) < 0.0)
        {
            std::swap(b, c);
        }
        const libMesh::dof_id_type mab = getMid(a, b);
        const libMesh::dof_id_type mbc = getMid(b, c);
        const libMesh::dof_id_type mca = getMid(c, a);
        auto* tri = new libMesh::Tri6();
        tri->set_id() = static_cast<libMesh::dof_id_type>(ti);
        const std::array<libMesh::dof_id_type, 6> conn = { a, b, c, mab, mbc, mca };
        for (unsigned int k = 0; k < 6; ++k) tri->set_node(k) = mesh->node_ptr(conn[k]);
        mesh->add_elem(tri);
    }
    mesh->prepare_for_use();
    return mesh;
}

struct HelmholtzMetrics
{
    double const_max_err = 0.0;
    double y20_int_rel_err = 0.0;
    double y20_linf_err = 0.0;
    double spike_ybar_peak = 0.0;
};

HelmholtzMetrics
runHelmholtzSuiteOnMesh(const libMesh::LibMeshInit& init, const double radius, const int n_subdiv, const double l_d)
{
    auto mesh = buildSphereTri6Mesh(init, radius, n_subdiv);
    TBOX_ASSERT(IBTK::IBTK_MPI::minReduction(static_cast<int>(mesh->n_active_local_elem())) > 0);

    const libMesh::FEType fe_type(libMesh::SECOND, libMesh::LAGRANGE);
    SurfaceNonlocalHelmholtzFilter filter;
    filter.ensureInitialized(*mesh, fe_type, libMesh::QGAUSS, libMesh::FIFTH, l_d);

    std::unique_ptr<libMesh::QBase> qrule = libMesh::QBase::build(libMesh::QGAUSS, 2, libMesh::FIFTH);
    std::unique_ptr<libMesh::FEBase> fe = libMesh::FEBase::build(2, fe_type);
    fe->attach_quadrature_rule(qrule.get());
    const auto& xyz_qp = fe->get_xyz();
    const auto& dX_dxi = fe->get_dxyzdxi();
    const auto& dX_deta = fe->get_dxyzdeta();
    const auto& qw = qrule->get_weights();

    const std::size_t n_loc = static_cast<std::size_t>(mesh->n_active_local_elem());
    std::vector<std::vector<double>> y_const(n_loc), y_y20(n_loc), y_spike(n_loc);
    std::vector<std::vector<double>> y20_exact(n_loc), jxw0(n_loc);

    const double eig_factor = 1.0 / (1.0 + 6.0 * (l_d * l_d) / (radius * radius));
    const double sigma = 0.70;
    const double y_bg = 0.10;
    const double y_amp = 3.0;

    std::size_t idx = 0;
    for (auto el_it = mesh->active_local_elements_begin(); el_it != mesh->active_local_elements_end(); ++el_it, ++idx)
    {
        libMesh::Elem* const elem = *el_it;
        fe->reinit(elem);
        const unsigned int n_qp = qrule->n_points();
        y_const[idx].assign(n_qp, 3.25);
        y_y20[idx].resize(n_qp);
        y_spike[idx].resize(n_qp);
        y20_exact[idx].resize(n_qp);
        jxw0[idx].resize(n_qp);

        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            SurfaceTangents G0;
            for (unsigned int i = 0; i < 3; ++i)
            {
                G0(i, 0) = dX_dxi[qp](i);
                G0(i, 1) = dX_deta[qp](i);
            }
            const Metric2x2 M0 = computeSurfaceMetric(G0);
            jxw0[idx][qp] = std::sqrt(M0.det) * qw[qp];

            const libMesh::Point& p = xyz_qp[qp];
            const double uz = p(2) / radius;
            const double y20 = 0.5 * (3.0 * uz * uz - 1.0);
            y_y20[idx][qp] = 2.0 + y20;
            y20_exact[idx][qp] = 2.0 + eig_factor * y20;

            const double r2_xy = p(0) * p(0) + p(1) * p(1);
            y_spike[idx][qp] = y_bg + ((p(2) > 0.0) ? (y_amp * std::exp(-r2_xy / (2.0 * sigma * sigma))) : 0.0);
        }
    }

    HelmholtzMetrics hm;
    std::vector<std::vector<double>> yb_const, yb_y20, yb_spike;

    filter.filterQuadratureField(y_const, yb_const);
    double loc_const_err = 0.0;
    for (const auto& vec : yb_const)
    {
        for (double v : vec) loc_const_err = std::max(loc_const_err, std::abs(v - 3.25));
    }
    hm.const_max_err = IBTK::IBTK_MPI::maxReduction(loc_const_err);

    filter.filterQuadratureField(y_y20, yb_y20);
    double loc_int_in = 0.0, loc_int_out = 0.0, loc_linf = 0.0;
    for (std::size_t ei = 0; ei < n_loc; ++ei)
    {
        for (std::size_t qp = 0; qp < yb_y20[ei].size(); ++qp)
        {
            loc_int_in += y_y20[ei][qp] * jxw0[ei][qp];
            loc_int_out += yb_y20[ei][qp] * jxw0[ei][qp];
            loc_linf = std::max(loc_linf, std::abs(yb_y20[ei][qp] - y20_exact[ei][qp]));
        }
    }
    const double glob_in = IBTK::IBTK_MPI::sumReduction(loc_int_in);
    const double glob_out = IBTK::IBTK_MPI::sumReduction(loc_int_out);
    hm.y20_int_rel_err = std::abs(glob_out - glob_in) / glob_in;
    hm.y20_linf_err = IBTK::IBTK_MPI::maxReduction(loc_linf);

    filter.filterQuadratureField(y_spike, yb_spike);
    double loc_peak = 0.0;
    for (const auto& vec : yb_spike)
    {
        for (double v : vec) loc_peak = std::max(loc_peak, v);
    }
    hm.spike_ybar_peak = IBTK::IBTK_MPI::maxReduction(loc_peak);
    return hm;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "surface_damage_01.log");

        // 1. Analytic thin-sheet strain energy & finite-difference PK1 patch check
        SurfaceTangents G_0;
        G_0(0, 0) = 1.0;
        G_0(1, 1) = 1.0;
        SurfaceTangents G_t;
        const double l1 = 1.35;
        const double l2 = 1.15;
        G_t(0, 0) = l1;
        G_t(0, 1) = 0.25;
        G_t(1, 1) = l2;

        const double G_s = 0.50;
        const double K_s = 3.0 * G_s;
        SurfaceDamageParams dp_nh;
        dp_nh.enabled = true;
        dp_nh.law = SurfaceDamageLawType::BARTHES_BIESEL_NEO_HOOKEAN;
        dp_nh.G_s = G_s;
        dp_nh.T_0 = 0.20;
        dp_nh.Y_D = 0.20;
        dp_nh.Y_C = 2.0;

        // Exact equibiaxial stretch check: G_eq = diag(lam, lam)
        SurfaceTangents G_eq;
        const double lam = 1.30;
        G_eq(0, 0) = lam;
        G_eq(1, 1) = lam;
        const double J_eq = lam * lam;
        const double I1_eq = 2.0 * lam * lam - 2.0;
        const double I2_eq = J_eq * J_eq - 1.0;
        const double W_nh_exact = 0.5 * G_s * (I1_eq + 1.0 / (J_eq * J_eq) - 1.0);
        const double nh_energy_err = std::abs(evalUndamagedSurfaceStrainEnergy(dp_nh, G_0, G_eq) - W_nh_exact);
        TBOX_ASSERT(nh_energy_err < 1.0e-14);

        SurfaceDamageParams dp_sk = dp_nh;
        dp_sk.law = SurfaceDamageLawType::SKALAK;
        dp_sk.K_s = K_s;
        const double W_sk_exact =
            0.25 * G_s * (I1_eq * I1_eq + 2.0 * I1_eq - 2.0 * I2_eq) + 0.125 * (K_s - G_s) * (I2_eq * I2_eq);
        const double sk_energy_err = std::abs(evalUndamagedSurfaceStrainEnergy(dp_sk, G_0, G_eq) - W_sk_exact);
        TBOX_ASSERT(sk_energy_err < 1.0e-13);

        // Central finite-difference check of damaged P_tilde = d((1-d)*W_0 + T_0 * J_s) / dG_t
        const double d_test = 0.35;
        const double mult = 1.0 - d_test;
        SurfaceTangents P_ana;
        evalDamagedSurfacePK1(dp_nh, G_0, G_t, P_ana, mult);
        const double eps = 1.0e-6;
        double pk1_fd_err = 0.0;
        auto totalDegradedEnergy = [&](const SurfaceTangents& G) -> double
        {
            double js = 1.0;
            const double W0 = evalUndamagedSurfaceStrainEnergy(dp_nh, G_0, G, &js);
            return mult * W0 + dp_nh.T_0 * js;
        };
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int a = 0; a < 2; ++a)
            {
                SurfaceTangents Gp = G_t, Gm = G_t;
                Gp(i, a) += eps;
                Gm(i, a) -= eps;
                const double fd = (totalDegradedEnergy(Gp) - totalDegradedEnergy(Gm)) / (2.0 * eps);
                pk1_fd_err = std::max(pk1_fd_err, std::abs(P_ana(i, a) - fd));
            }
        }
        TBOX_ASSERT(pk1_fd_err < 1.0e-7);

        // 2. Irreversibility, below-threshold d=0, uncapped rupture d=1, and sticky fold mask
        SurfaceDamageTracker tr_local;
        tr_local.configure(dp_nh);
        const auto s0 = tr_local.updateSlotFromEnergy(1U, 0U, 0.15 * G_s, false);
        const auto s1 = tr_local.updateSlotFromEnergy(1U, 0U, 0.80 * G_s, false);
        const auto s2 = tr_local.updateSlotFromEnergy(1U, 0U, 0.30 * G_s, false);
        const auto s_rup = tr_local.updateSlotFromEnergy(2U, 0U, 2.40 * G_s, false);
        tr_local.updateSlotFromEnergy(3U, 0U, 0.80 * G_s, false);
        const auto s_fold = tr_local.updateSlotFromEnergy(3U, 0U, 2.50 * G_s, true);
        TBOX_ASSERT(s0.d_eff == 0.0);
        TBOX_ASSERT(std::abs(s1.d_eff - 0.30) < 1.0e-14 && std::abs(s2.d_eff - 0.30) < 1.0e-14);
        TBOX_ASSERT(s_rup.ruptured && s_rup.d_eff == 1.0);
        TBOX_ASSERT(s_fold.fold_invalid && std::abs(s_fold.d_eff - 0.30) < 1.0e-14);

        // 3. Peerlings Helmholtz filter convergence on 1280 and 5120 TRI6 spheres (R = 6.0, l_d = 1.0)
        const HelmholtzMetrics hm1 = runHelmholtzSuiteOnMesh(init, 6.0, 3 /* 1280 TRI6 */, 1.0);
        const HelmholtzMetrics hm2 = runHelmholtzSuiteOnMesh(init, 6.0, 4 /* 5120 TRI6 */, 1.0);
        const double helm_const_err = std::max(hm1.const_max_err, hm2.const_max_err);
        const double helm_int_err = std::max(hm1.y20_int_rel_err, hm2.y20_int_rel_err);
        const double helm_spike_rel = std::abs(hm2.spike_ybar_peak - hm1.spike_ybar_peak) / hm2.spike_ybar_peak;
        TBOX_ASSERT(helm_const_err <= 1.0e-11);
        TBOX_ASSERT(helm_int_err <= 1.0e-11);
        TBOX_ASSERT(hm1.y20_linf_err <= 5.0e-3);
        TBOX_ASSERT(hm2.y20_linf_err < hm1.y20_linf_err);
        TBOX_ASSERT(helm_spike_rel <= 0.02);

        // 4. Collective MPI v2 SAMRAI restart gather, deduplicate, swapped prune, & v1 compatibility
        const int rank = IBTK_MPI::getRank();
        const int nprocs = IBTK_MPI::getNodes();
        SurfaceDamageTracker mpi_tr;
        mpi_tr.configure(dp_nh);
        if (nprocs == 1 || rank == 0)
        {
            mpi_tr.updateSlotFromEnergy(10U, 2U, 0.80 * G_s, false);
            mpi_tr.updateSlotFromEnergy(55U, 0U, 0.30 * G_s, false);
            mpi_tr.setCommittedTime(10.0);
        }
        if (nprocs == 1 || rank == nprocs - 1)
        {
            mpi_tr.updateSlotFromEnergy(77U, 5U, 2.40 * G_s, false);
            mpi_tr.updateSlotFromEnergy(88U, 1U, 0.0, true);
            mpi_tr.updateSlotFromEnergy(55U, 0U, 0.60 * G_s, false);
            mpi_tr.setCommittedTime(12.5);
        }

        Pointer<MemoryDatabase> db = new MemoryDatabase("dmg_mpi_restart");
        mpi_tr.putToDatabase(db);
        TBOX_ASSERT(db->getInteger("damage_tracker_version") == 2);
        TBOX_ASSERT(db->getInteger("damage_n_slots") == 4);
        const int global_slots = db->getInteger("damage_n_slots");

        SurfaceDamageTracker mpi_rest;
        mpi_rest.configure(dp_nh);
        mpi_rest.getFromDatabase(db);
        std::unordered_set<std::uint64_t> owned_after;
        if (nprocs == 1)
        {
            owned_after = { 10U, 55U, 77U, 88U };
        }
        else if (rank == 0)
        {
            owned_after = { 55U, 77U };
        }
        else if (rank == 1)
        {
            owned_after = { 10U, 88U };
        }
        mpi_rest.pruneToLocalElements(owned_after);
        const auto rest_loc = mpi_rest.summarizeLocal();
        TBOX_ASSERT(IBTK_MPI::sumReduction(rest_loc.n_total_qp) == 4);
        TBOX_ASSERT(IBTK_MPI::sumReduction(rest_loc.n_ruptured_qp) == 1);
        TBOX_ASSERT(IBTK_MPI::sumReduction(rest_loc.n_fold_invalid_qp) == 1);
        const double restart_mpi_diff = IBTK_MPI::maxReduction(
            std::abs(mpi_rest.getSlot(55U, 0U).W0_qp_max - (owned_after.count(55U) ? 0.60 * G_s : 0.0)));
        TBOX_ASSERT(restart_mpi_diff <= 1.0e-14);

        Pointer<MemoryDatabase> v1_db = new MemoryDatabase("dmg_v1_compat");
        v1_db->putInteger("damage_tracker_version", 1);
        v1_db->putDouble("damage_committed_time", 4.5);
        v1_db->putInteger("damage_n_slots", 1);
        double kh = 0.0, kl = static_cast<double>(SurfaceDamageTracker::makeKey(5U, 3U));
        double w0 = 0.4, du = 0.3, de = 0.3, pm = 0.0;
        int rp = 0, fi = 0;
        v1_db->putDoubleArray("damage_keys_hi", &kh, 1);
        v1_db->putDoubleArray("damage_keys_lo", &kl, 1);
        v1_db->putDoubleArray("damage_W0_max", &w0, 1);
        v1_db->putDoubleArray("damage_d_uncapped", &du, 1);
        v1_db->putDoubleArray("damage_d_eff", &de, 1);
        v1_db->putDoubleArray("damage_P_m", &pm, 1);
        v1_db->putIntegerArray("damage_ruptured", &rp, 1);
        v1_db->putIntegerArray("damage_fold_invalid", &fi, 1);
        SurfaceDamageTracker rest_v1;
        rest_v1.configure(dp_nh);
        rest_v1.getFromDatabase(v1_db);
        TBOX_ASSERT(rest_v1.committedTime() == 4.5 && rest_v1.getSlot(5U, 3U).W0_qp_max == 0.4);

        if (rank == 0)
        {
            const auto fmt_tol = [](const double val) { return (val < 1.0e-12) ? 0.0 : val; };
            std::ostringstream oss;
            oss << "SURFACE_DAMAGE_CHECKS PASS\n"
                << std::scientific << std::setprecision(5) << "nh_energy_err=" << fmt_tol(nh_energy_err) << "\n"
                << "sk_energy_err=" << fmt_tol(sk_energy_err) << "\n"
                << "pk1_fd_lt_1e7=" << (pk1_fd_err < 1.0e-7 ? 1 : 0) << "\n"
                << "helm_const_err=" << fmt_tol(helm_const_err) << "\n"
                << "helm_int_err=" << fmt_tol(helm_int_err) << "\n"
                << "helm_y20_err_1x=" << hm1.y20_linf_err << "\n"
                << "helm_y20_err_2x=" << hm2.y20_linf_err << "\n"
                << "helm_spike_rel=" << std::fixed << std::setprecision(5) << helm_spike_rel << "\n"
                << std::scientific << std::setprecision(5) << "restart_mpi_diff=" << fmt_tol(restart_mpi_diff) << "\n"
                << "global_slots=" << global_slots << "\n";
            pout << oss.str();
            std::ofstream out("output");
            out << oss.str();
            out.flush();
            TBOX_ASSERT(out.good());
        }
    }
    return 0;
}
