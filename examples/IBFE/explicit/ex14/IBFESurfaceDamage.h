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

// Irreversible Marigo scalar continuum surface damage model paired with a
// Peerlings implicit-gradient (Helmholtz) nonlocal damage energy filter on
// curved codimension-1 TRI3 / TRI6 reference surfaces (`Gamma_0`):
//
//   (I - l_d^2 * Delta_{s,0}) Y_bar = Y   on Gamma_0,
//   d(X, t) = clip((max_{0 <= tau <= t} Y_bar(X, tau) / G_s - Y_D) / Y_C, 0, 1),
//   d_eff   = min(d, d_cap),
//   P_tilde_elastic_degraded = (1 - d_eff) * P_tilde_elastic + P_tilde_prestress.
//
// Constant surface Helmholtz operator (M_0 + l_d^2 K_0) is assembled once on
// mesh.comm() via PETSc MATAIJ + KSPCG + PCJACOBI without altering libMesh
// EquationSystems schemas, and restart state is gathered collectively across
// MPI ranks via sorted-key SAMRAI databases.
//
// Sources:
//   - Marigo, J.-J. (1981), C. R. Acad. Sci. Paris 292(19):1309-1312.
//   - Peerlings et al. (1996), Int. J. Numer. Meth. Engng. 39:3391-3403.
//   - Barthes-Biesel, Diaz & Dhenin (2002), J. Fluid Mech. 460:211-222.
//   - Grandmaison, Brancherie & Salsac (2021), J. Fluid Mech. 914:A25.
//   - Lu, R., Yu, P. & Sui, Y. (2024), Soft Matter 20:4057-4071.

#ifndef included_IBAMR_IBFESurfaceDamage
#define included_IBAMR_IBFESurfaceDamage

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <libmesh/elem.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe_base.h>
#include <libmesh/fe_type.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include <libmesh/quadrature.h>
#include <libmesh/tensor_value.h>

#include <petscis.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace IBAMR
{
namespace FESurfaceDamage
{
/*!
 * 3x2 matrix holding covariant surface tangent vectors G = [dx/dxi^1, dx/dxi^2]
 * or parametric first Piola-Kirchhoff surface tractions P_tilde = dW/dG_t.
 */
struct SurfaceTangents
{
    double m[3][2] = { { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } };

    inline double& operator()(unsigned int i, unsigned int a)
    {
        return m[i][a];
    }

    inline double operator()(unsigned int i, unsigned int a) const
    {
        return m[i][a];
    }

    inline void setZero()
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            m[i][0] = 0.0;
            m[i][1] = 0.0;
        }
    }
};

struct Metric2x2
{
    double m00 = 0.0;
    double m01 = 0.0;
    double m11 = 0.0;
    double det = 0.0;
};

inline Metric2x2
computeSurfaceMetric(const SurfaceTangents& G)
{
    Metric2x2 M;
    for (unsigned int i = 0; i < 3; ++i)
    {
        M.m00 += G(i, 0) * G(i, 0);
        M.m01 += G(i, 0) * G(i, 1);
        M.m11 += G(i, 1) * G(i, 1);
    }
    M.det = M.m00 * M.m11 - M.m01 * M.m01;
    return M;
}

inline SurfaceTangents
multiplyTangentByInvMetric(const SurfaceTangents& G_t, const Metric2x2& M)
{
    const double inv00 = M.m11 / M.det;
    const double inv01 = -M.m01 / M.det;
    const double inv11 = M.m00 / M.det;
    SurfaceTangents out;
    for (unsigned int i = 0; i < 3; ++i)
    {
        out(i, 0) = G_t(i, 0) * inv00 + G_t(i, 1) * inv01;
        out(i, 1) = G_t(i, 0) * inv01 + G_t(i, 1) * inv11;
    }
    return out;
}

inline double
traceInverseMetricProduct(const Metric2x2& M_0, const Metric2x2& M_t)
{
    return (M_0.m11 * M_t.m00 - 2.0 * M_0.m01 * M_t.m01 + M_0.m00 * M_t.m11) / M_0.det;
}

inline void
computeSurfaceDeformationGradient(const SurfaceTangents& G_0,
                                  const SurfaceTangents& G_t,
                                  libMesh::TensorValue<double>& F_s)
{
    const Metric2x2 M_0 = computeSurfaceMetric(G_0);
    F_s.zero();
    if (!(std::isfinite(M_0.det) && M_0.det > 1.0e-30)) return;
    const SurfaceTangents Gt_M0_inv = multiplyTangentByInvMetric(G_t, M_0);
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            F_s(i, j) = Gt_M0_inv(i, 0) * G_0(j, 0) + Gt_M0_inv(i, 1) * G_0(j, 1);
        }
    }
}

enum class SurfaceDamageLawType
{
    BARTHES_BIESEL_NEO_HOOKEAN = 0,
    SKALAK = 1
};

struct SurfaceDamageParams
{
    bool enabled = false;
    bool two_way_coupling = true;
    SurfaceDamageLawType law = SurfaceDamageLawType::BARTHES_BIESEL_NEO_HOOKEAN;
    double G_s = 0.0;             // Surface shear modulus [F / L]
    double K_s = 0.0;             // Area-dilation modulus [F / L] (SKALAK: K_s = 3*G_s -> C = 1)
    double T_0 = 0.0;             // Undegraded isotropic prestress [F / L]
    double Y_D = 0.20;            // Dimensionless initiation threshold (W_0 / G_s > Y_D)
    double Y_C = 2.0;             // Dimensionless linear hardening slope Y_C > 0
    double d_cap = 1.0;           // 1.0 = uncapped damage evolution
    double k_p = 0.0;             // Solute permeability coefficient P_m = k_p * d [L / T]
    double nonlocal_length = 0.0; // Peerlings Helmholtz length l_d [L] (<= 0 disables filter)
    double Js_floor = 1.0e-4;     // Positive lower floor on surface area stretch J_s
};

inline double
evalUndamagedSurfaceStrainEnergy(const SurfaceDamageParams& p,
                                 const SurfaceTangents& G_0,
                                 const SurfaceTangents& G_t,
                                 double* J_s_raw_out = nullptr)
{
    TBOX_ASSERT(!p.enabled || (p.G_s > 0.0 && p.Y_C > 0.0));
    const Metric2x2 M_0 = computeSurfaceMetric(G_0);
    const Metric2x2 M_t = computeSurfaceMetric(G_t);
    if (!(std::isfinite(M_0.det) && std::isfinite(M_t.det) && M_0.det > 1.0e-30 && M_t.det >= 0.0))
    {
        if (J_s_raw_out) *J_s_raw_out = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double dA_0 = std::sqrt(M_0.det);
    const double da = std::sqrt(M_t.det);
    const double J_s_raw = da / dA_0;
    if (J_s_raw_out) *J_s_raw_out = J_s_raw;
    const double J_s = std::max(p.Js_floor, J_s_raw);

    const double tr_C = traceInverseMetricProduct(M_0, M_t);
    const double I_1 = tr_C - 2.0;

    if (p.law == SurfaceDamageLawType::SKALAK)
    {
        const double I_2 = J_s * J_s - 1.0;
        const double W_sk = 0.25 * p.G_s * (I_1 * I_1 + 2.0 * I_1 - 2.0 * I_2) + 0.125 * (p.K_s - p.G_s) * (I_2 * I_2);
        return std::max(0.0, W_sk);
    }

    const double W_nh = 0.5 * p.G_s * (I_1 + 1.0 / (J_s * J_s) - 1.0);
    return std::max(0.0, W_nh);
}

inline void
evalDamagedSurfacePK1(const SurfaceDamageParams& p,
                      const SurfaceTangents& G_0,
                      const SurfaceTangents& G_t,
                      SurfaceTangents& P_tilde_out,
                      double damage_elastic_mult = 1.0,
                      double* J_s_raw_out = nullptr)
{
    P_tilde_out.setZero();
    const Metric2x2 M_0 = computeSurfaceMetric(G_0);
    const Metric2x2 M_t = computeSurfaceMetric(G_t);
    if (!(std::isfinite(M_0.det) && std::isfinite(M_t.det) && M_0.det > 1.0e-30 && M_t.det > 1.0e-30))
    {
        if (J_s_raw_out) *J_s_raw_out = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const double J_s_raw = std::sqrt(M_t.det) / std::sqrt(M_0.det);
    if (J_s_raw_out) *J_s_raw_out = J_s_raw;
    const double J_s = std::max(p.Js_floor, J_s_raw);

    const SurfaceTangents Gt_M0_inv = multiplyTangentByInvMetric(G_t, M_0);
    const SurfaceTangents Gt_Mt_inv = multiplyTangentByInvMetric(G_t, M_t);
    const double tr_C = traceInverseMetricProduct(M_0, M_t);
    const double I_1 = tr_C - 2.0;

    double coef_ref = 0.0;
    double coef_elastic_cur = 0.0;
    if (p.law == SurfaceDamageLawType::SKALAK)
    {
        const double I_2 = J_s * J_s - 1.0;
        coef_ref = p.G_s * (I_1 + 1.0);
        coef_elastic_cur = (J_s * J_s) * (0.5 * (p.K_s - p.G_s) * I_2 - p.G_s);
    }
    else
    {
        coef_ref = p.G_s;
        coef_elastic_cur = -p.G_s / (J_s * J_s);
    }

    const double d_mult = std::clamp(damage_elastic_mult, 0.0, 1.0);
    const double coef_prestress_cur = p.T_0 * J_s;
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int a = 0; a < 2; ++a)
        {
            P_tilde_out(i, a) = d_mult * (coef_ref * Gt_M0_inv(i, a) + coef_elastic_cur * Gt_Mt_inv(i, a)) +
                                coef_prestress_cur * Gt_Mt_inv(i, a);
        }
    }
}

inline double
evalMarigoDamage(const SurfaceDamageParams& p, double W0_max)
{
    if (!p.enabled || !std::isfinite(W0_max)) return 0.0;
    TBOX_ASSERT(p.G_s > 0.0);
    TBOX_ASSERT(p.Y_C > 0.0);
    const double y_norm = W0_max / p.G_s;
    if (y_norm <= p.Y_D) return 0.0;
    const double d = (y_norm - p.Y_D) / p.Y_C;
    return std::clamp(d, 0.0, 1.0);
}

struct SurfaceDamageQPSlot
{
    double W0_qp_max = 0.0;
    double d_uncapped = 0.0;
    double d_eff = 0.0;
    double P_m = 0.0;
    bool ruptured = false;
    bool fold_invalid = false;
};

struct SurfaceDamageSummary
{
    double d_uncapped_max = 0.0;
    double d_eff_max = 0.0;
    double y_norm_max = 0.0;
    double P_m_max = 0.0;
    int n_ruptured_qp = 0;
    int n_fold_invalid_qp = 0;
    int n_damaged_qp = 0;
    int n_total_qp = 0;
};

class SurfaceDamageTracker
{
public:
    SurfaceDamageTracker() = default;

    void configure(const SurfaceDamageParams& params)
    {
        TBOX_ASSERT(!params.enabled || (params.G_s > 0.0 && params.Y_C > 0.0));
        d_params = params;
    }

    const SurfaceDamageParams& params() const
    {
        return d_params;
    }

    void clear()
    {
        d_slots.clear();
        d_committed_time = 0.0;
    }

    static std::uint64_t makeKey(std::uint64_t elem_id, unsigned int qp)
    {
        return (elem_id << 8) | static_cast<std::uint64_t>(qp & 0xFFu);
    }

    double getElasticMultiplier(std::uint64_t elem_id, unsigned int qp) const
    {
        if (!d_params.enabled || !d_params.two_way_coupling) return 1.0;
        const auto it = d_slots.find(makeKey(elem_id, qp));
        if (it == d_slots.end()) return 1.0;
        return std::max(0.0, 1.0 - it->second.d_eff);
    }

    SurfaceDamageQPSlot getSlot(std::uint64_t elem_id, unsigned int qp) const
    {
        const auto it = d_slots.find(makeKey(elem_id, qp));
        return (it != d_slots.end()) ? it->second : SurfaceDamageQPSlot{};
    }

    SurfaceDamageQPSlot
    updateSlotFromEnergy(std::uint64_t elem_id, unsigned int qp, double W0_effective, bool elem_fold_invalid)
    {
        SurfaceDamageQPSlot& slot = d_slots[makeKey(elem_id, qp)];
        if (elem_fold_invalid || !std::isfinite(W0_effective))
        {
            slot.fold_invalid = true;
            return slot;
        }
        if (slot.fold_invalid)
        {
            return slot;
        }
        const double w_eff = std::max(0.0, W0_effective);
        if (w_eff > slot.W0_qp_max)
        {
            slot.W0_qp_max = w_eff;
        }
        slot.d_uncapped = evalMarigoDamage(d_params, slot.W0_qp_max);
        slot.d_eff = std::min(slot.d_uncapped, std::clamp(d_params.d_cap, 0.0, 1.0));
        slot.P_m = (d_params.k_p > 0.0) ? (d_params.k_p * slot.d_uncapped) : 0.0;
        slot.ruptured = d_params.enabled && (d_params.G_s > 0.0) &&
                        ((slot.W0_qp_max / d_params.G_s) >= (d_params.Y_D + d_params.Y_C));
        return slot;
    }

    SurfaceDamageQPSlot updateSlot(std::uint64_t elem_id,
                                   unsigned int qp,
                                   const SurfaceTangents& G_0,
                                   const SurfaceTangents& G_t,
                                   bool elem_fold_invalid)
    {
        if (elem_fold_invalid)
        {
            return updateSlotFromEnergy(elem_id, qp, 0.0, true);
        }
        double js_raw = 1.0;
        const double W0 = evalUndamagedSurfaceStrainEnergy(d_params, G_0, G_t, &js_raw);
        const bool invalid = !std::isfinite(W0) || !std::isfinite(js_raw) || js_raw <= 0.0;
        return updateSlotFromEnergy(elem_id, qp, W0, invalid);
    }

    std::size_t pruneToLocalElements(const std::unordered_set<std::uint64_t>& local_elem_ids)
    {
        std::size_t removed = 0;
        for (auto it = d_slots.begin(); it != d_slots.end();)
        {
            const std::uint64_t elem_id = it->first >> 8;
            if (local_elem_ids.find(elem_id) == local_elem_ids.end())
            {
                it = d_slots.erase(it);
                ++removed;
            }
            else
            {
                ++it;
            }
        }
        return removed;
    }

    void setCommittedTime(double t)
    {
        d_committed_time = t;
    }

    double committedTime() const
    {
        return d_committed_time;
    }

    SurfaceDamageSummary summarizeLocal() const
    {
        SurfaceDamageSummary s;
        for (const auto& kv : d_slots)
        {
            const SurfaceDamageQPSlot& slot = kv.second;
            ++s.n_total_qp;
            if (slot.fold_invalid) ++s.n_fold_invalid_qp;
            if (slot.ruptured) ++s.n_ruptured_qp;
            if (slot.d_uncapped > 0.0) ++s.n_damaged_qp;
            s.d_uncapped_max = std::max(s.d_uncapped_max, slot.d_uncapped);
            s.d_eff_max = std::max(s.d_eff_max, slot.d_eff);
            if (d_params.G_s > 0.0)
            {
                s.y_norm_max = std::max(s.y_norm_max, slot.W0_qp_max / d_params.G_s);
            }
            s.P_m_max = std::max(s.P_m_max, slot.P_m);
        }
        return s;
    }

    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) const
    {
        TBOX_ASSERT(db);

        struct PackedSlot
        {
            std::uint64_t key = 0;
            SurfaceDamageQPSlot slot;
        };
        std::vector<PackedSlot> packed;
        packed.reserve(d_slots.size());
        for (const auto& kv : d_slots)
        {
            packed.push_back({ kv.first, kv.second });
        }

        double committed_time_out = d_committed_time;
        int mpi_init = 0;
        int mpi_fini = 0;
        MPI_Initialized(&mpi_init);
        MPI_Finalized(&mpi_fini);
        if (mpi_init && !mpi_fini)
        {
            const MPI_Comm comm = IBTK::IBTK_MPI::getCommunicator();
            int nprocs = 1;
            MPI_Comm_size(comm, &nprocs);
            if (nprocs > 1)
            {
                MPI_Allreduce(&d_committed_time, &committed_time_out, 1, MPI_DOUBLE, MPI_MAX, comm);

                const int local_n = static_cast<int>(packed.size());
                std::vector<int> counts(static_cast<std::size_t>(nprocs), 0);
                MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

                std::vector<int> displs(static_cast<std::size_t>(nprocs), 0);
                int total_n = 0;
                for (int r = 0; r < nprocs; ++r)
                {
                    displs[r] = total_n;
                    total_n += counts[r];
                }

                std::vector<unsigned long long> loc_keys(static_cast<std::size_t>(local_n));
                std::vector<double> loc_vals(static_cast<std::size_t>(6 * local_n));
                for (int i = 0; i < local_n; ++i)
                {
                    loc_keys[i] = static_cast<unsigned long long>(packed[i].key);
                    loc_vals[6 * i + 0] = packed[i].slot.W0_qp_max;
                    loc_vals[6 * i + 1] = packed[i].slot.d_uncapped;
                    loc_vals[6 * i + 2] = packed[i].slot.d_eff;
                    loc_vals[6 * i + 3] = packed[i].slot.P_m;
                    loc_vals[6 * i + 4] = packed[i].slot.ruptured ? 1.0 : 0.0;
                    loc_vals[6 * i + 5] = packed[i].slot.fold_invalid ? 1.0 : 0.0;
                }

                std::vector<int> val_counts(static_cast<std::size_t>(nprocs), 0);
                std::vector<int> val_displs(static_cast<std::size_t>(nprocs), 0);
                for (int r = 0; r < nprocs; ++r)
                {
                    val_counts[r] = 6 * counts[r];
                    val_displs[r] = 6 * displs[r];
                }

                std::vector<unsigned long long> glob_keys(static_cast<std::size_t>(total_n));
                std::vector<double> glob_vals(static_cast<std::size_t>(6 * total_n));
                MPI_Allgatherv(local_n > 0 ? loc_keys.data() : nullptr,
                               local_n,
                               MPI_UNSIGNED_LONG_LONG,
                               total_n > 0 ? glob_keys.data() : nullptr,
                               counts.data(),
                               displs.data(),
                               MPI_UNSIGNED_LONG_LONG,
                               comm);
                MPI_Allgatherv(local_n > 0 ? loc_vals.data() : nullptr,
                               6 * local_n,
                               MPI_DOUBLE,
                               total_n > 0 ? glob_vals.data() : nullptr,
                               val_counts.data(),
                               val_displs.data(),
                               MPI_DOUBLE,
                               comm);

                packed.clear();
                packed.reserve(static_cast<std::size_t>(total_n));
                for (int i = 0; i < total_n; ++i)
                {
                    SurfaceDamageQPSlot slot;
                    slot.W0_qp_max = glob_vals[6 * i + 0];
                    slot.d_uncapped = glob_vals[6 * i + 1];
                    slot.d_eff = glob_vals[6 * i + 2];
                    slot.P_m = glob_vals[6 * i + 3];
                    slot.ruptured = (glob_vals[6 * i + 4] >= 0.5);
                    slot.fold_invalid = (glob_vals[6 * i + 5] >= 0.5);
                    packed.push_back({ static_cast<std::uint64_t>(glob_keys[i]), slot });
                }
            }
        }

        std::sort(packed.begin(), packed.end(), [](const PackedSlot& a, const PackedSlot& b) { return a.key < b.key; });

        std::vector<PackedSlot> merged;
        merged.reserve(packed.size());
        for (const auto& item : packed)
        {
            if (merged.empty() || merged.back().key != item.key)
            {
                merged.push_back(item);
            }
            else
            {
                auto& dst = merged.back().slot;
                dst.W0_qp_max = std::max(dst.W0_qp_max, item.slot.W0_qp_max);
                dst.d_uncapped = std::max(dst.d_uncapped, item.slot.d_uncapped);
                dst.d_eff = std::max(dst.d_eff, item.slot.d_eff);
                dst.P_m = std::max(dst.P_m, item.slot.P_m);
                dst.ruptured = dst.ruptured || item.slot.ruptured;
                dst.fold_invalid = dst.fold_invalid || item.slot.fold_invalid;
            }
        }

        db->putInteger("damage_tracker_version", 2);
        db->putDouble("damage_committed_time", committed_time_out);
        const int n = static_cast<int>(merged.size());
        db->putInteger("damage_n_slots", n);
        if (n == 0) return;

        std::vector<double> keys_hi(n), keys_lo(n), w0_max(n), d_uncapped(n), d_eff(n), p_m(n);
        std::vector<int> ruptured(n), fold_invalid(n);
        for (int idx = 0; idx < n; ++idx)
        {
            const std::uint64_t key = merged[idx].key;
            const auto& slot = merged[idx].slot;
            keys_hi[idx] = static_cast<double>(key >> 32);
            keys_lo[idx] = static_cast<double>(key & 0xFFFFFFFFull);
            w0_max[idx] = slot.W0_qp_max;
            d_uncapped[idx] = slot.d_uncapped;
            d_eff[idx] = slot.d_eff;
            p_m[idx] = slot.P_m;
            ruptured[idx] = slot.ruptured ? 1 : 0;
            fold_invalid[idx] = slot.fold_invalid ? 1 : 0;
        }
        db->putDoubleArray("damage_keys_hi", keys_hi.data(), n);
        db->putDoubleArray("damage_keys_lo", keys_lo.data(), n);
        db->putDoubleArray("damage_W0_max", w0_max.data(), n);
        db->putDoubleArray("damage_d_uncapped", d_uncapped.data(), n);
        db->putDoubleArray("damage_d_eff", d_eff.data(), n);
        db->putDoubleArray("damage_P_m", p_m.data(), n);
        db->putIntegerArray("damage_ruptured", ruptured.data(), n);
        db->putIntegerArray("damage_fold_invalid", fold_invalid.data(), n);
    }

    void getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
    {
        if (!db) return;
        d_slots.clear();
        d_committed_time = db->getDoubleWithDefault("damage_committed_time", 0.0);
        const int n = db->getIntegerWithDefault("damage_n_slots", 0);
        if (n <= 0) return;

        std::vector<double> keys_hi(n), keys_lo(n), w0_max(n), d_uncapped(n), d_eff(n), p_m(n);
        std::vector<int> ruptured(n), fold_invalid(n);
        db->getDoubleArray("damage_keys_hi", keys_hi.data(), n);
        db->getDoubleArray("damage_keys_lo", keys_lo.data(), n);
        db->getDoubleArray("damage_W0_max", w0_max.data(), n);
        db->getDoubleArray("damage_d_uncapped", d_uncapped.data(), n);
        db->getDoubleArray("damage_d_eff", d_eff.data(), n);
        db->getDoubleArray("damage_P_m", p_m.data(), n);
        db->getIntegerArray("damage_ruptured", ruptured.data(), n);
        db->getIntegerArray("damage_fold_invalid", fold_invalid.data(), n);

        for (int idx = 0; idx < n; ++idx)
        {
            const std::uint64_t hi = static_cast<std::uint64_t>(std::llround(keys_hi[idx]));
            const std::uint64_t lo = static_cast<std::uint64_t>(std::llround(keys_lo[idx]));
            const std::uint64_t key = (hi << 32) | (lo & 0xFFFFFFFFull);
            SurfaceDamageQPSlot slot;
            slot.W0_qp_max = w0_max[idx];
            slot.d_uncapped = d_uncapped[idx];
            slot.d_eff = d_eff[idx];
            slot.P_m = p_m[idx];
            slot.ruptured = (ruptured[idx] != 0);
            slot.fold_invalid = (fold_invalid[idx] != 0);
            d_slots[key] = slot;
        }
    }

private:
    SurfaceDamageParams d_params;
    double d_committed_time = 0.0;
    std::unordered_map<std::uint64_t, SurfaceDamageQPSlot> d_slots;
};

class SurfaceNonlocalHelmholtzFilter
{
public:
    SurfaceNonlocalHelmholtzFilter() = default;

    ~SurfaceNonlocalHelmholtzFilter()
    {
        destroy();
    }

    SurfaceNonlocalHelmholtzFilter(const SurfaceNonlocalHelmholtzFilter&) = delete;
    SurfaceNonlocalHelmholtzFilter& operator=(const SurfaceNonlocalHelmholtzFilter&) = delete;

    void destroy()
    {
        PetscBool petsc_init = PETSC_FALSE;
        PetscBool petsc_fini = PETSC_TRUE;
        PetscInitialized(&petsc_init);
        PetscFinalized(&petsc_fini);
        if (petsc_init && !petsc_fini)
        {
            if (d_scatter)
            {
                VecScatterDestroy(&d_scatter);
                d_scatter = nullptr;
            }
            if (d_ybar_loc)
            {
                VecDestroy(&d_ybar_loc);
                d_ybar_loc = nullptr;
            }
            if (d_ybar)
            {
                VecDestroy(&d_ybar);
                d_ybar = nullptr;
            }
            if (d_b)
            {
                VecDestroy(&d_b);
                d_b = nullptr;
            }
            if (d_ksp)
            {
                KSPDestroy(&d_ksp);
                d_ksp = nullptr;
            }
            if (d_A)
            {
                MatDestroy(&d_A);
                d_A = nullptr;
            }
        }
        else
        {
            d_scatter = nullptr;
            d_ybar_loc = nullptr;
            d_ybar = nullptr;
            d_b = nullptr;
            d_ksp = nullptr;
            d_A = nullptr;
        }
        d_local_elems.clear();
        d_local_surf_indices.clear();
        d_valid = false;
    }

    void ensureInitialized(const libMesh::MeshBase& mesh,
                           const libMesh::FEType& fe_type,
                           const libMesh::QuadratureType quad_type,
                           const libMesh::Order quad_order,
                           const double l_d)
    {
        TBOX_ASSERT(l_d > 0.0);
        const std::size_t n_nodes = static_cast<std::size_t>(mesh.n_nodes());
        const std::size_t n_active = static_cast<std::size_t>(mesh.n_active_elem());
        const std::size_t n_active_local = static_cast<std::size_t>(mesh.n_active_local_elem());
        if (d_valid && d_cached_n_nodes == n_nodes && d_cached_n_active_elem == n_active &&
            d_cached_n_active_local_elem == n_active_local && d_cached_quad_type == quad_type &&
            d_cached_order == quad_order && d_cached_fe_family == fe_type.family &&
            d_cached_fe_order == fe_type.order && std::abs(d_cached_ld - l_d) <= 1.0e-15 * std::max(1.0, l_d))
        {
            return;
        }

        destroy();
        d_cached_n_nodes = n_nodes;
        d_cached_n_active_elem = n_active;
        d_cached_n_active_local_elem = n_active_local;
        d_cached_ld = l_d;
        d_cached_quad_type = quad_type;
        d_cached_order = quad_order;
        d_cached_fe_family = fe_type.family;
        d_cached_fe_order = fe_type.order;

        const MPI_Comm comm = mesh.comm().get();
        int nprocs = 1;
        MPI_Comm_size(comm, &nprocs);

        std::vector<unsigned long long> loc_node_ids;
        for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
        {
            const libMesh::Elem* const elem = *el_it;
            for (unsigned int a = 0; a < elem->n_nodes(); ++a)
            {
                loc_node_ids.push_back(static_cast<unsigned long long>(elem->node_id(a)));
            }
        }
        std::sort(loc_node_ids.begin(), loc_node_ids.end());
        loc_node_ids.erase(std::unique(loc_node_ids.begin(), loc_node_ids.end()), loc_node_ids.end());

        std::vector<unsigned long long> all_node_ids;
        if (nprocs > 1)
        {
            const int loc_cnt = static_cast<int>(loc_node_ids.size());
            std::vector<int> counts(static_cast<std::size_t>(nprocs), 0);
            MPI_Allgather(&loc_cnt, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
            std::vector<int> displs(static_cast<std::size_t>(nprocs), 0);
            int total = 0;
            for (int r = 0; r < nprocs; ++r)
            {
                displs[r] = total;
                total += counts[r];
            }
            all_node_ids.resize(static_cast<std::size_t>(total));
            MPI_Allgatherv(loc_cnt > 0 ? loc_node_ids.data() : nullptr,
                           loc_cnt,
                           MPI_UNSIGNED_LONG_LONG,
                           total > 0 ? all_node_ids.data() : nullptr,
                           counts.data(),
                           displs.data(),
                           MPI_UNSIGNED_LONG_LONG,
                           comm);
            std::sort(all_node_ids.begin(), all_node_ids.end());
            all_node_ids.erase(std::unique(all_node_ids.begin(), all_node_ids.end()), all_node_ids.end());
        }
        else
        {
            all_node_ids = std::move(loc_node_ids);
        }

        const PetscInt N_surf = static_cast<PetscInt>(all_node_ids.size());
        if (N_surf == 0)
        {
            d_valid = true;
            return;
        }

        std::unordered_map<std::uint64_t, PetscInt> nid_to_idx;
        nid_to_idx.reserve(all_node_ids.size() * 2 + 1);
        for (PetscInt i = 0; i < N_surf; ++i)
        {
            nid_to_idx.emplace(static_cast<std::uint64_t>(all_node_ids[static_cast<std::size_t>(i)]), i);
        }

        PetscErrorCode ierr = 0;
        ierr = VecCreateMPI(comm, PETSC_DECIDE, N_surf, &d_b);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(d_b, &d_ybar);
        IBTK_CHKERRQ(ierr);

        PetscInt n_loc = 0;
        ierr = VecGetLocalSize(d_b, &n_loc);
        IBTK_CHKERRQ(ierr);

        ierr = MatCreate(comm, &d_A);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetSizes(d_A, n_loc, n_loc, N_surf, N_surf);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetType(d_A, MATAIJ);
        IBTK_CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(d_A, 48, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatMPIAIJSetPreallocation(d_A, 48, nullptr, 48, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetOption(d_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetOption(d_A, MAT_SYMMETRIC, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetOption(d_A, MAT_SPD, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);

        const unsigned int dim = mesh.mesh_dimension();
        std::unique_ptr<libMesh::QBase> qrule = libMesh::QBase::build(quad_type, dim, quad_order);
        std::unique_ptr<libMesh::FEBase> fe = libMesh::FEBase::build(dim, fe_type);
        fe->attach_quadrature_rule(qrule.get());
        const auto& phi = fe->get_phi();
        const auto& dphi_dxi = fe->get_dphidxi();
        const auto& dphi_deta = fe->get_dphideta();
        const auto& dX_dxi = fe->get_dxyzdxi();
        const auto& dX_deta = fe->get_dxyzdeta();
        const auto& qweights = qrule->get_weights();

        std::vector<PetscInt> ghost_indices;
        const double ld2 = l_d * l_d;

        for (auto el_it = mesh.active_local_elements_begin(); el_it != mesh.active_local_elements_end(); ++el_it)
        {
            libMesh::Elem* const elem = *el_it;
            fe->reinit(elem);
            const unsigned int n_basis = elem->n_nodes();
            const unsigned int n_qp = qrule->n_points();

            CachedLocalElem rec;
            rec.elem_id = static_cast<std::uint64_t>(elem->id());
            rec.global_dofs.resize(n_basis);
            for (unsigned int a = 0; a < n_basis; ++a)
            {
                const PetscInt gidx = nid_to_idx.at(static_cast<std::uint64_t>(elem->node_id(a)));
                rec.global_dofs[a] = gidx;
                ghost_indices.push_back(gidx);
            }
            rec.phi_qp_vals.resize(static_cast<std::size_t>(n_qp) * n_basis, 0.0);
            rec.jxw0.resize(n_qp, 0.0);
            rec.n_basis = n_basis;
            rec.n_qp = n_qp;

            std::vector<double> ke(static_cast<std::size_t>(n_basis) * n_basis, 0.0);
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                SurfaceTangents G0;
                for (unsigned int i = 0; i < 3; ++i)
                {
                    G0(i, 0) = dX_dxi[qp](i);
                    G0(i, 1) = dX_deta[qp](i);
                }
                const Metric2x2 M0 = computeSurfaceMetric(G0);
                TBOX_ASSERT(M0.det > 1.0e-30);
                const double w0 = std::sqrt(M0.det) * qweights[qp];
                rec.jxw0[qp] = w0;

                const double inv00 = M0.m11 / M0.det;
                const double inv01 = -M0.m01 / M0.det;
                const double inv11 = M0.m00 / M0.det;

                for (unsigned int a = 0; a < n_basis; ++a)
                {
                    rec.phi_qp_vals[static_cast<std::size_t>(qp) * n_basis + a] = phi[a][qp];
                    const double ga0 = dphi_dxi[a][qp];
                    const double ga1 = dphi_deta[a][qp];
                    const double minv_ga0 = inv00 * ga0 + inv01 * ga1;
                    const double minv_ga1 = inv01 * ga0 + inv11 * ga1;
                    for (unsigned int b = 0; b < n_basis; ++b)
                    {
                        const double grad_dot = minv_ga0 * dphi_dxi[b][qp] + minv_ga1 * dphi_deta[b][qp];
                        ke[static_cast<std::size_t>(a) * n_basis + b] +=
                            w0 * (phi[a][qp] * phi[b][qp] + ld2 * grad_dot);
                    }
                }
            }

            ierr = MatSetValues(d_A,
                                static_cast<PetscInt>(n_basis),
                                rec.global_dofs.data(),
                                static_cast<PetscInt>(n_basis),
                                rec.global_dofs.data(),
                                ke.data(),
                                ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            d_local_elems.push_back(std::move(rec));
        }

        ierr = MatAssemblyBegin(d_A, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(d_A, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);

        ierr = KSPCreate(comm, &d_ksp);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(d_ksp, d_A, d_A);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(d_ksp, KSPCG);
        IBTK_CHKERRQ(ierr);
        PC pc = nullptr;
        ierr = KSPGetPC(d_ksp, &pc);
        IBTK_CHKERRQ(ierr);
        ierr = PCSetType(pc, PCJACOBI);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetTolerances(d_ksp, 1.0e-12, 1.0e-18, PETSC_DEFAULT, 2000);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetUp(d_ksp);
        IBTK_CHKERRQ(ierr);

        std::sort(ghost_indices.begin(), ghost_indices.end());
        ghost_indices.erase(std::unique(ghost_indices.begin(), ghost_indices.end()), ghost_indices.end());
        d_local_surf_indices = std::move(ghost_indices);

        std::unordered_map<PetscInt, PetscInt> glob_to_loc;
        glob_to_loc.reserve(d_local_surf_indices.size() * 2 + 1);
        for (PetscInt li = 0; li < static_cast<PetscInt>(d_local_surf_indices.size()); ++li)
        {
            glob_to_loc.emplace(d_local_surf_indices[static_cast<std::size_t>(li)], li);
        }
        for (auto& rec : d_local_elems)
        {
            rec.local_dofs.resize(rec.n_basis);
            for (unsigned int a = 0; a < rec.n_basis; ++a)
            {
                rec.local_dofs[a] = glob_to_loc.at(rec.global_dofs[a]);
            }
        }

        const PetscInt n_loc_surf = static_cast<PetscInt>(d_local_surf_indices.size());
        ierr = VecCreateSeq(PETSC_COMM_SELF, n_loc_surf, &d_ybar_loc);
        IBTK_CHKERRQ(ierr);
        IS is_glob = nullptr;
        IS is_loc = nullptr;
        ierr = ISCreateGeneral(
            comm, n_loc_surf, n_loc_surf > 0 ? d_local_surf_indices.data() : nullptr, PETSC_COPY_VALUES, &is_glob);
        IBTK_CHKERRQ(ierr);
        ierr = ISCreateStride(PETSC_COMM_SELF, n_loc_surf, 0, 1, &is_loc);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterCreate(d_ybar, is_glob, d_ybar_loc, is_loc, &d_scatter);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&is_glob);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&is_loc);
        IBTK_CHKERRQ(ierr);

        d_valid = true;
    }

    void filterQuadratureField(const std::vector<std::vector<double>>& y_local_qp,
                               std::vector<std::vector<double>>& ybar_local_qp) const
    {
        TBOX_ASSERT(d_valid);
        TBOX_ASSERT(y_local_qp.size() == d_local_elems.size());
        ybar_local_qp.resize(d_local_elems.size());
        if (d_local_surf_indices.empty() && !d_b)
        {
            return;
        }

        PetscErrorCode ierr = 0;
        ierr = VecZeroEntries(d_b);
        IBTK_CHKERRQ(ierr);

        std::vector<double> be;
        for (std::size_t ei = 0; ei < d_local_elems.size(); ++ei)
        {
            const auto& rec = d_local_elems[ei];
            const auto& yq = y_local_qp[ei];
            TBOX_ASSERT(yq.size() == rec.n_qp);
            be.assign(rec.n_basis, 0.0);
            for (unsigned int qp = 0; qp < rec.n_qp; ++qp)
            {
                const double w0_yq = rec.jxw0[qp] * yq[qp];
                const double* phi_row = &rec.phi_qp_vals[static_cast<std::size_t>(qp) * rec.n_basis];
                for (unsigned int a = 0; a < rec.n_basis; ++a)
                {
                    be[a] += phi_row[a] * w0_yq;
                }
            }
            ierr = VecSetValues(d_b, static_cast<PetscInt>(rec.n_basis), rec.global_dofs.data(), be.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }

        ierr = VecAssemblyBegin(d_b);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(d_b);
        IBTK_CHKERRQ(ierr);

        ierr = KSPSolve(d_ksp, d_b, d_ybar);
        IBTK_CHKERRQ(ierr);
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        ierr = KSPGetConvergedReason(d_ksp, &reason);
        IBTK_CHKERRQ(ierr);
        if (reason <= 0)
        {
            TBOX_ERROR("SurfaceNonlocalHelmholtzFilter: KSPSolve failed to converge (reason=" << reason << ")\n");
        }

        ierr = VecScatterBegin(d_scatter, d_ybar, d_ybar_loc, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(d_scatter, d_ybar, d_ybar_loc, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);

        const PetscScalar* ybar_arr = nullptr;
        ierr = VecGetArrayRead(d_ybar_loc, &ybar_arr);
        IBTK_CHKERRQ(ierr);

        for (std::size_t ei = 0; ei < d_local_elems.size(); ++ei)
        {
            const auto& rec = d_local_elems[ei];
            auto& out_q = ybar_local_qp[ei];
            out_q.assign(rec.n_qp, 0.0);
            for (unsigned int qp = 0; qp < rec.n_qp; ++qp)
            {
                const double* phi_row = &rec.phi_qp_vals[static_cast<std::size_t>(qp) * rec.n_basis];
                double val = 0.0;
                for (unsigned int a = 0; a < rec.n_basis; ++a)
                {
                    val += phi_row[a] * static_cast<double>(ybar_arr[rec.local_dofs[a]]);
                }
                out_q[qp] = std::max(0.0, val);
            }
        }

        ierr = VecRestoreArrayRead(d_ybar_loc, &ybar_arr);
        IBTK_CHKERRQ(ierr);
    }

    int lastSolverIterations() const
    {
        if (!d_ksp) return 0;
        PetscInt its = 0;
        KSPGetIterationNumber(d_ksp, &its);
        return static_cast<int>(its);
    }

    std::size_t numLocalElems() const
    {
        return d_local_elems.size();
    }

private:
    struct CachedLocalElem
    {
        std::uint64_t elem_id = 0;
        unsigned int n_basis = 0;
        unsigned int n_qp = 0;
        std::vector<PetscInt> global_dofs;
        std::vector<PetscInt> local_dofs;
        std::vector<double> phi_qp_vals;
        std::vector<double> jxw0;
    };

    bool d_valid = false;
    std::size_t d_cached_n_nodes = 0;
    std::size_t d_cached_n_active_elem = 0;
    std::size_t d_cached_n_active_local_elem = 0;
    double d_cached_ld = 0.0;
    libMesh::QuadratureType d_cached_quad_type = libMesh::INVALID_Q_RULE;
    libMesh::Order d_cached_order = libMesh::INVALID_ORDER;
    libMesh::FEFamily d_cached_fe_family = libMesh::INVALID_FE;
    libMesh::Order d_cached_fe_order = libMesh::INVALID_ORDER;

    Mat d_A = nullptr;
    Vec d_b = nullptr;
    Vec d_ybar = nullptr;
    Vec d_ybar_loc = nullptr;
    KSP d_ksp = nullptr;
    VecScatter d_scatter = nullptr;

    std::vector<PetscInt> d_local_surf_indices;
    std::vector<CachedLocalElem> d_local_elems;
};

} // namespace FESurfaceDamage
} // namespace IBAMR

#endif // #ifdef IBAMR_HAVE_LIBMESH
#endif // #ifndef included_IBAMR_IBFESurfaceDamage
