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

// Objective finite-strain isochoric Standard Linear Solid (SLS) /
// generalized Maxwell internal-variable model (C++17, promotes cleanly
// from example helper to include/ibamr/):
//   S_{m,iso}(C) = G_m J^(-2/3) (I - (I_1 / 3) C^(-1)),  where I is the
//   3x3 identity tensor, paired with exact exponential linear-ramp
//   integrator phi_1(x), psi_1(x) and accepted-step SAMRAI restart state.
//
// Sources:
//   - Simo, J. C. (1987), Comput. Methods Appl. Mech. Eng. 60:153-173.
//   - Simo, J. C. & Hughes, T. J. R. (1998), Computational Inelasticity, Ch. 10.
//   - Holzapfel, G. A. (2000), Nonlinear Solid Mechanics, Ch. 6.
//
// Quadrature-slot lookup note:
//   Because IBAMR FEMechanicsBase PK1 callbacks pass reference coordinate s
//   without a quadrature-point index qp, live IBFEMethod callers match s to
//   stored Gauss coordinates via nearest squared distance <= elem_qp_tol2
//   (requiring fixed non-adaptive quadrature, with centroid fallback for
//   non-Gauss evaluations).

#ifndef included_IBAMR_IBFEViscoelasticity
#define included_IBAMR_IBFEViscoelasticity

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <libmesh/tensor_value.h>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace IBAMR
{
namespace FEViscoelasticity
{
enum class InitialBranchState
{
    RELAXED = 0,
    LOADED = 1
};

struct SymmetricTensor3
{
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;

    inline double operator()(unsigned int i, unsigned int j) const
    {
        if (i > j) std::swap(i, j);
        if (i == 0 && j == 0) return xx;
        if (i == 1 && j == 1) return yy;
        if (i == 2 && j == 2) return zz;
        if (i == 0 && j == 1) return xy;
        if (i == 0 && j == 2) return xz;
        return yz;
    }

    inline SymmetricTensor3 operator+(const SymmetricTensor3& o) const
    {
        return { xx + o.xx, yy + o.yy, zz + o.zz, xy + o.xy, xz + o.xz, yz + o.yz };
    }

    inline SymmetricTensor3 operator-(const SymmetricTensor3& o) const
    {
        return { xx - o.xx, yy - o.yy, zz - o.zz, xy - o.xy, xz - o.xz, yz - o.yz };
    }
};

struct ViscoelasticBranchParams
{
    double G = 0.0;
    double tau_R = 1.0;
    InitialBranchState initial_state = InitialBranchState::RELAXED;
};

// Exponential linear-ramp integration weights phi_1(x) = (1 - exp(-x)) / x and psi_1(x) = 1 - phi_1(x).
// Uses a 5-term Taylor series for |x| < 1e-4 to avoid subtraction cancellation (Simo & Hughes 1998, Ch. 10).
inline double
phi1(const double x)
{
    if (std::abs(x) < 1.0e-4)
    {
        const double x2 = x * x;
        return 1.0 - 0.5 * x + (x2 / 6.0) - (x2 * x / 24.0) + (x2 * x2 / 120.0);
    }
    return -std::expm1(-x) / x;
}

inline double
psi1(const double x)
{
    if (std::abs(x) < 1.0e-4)
    {
        const double x2 = x * x;
        return 0.5 * x - (x2 / 6.0) + (x2 * x / 24.0) - (x2 * x2 / 120.0);
    }
    return 1.0 - phi1(x);
}

inline void
embed3x3(const libMesh::TensorValue<double>& FF, double (&F3)[3][3])
{
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            F3[i][j] = (i < NDIM && j < NDIM) ? FF(i, j) : 0.0;
        }
    }
#if (NDIM == 2)
    F3[2][2] = 1.0;
#endif
}

inline double
detDeformationGradient(const libMesh::TensorValue<double>& FF)
{
#if (NDIM == 2)
    return FF(0, 0) * FF(1, 1) - FF(0, 1) * FF(1, 0);
#else
    return FF.det();
#endif
}

// Compute driving symmetric isochoric 2nd Piola-Kirchhoff branch target:
//   S_{m,iso}(C) = G * J^(-2/3) * (I - (I_1 / 3) * C^-1).
inline void
computeIsochoric2PK(const libMesh::TensorValue<double>& FF,
                    const double G,
                    SymmetricTensor3& S_iso,
                    SymmetricTensor3& C_out,
                    SymmetricTensor3& C_inv_out)
{
    const double J = detDeformationGradient(FF);
    TBOX_ASSERT(std::isfinite(J) && J > 0.0);
    double F3[3][3];
    embed3x3(FF, F3);

    C_out.xx = F3[0][0] * F3[0][0] + F3[1][0] * F3[1][0] + F3[2][0] * F3[2][0];
    C_out.yy = F3[0][1] * F3[0][1] + F3[1][1] * F3[1][1] + F3[2][1] * F3[2][1];
    C_out.zz = F3[0][2] * F3[0][2] + F3[1][2] * F3[1][2] + F3[2][2] * F3[2][2];
    C_out.xy = F3[0][0] * F3[0][1] + F3[1][0] * F3[1][1] + F3[2][0] * F3[2][1];
    C_out.xz = F3[0][0] * F3[0][2] + F3[1][0] * F3[1][2] + F3[2][0] * F3[2][2];
    C_out.yz = F3[0][1] * F3[0][2] + F3[1][1] * F3[1][2] + F3[2][1] * F3[2][2];

    const double inv_J2 = 1.0 / (J * J);
    C_inv_out.xx = (C_out.yy * C_out.zz - C_out.yz * C_out.yz) * inv_J2;
    C_inv_out.yy = (C_out.xx * C_out.zz - C_out.xz * C_out.xz) * inv_J2;
    C_inv_out.zz = (C_out.xx * C_out.yy - C_out.xy * C_out.xy) * inv_J2;
    C_inv_out.xy = (C_out.xz * C_out.yz - C_out.xy * C_out.zz) * inv_J2;
    C_inv_out.xz = (C_out.xy * C_out.yz - C_out.yy * C_out.xz) * inv_J2;
    C_inv_out.yz = (C_out.xy * C_out.xz - C_out.xx * C_out.yz) * inv_J2;

    const double mean_I1 = (C_out.xx + C_out.yy + C_out.zz) / 3.0;
    const double coeff = G * std::pow(J, -2.0 / 3.0);
    S_iso.xx = coeff * (1.0 - mean_I1 * C_inv_out.xx);
    S_iso.yy = coeff * (1.0 - mean_I1 * C_inv_out.yy);
    S_iso.zz = coeff * (1.0 - mean_I1 * C_inv_out.zz);
    S_iso.xy = coeff * (-mean_I1 * C_inv_out.xy);
    S_iso.xz = coeff * (-mean_I1 * C_inv_out.xz);
    S_iso.yz = coeff * (-mean_I1 * C_inv_out.yz);
}

inline SymmetricTensor3
projectDevC(const SymmetricTensor3& H, const SymmetricTensor3& C, const SymmetricTensor3& C_inv)
{
    const double tr_HC = H.xx * C.xx + H.yy * C.yy + H.zz * C.zz + 2.0 * (H.xy * C.xy + H.xz * C.xz + H.yz * C.yz);
    const double c3 = tr_HC / 3.0;
    return { H.xx - c3 * C_inv.xx, H.yy - c3 * C_inv.yy, H.zz - c3 * C_inv.zz,
             H.xy - c3 * C_inv.xy, H.xz - c3 * C_inv.xz, H.yz - c3 * C_inv.yz };
}

inline void
pushForward2PKToPK1(const libMesh::TensorValue<double>& FF, const SymmetricTensor3& Q, libMesh::TensorValue<double>& PP)
{
    PP.zero();
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int A = 0; A < NDIM; ++A)
        {
            double s = 0.0;
            for (unsigned int B = 0; B < NDIM; ++B)
            {
                s += FF(i, B) * Q(B, A);
            }
            PP(i, A) = s;
        }
    }
}

// Read-only trial update of one Prony/Maxwell branch from committed (q_n, S_n) over dt_sub >= 0.
inline void
evaluateBranchFromCommitted(SymmetricTensor3& q_trial,
                            SymmetricTensor3& S_trial,
                            SymmetricTensor3& Q_overstress,
                            libMesh::TensorValue<double>& P_branch,
                            const libMesh::TensorValue<double>& FF,
                            const SymmetricTensor3& q_committed,
                            const SymmetricTensor3& S_committed,
                            const double dt_sub,
                            const ViscoelasticBranchParams& params)
{
    SymmetricTensor3 C_cur, C_inv;
    computeIsochoric2PK(FF, params.G, S_trial, C_cur, C_inv);
    if (params.G <= 0.0)
    {
        q_trial = SymmetricTensor3{};
        Q_overstress = SymmetricTensor3{};
        P_branch.zero();
        return;
    }

    const double dt_pos = std::max(0.0, dt_sub);
    SymmetricTensor3 H;
    if (dt_pos == 0.0)
    {
        q_trial = q_committed;
        H = S_trial - q_committed;
    }
    else
    {
        TBOX_ASSERT(params.tau_R > 0.0);
        const double x = dt_pos / params.tau_R;
        const double decay = std::exp(-x);
        const double one_minus_decay = -std::expm1(-x);
        const double phi = phi1(x);
        const double psi = psi1(x);
        const SymmetricTensor3 dSq = S_committed - q_committed;
        const SymmetricTensor3 dS = S_trial - S_committed;
        q_trial = { q_committed.xx + one_minus_decay * dSq.xx + psi * dS.xx,
                    q_committed.yy + one_minus_decay * dSq.yy + psi * dS.yy,
                    q_committed.zz + one_minus_decay * dSq.zz + psi * dS.zz,
                    q_committed.xy + one_minus_decay * dSq.xy + psi * dS.xy,
                    q_committed.xz + one_minus_decay * dSq.xz + psi * dS.xz,
                    q_committed.yz + one_minus_decay * dSq.yz + psi * dS.yz };
        H = { decay * dSq.xx + phi * dS.xx, decay * dSq.yy + phi * dS.yy, decay * dSq.zz + phi * dS.zz,
              decay * dSq.xy + phi * dS.xy, decay * dSq.xz + phi * dS.xz, decay * dSq.yz + phi * dS.yz };
    }
    Q_overstress = projectDevC(H, C_cur, C_inv);
    pushForward2PKToPK1(FF, Q_overstress, P_branch);
}

// Accepted-step multi-branch quadrature history container with SAMRAI restart serialization.
class GeneralizedMaxwellHistoryStore
{
public:
    static constexpr int s_restart_version = 1;

    GeneralizedMaxwellHistoryStore() = default;

    void configure(const std::vector<ViscoelasticBranchParams>& branches)
    {
        for (const auto& branch : branches)
        {
            TBOX_ASSERT(branch.tau_R > 0.0);
        }
        d_branches = branches;
    }

    void clear()
    {
        d_slots.clear();
        d_committed_time = 0.0;
    }

    size_t numSlots() const
    {
        return d_slots.size();
    }

    double committedTime() const
    {
        return d_committed_time;
    }

    // Nearest-reference-coordinate lookup helper for FEMechanicsBase PK1 callbacks
    // that receive s without a quadrature-point index (fixed non-adaptive rule).
    static int findNearestSlotIndex(const std::vector<std::array<double, 3>>& qp_ref_coords,
                                    const double s0,
                                    const double s1,
                                    const double s2,
                                    const double tol2 = 1.0e-12)
    {
        int best = -1;
        double best_d2 = tol2;
        for (int idx = 0; idx < static_cast<int>(qp_ref_coords.size()); ++idx)
        {
            const double dx = s0 - qp_ref_coords[idx][0];
            const double dy = s1 - qp_ref_coords[idx][1];
            const double dz = s2 - qp_ref_coords[idx][2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 <= best_d2)
            {
                best_d2 = d2;
                best = idx;
            }
        }
        return best;
    }

    static std::uint64_t makeKey(std::uint64_t elem_id, unsigned int qp)
    {
        return (elem_id << 8) | static_cast<std::uint64_t>(qp & 0xFFu);
    }

    void initializeSlot(std::uint64_t elem_id, unsigned int qp, const libMesh::TensorValue<double>& FF_0)
    {
        auto& vec = d_slots[makeKey(elem_id, qp)];
        vec.resize(d_branches.size());
        SymmetricTensor3 C, C_inv;
        for (size_t b = 0; b < d_branches.size(); ++b)
        {
            computeIsochoric2PK(FF_0, d_branches[b].G, vec[b].second, C, C_inv);
            if (d_branches[b].initial_state == InitialBranchState::LOADED)
            {
                vec[b].first = SymmetricTensor3{};
            }
            else
            {
                vec[b].first = vec[b].second;
            }
        }
    }

    // Read-only substep PK1 evaluation using committed step history.
    void evaluatePK1(libMesh::TensorValue<double>& PP_total,
                     std::uint64_t elem_id,
                     unsigned int qp,
                     const libMesh::TensorValue<double>& FF,
                     const double eval_time) const
    {
        PP_total.zero();
        const auto it = d_slots.find(makeKey(elem_id, qp));
        if (it == d_slots.end()) return;
        const double dt_sub = std::max(0.0, eval_time - d_committed_time);
        const auto& vec = it->second;
        for (size_t b = 0; b < d_branches.size(); ++b)
        {
            if (d_branches[b].G <= 0.0) continue;
            SymmetricTensor3 q_t, S_t, Q_t;
            libMesh::TensorValue<double> P_b;
            evaluateBranchFromCommitted(q_t, S_t, Q_t, P_b, FF, vec[b].first, vec[b].second, dt_sub, d_branches[b]);
            PP_total += P_b;
        }
    }

    // Advance committed state strictly once per accepted time step (rejects backward timestamps).
    void commitQuadraturePoint(std::uint64_t elem_id,
                               unsigned int qp,
                               const libMesh::TensorValue<double>& FF_new,
                               const double new_time)
    {
        TBOX_ASSERT(new_time >= d_committed_time - 1.0e-12);
        auto it = d_slots.find(makeKey(elem_id, qp));
        if (it == d_slots.end())
        {
            initializeSlot(elem_id, qp, FF_new);
            return;
        }
        const double dt_step = std::max(0.0, new_time - d_committed_time);
        auto& vec = it->second;
        for (size_t b = 0; b < d_branches.size(); ++b)
        {
            SymmetricTensor3 q_new, S_new, Q_unused;
            libMesh::TensorValue<double> P_unused;
            evaluateBranchFromCommitted(
                q_new, S_new, Q_unused, P_unused, FF_new, vec[b].first, vec[b].second, dt_step, d_branches[b]);
            vec[b].first = q_new;
            vec[b].second = S_new;
        }
    }

    void finishCommitStep(const double new_time)
    {
        TBOX_ASSERT(new_time >= d_committed_time - 1.0e-12);
        d_committed_time = new_time;
    }

    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) const
    {
        TBOX_ASSERT(db);
        int mpi_init = 0;
        MPI_Initialized(&mpi_init);
        const int nprocs = mpi_init ? IBTK::IBTK_MPI::getNodes() : 1;
        double committed_time_out = d_committed_time;
        if (mpi_init && nprocs > 1)
        {
            const MPI_Comm comm = IBTK::IBTK_MPI::getCommunicator();
            MPI_Allreduce(&d_committed_time, &committed_time_out, 1, MPI_DOUBLE, MPI_MAX, comm);
        }

        db->putInteger("sls_restart_version", s_restart_version);
        db->putDouble("sls_committed_time", committed_time_out);
        const int n_branches = static_cast<int>(d_branches.size());
        db->putInteger("sls_n_branches", n_branches);
        if (n_branches > 0)
        {
            std::vector<double> G_vec(n_branches), tau_vec(n_branches);
            std::vector<int> init_vec(n_branches);
            for (int b = 0; b < n_branches; ++b)
            {
                G_vec[b] = d_branches[b].G;
                tau_vec[b] = d_branches[b].tau_R;
                init_vec[b] = static_cast<int>(d_branches[b].initial_state);
            }
            db->putDoubleArray("sls_branch_G", G_vec.data(), n_branches);
            db->putDoubleArray("sls_branch_tau", tau_vec.data(), n_branches);
            db->putIntegerArray("sls_branch_init", init_vec.data(), n_branches);
        }

        if (n_branches == 0)
        {
            db->putInteger("sls_n_slots", 0);
            return;
        }

        const size_t vals_per_slot = static_cast<size_t>(n_branches) * 12u;
        struct PackedSlot
        {
            std::uint64_t key;
            std::vector<double> vals;
        };
        std::vector<PackedSlot> packed;
        packed.reserve(d_slots.size());
        for (const auto& kv : d_slots)
        {
            PackedSlot item;
            item.key = kv.first;
            item.vals.resize(vals_per_slot, 0.0);
            for (int b = 0; b < n_branches; ++b)
            {
                const auto& q = kv.second[b].first;
                const auto& S = kv.second[b].second;
                const size_t base = static_cast<size_t>(b) * 12u;
                item.vals[base + 0] = q.xx;
                item.vals[base + 1] = q.yy;
                item.vals[base + 2] = q.zz;
                item.vals[base + 3] = q.xy;
                item.vals[base + 4] = q.xz;
                item.vals[base + 5] = q.yz;
                item.vals[base + 6] = S.xx;
                item.vals[base + 7] = S.yy;
                item.vals[base + 8] = S.zz;
                item.vals[base + 9] = S.xy;
                item.vals[base + 10] = S.xz;
                item.vals[base + 11] = S.yz;
            }
            packed.push_back(std::move(item));
        }
        std::sort(packed.begin(), packed.end(), [](const PackedSlot& a, const PackedSlot& b) { return a.key < b.key; });

        if (mpi_init && nprocs > 1)
        {
            const MPI_Comm comm = IBTK::IBTK_MPI::getCommunicator();
            const int local_n = static_cast<int>(packed.size());
            std::vector<int> counts(static_cast<size_t>(nprocs), 0);
            MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

            std::vector<int> displs(static_cast<size_t>(nprocs), 0);
            int total_n = 0;
            for (int r = 0; r < nprocs; ++r)
            {
                displs[r] = total_n;
                total_n += counts[r];
            }

            if (total_n > 0)
            {
                const int stride = static_cast<int>(vals_per_slot);
                std::vector<unsigned long long> loc_keys(static_cast<size_t>(local_n));
                std::vector<double> loc_vals(static_cast<size_t>(local_n) * vals_per_slot);
                for (int i = 0; i < local_n; ++i)
                {
                    loc_keys[i] = static_cast<unsigned long long>(packed[i].key);
                    std::copy(packed[i].vals.begin(),
                              packed[i].vals.end(),
                              loc_vals.begin() + static_cast<size_t>(i) * vals_per_slot);
                }

                std::vector<int> val_counts(static_cast<size_t>(nprocs), 0);
                std::vector<int> val_displs(static_cast<size_t>(nprocs), 0);
                for (int r = 0; r < nprocs; ++r)
                {
                    val_counts[r] = stride * counts[r];
                    val_displs[r] = stride * displs[r];
                }

                std::vector<unsigned long long> glob_keys(static_cast<size_t>(total_n));
                std::vector<double> glob_vals(static_cast<size_t>(total_n) * vals_per_slot);
                MPI_Allgatherv(local_n > 0 ? loc_keys.data() : nullptr,
                               local_n,
                               MPI_UNSIGNED_LONG_LONG,
                               glob_keys.data(),
                               counts.data(),
                               displs.data(),
                               MPI_UNSIGNED_LONG_LONG,
                               comm);
                MPI_Allgatherv(local_n > 0 ? loc_vals.data() : nullptr,
                               stride * local_n,
                               MPI_DOUBLE,
                               glob_vals.data(),
                               val_counts.data(),
                               val_displs.data(),
                               MPI_DOUBLE,
                               comm);

                packed.clear();
                packed.reserve(static_cast<size_t>(total_n));
                for (int i = 0; i < total_n; ++i)
                {
                    PackedSlot item;
                    item.key = static_cast<std::uint64_t>(glob_keys[i]);
                    const auto src_beg = glob_vals.begin() + static_cast<size_t>(i) * vals_per_slot;
                    item.vals.assign(src_beg, src_beg + static_cast<std::ptrdiff_t>(vals_per_slot));
                    packed.push_back(std::move(item));
                }
                std::sort(packed.begin(),
                          packed.end(),
                          [](const PackedSlot& a, const PackedSlot& b) { return a.key < b.key; });
            }
            else
            {
                packed.clear();
            }
        }

        std::vector<PackedSlot> merged;
        merged.reserve(packed.size());
        for (auto& item : packed)
        {
            if (merged.empty() || merged.back().key != item.key)
            {
                merged.push_back(std::move(item));
            }
        }

        const int n_slots = static_cast<int>(merged.size());
        db->putInteger("sls_n_slots", n_slots);
        if (n_slots == 0) return;

        std::vector<double> keys_hi(n_slots), keys_lo(n_slots);
        std::vector<double> state_vals(static_cast<size_t>(n_slots) * vals_per_slot);
        for (int idx = 0; idx < n_slots; ++idx)
        {
            keys_hi[idx] = static_cast<double>(merged[idx].key >> 32);
            keys_lo[idx] = static_cast<double>(merged[idx].key & 0xFFFFFFFFull);
            std::copy(merged[idx].vals.begin(),
                      merged[idx].vals.end(),
                      state_vals.begin() + static_cast<size_t>(idx) * vals_per_slot);
        }
        db->putDoubleArray("sls_keys_hi", keys_hi.data(), n_slots);
        db->putDoubleArray("sls_keys_lo", keys_lo.data(), n_slots);
        db->putDoubleArray("sls_state_vals", state_vals.data(), static_cast<int>(state_vals.size()));
    }

    void getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
    {
        TBOX_ASSERT(db);
        const int version = db->getInteger("sls_restart_version");
        TBOX_ASSERT(version == s_restart_version);
        d_committed_time = db->getDouble("sls_committed_time");
        const int n_branches = db->getInteger("sls_n_branches");
        d_branches.clear();
        if (n_branches > 0)
        {
            std::vector<double> G_vec(n_branches), tau_vec(n_branches);
            std::vector<int> init_vec(n_branches);
            db->getDoubleArray("sls_branch_G", G_vec.data(), n_branches);
            db->getDoubleArray("sls_branch_tau", tau_vec.data(), n_branches);
            db->getIntegerArray("sls_branch_init", init_vec.data(), n_branches);
            d_branches.resize(n_branches);
            for (int b = 0; b < n_branches; ++b)
            {
                d_branches[b].G = G_vec[b];
                d_branches[b].tau_R = tau_vec[b];
                TBOX_ASSERT(d_branches[b].tau_R > 0.0);
                d_branches[b].initial_state = static_cast<InitialBranchState>(init_vec[b]);
            }
        }
        d_slots.clear();
        const int n_slots = db->getInteger("sls_n_slots");
        if (n_slots <= 0 || n_branches <= 0) return;
        std::vector<double> keys_hi(n_slots), keys_lo(n_slots);
        std::vector<double> state_vals(static_cast<size_t>(n_slots) * static_cast<size_t>(n_branches) * 12u);
        db->getDoubleArray("sls_keys_hi", keys_hi.data(), n_slots);
        db->getDoubleArray("sls_keys_lo", keys_lo.data(), n_slots);
        db->getDoubleArray("sls_state_vals", state_vals.data(), static_cast<int>(state_vals.size()));
        for (int idx = 0; idx < n_slots; ++idx)
        {
            const std::uint64_t key = (static_cast<std::uint64_t>(std::llround(keys_hi[idx])) << 32) |
                                      (static_cast<std::uint64_t>(std::llround(keys_lo[idx])) & 0xFFFFFFFFull);
            auto& vec = d_slots[key];
            vec.resize(n_branches);
            size_t base = static_cast<size_t>(idx) * static_cast<size_t>(n_branches) * 12u;
            for (int b = 0; b < n_branches; ++b)
            {
                vec[b].first = { state_vals[base + 0], state_vals[base + 1], state_vals[base + 2],
                                 state_vals[base + 3], state_vals[base + 4], state_vals[base + 5] };
                vec[b].second = { state_vals[base + 6], state_vals[base + 7],  state_vals[base + 8],
                                  state_vals[base + 9], state_vals[base + 10], state_vals[base + 11] };
                base += 12u;
            }
        }
    }

private:
    std::vector<ViscoelasticBranchParams> d_branches;
    double d_committed_time = 0.0;
    std::unordered_map<std::uint64_t, std::vector<std::pair<SymmetricTensor3, SymmetricTensor3>>> d_slots;
};

} // namespace FEViscoelasticity
} // namespace IBAMR

#endif // IBAMR_HAVE_LIBMESH
#endif // included_IBAMR_IBFEViscoelasticity
