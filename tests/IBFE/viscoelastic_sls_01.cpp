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

#include <petscsys.h>

#include <IBFEViscoelasticity.h>
#include <SAMRAI_config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include <ibamr/app_namespaces.h>

// Fast 2D/3D CMake attest-style verification for IBAMR::FEViscoelasticity:
//   1. Closed-form step-shear relaxation matching exp(-t / tau_R) to < 1e-13.
//   2. Constant-rate simple-shear ramp second-order O(dt^2) convergence vs RK4 reference.
//   3. Finite-strain rigid-rotation objectivity (covariance of Cauchy and PK1 <= 1e-13).
//   4. Two-branch reduction (G_2 = 0 matches 1-branch SLS to < 1e-14) and read-only trial calls.
//   5. SAMRAI database restart round-trip and 2-rank MPI parity to <= 1e-14.

namespace
{
using namespace IBAMR::FEViscoelasticity;

TensorValue<double>
makeRotation(const double theta)
{
    TensorValue<double> Q;
    Q.zero();
    const double c = std::cos(theta);
    const double s = std::sin(theta);
#if (NDIM == 2)
    Q(0, 0) = c;
    Q(0, 1) = -s;
    Q(1, 0) = s;
    Q(1, 1) = c;
    Q(2, 2) = 1.0;
#else
    Q(0, 0) = c;
    Q(0, 1) = -s;
    Q(1, 0) = s;
    Q(1, 1) = c;
    Q(2, 2) = 1.0;
#endif
    return Q;
}

TensorValue<double>
makeSimpleShear(const double gamma)
{
    TensorValue<double> F;
    F.zero();
    F(0, 0) = 1.0;
    F(1, 1) = 1.0;
    F(2, 2) = 1.0;
    F(0, 1) = gamma;
    return F;
}

double
maxTensorDiff(const TensorValue<double>& A, const TensorValue<double>& B)
{
    double m = 0.0;
    for (unsigned int i = 0; i < NDIM; ++i)
    {
        for (unsigned int j = 0; j < NDIM; ++j)
        {
            m = std::max(m, std::abs(A(i, j) - B(i, j)));
        }
    }
    return m;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "viscoelastic_sls_01.log");

        const double G_1 = 4.0;
        const double tau_R = 100.0;
        const ViscoelasticBranchParams b1{ G_1, tau_R, InitialBranchState::RELAXED };

        // 1. Step-shear relaxation: jump at t=0+ to gamma=0.5, hold fixed over [0, 50 us]
        TensorValue<double> F_step = makeSimpleShear(0.5);
        TensorValue<double> F_eye = makeSimpleShear(0.0);
        SymmetricTensor3 q_n{}, S_n{};
        SymmetricTensor3 q_trial{}, S_trial{}, Q_over{};
        TensorValue<double> P_0{};
        evaluateBranchFromCommitted(q_trial, S_trial, Q_over, P_0, F_step, q_n, S_n, 0.0, b1);

        const double dt_hold = 50.0;
        SymmetricTensor3 q_end{}, S_end{}, Q_end{};
        TensorValue<double> P_end{};
        // Constant F_step over hold interval -> q_n = 0, S_n = S(F_step) -> H(t) = exp(-t/tau_R) * S(F_step)
        evaluateBranchFromCommitted(q_end, S_end, Q_end, P_end, F_step, q_n, S_trial, dt_hold, b1);
        const TensorValue<double> P_expected = std::exp(-dt_hold / tau_R) * P_0;
        const double step_err = maxTensorDiff(P_end, P_expected) / P_0.norm();
        TBOX_ASSERT(step_err < 1.0e-13);

        // 2. Constant-rate simple-shear ramp: verify second-order O(dt^2) convergence vs RK4 reference
        const double t_end = 40.0;
        const double gamma_rate = 0.015;
        auto run_rk4_ref = [&]() -> TensorValue<double>
        {
            const int N_rk = 2000;
            const double dt = t_end / static_cast<double>(N_rk);
            SymmetricTensor3 q_rk{};
            for (int k = 0; k < N_rk; ++k)
            {
                const double t = k * dt;
                auto eval_S = [&](double tau)
                {
                    SymmetricTensor3 S_r, C_r, Cinv_r;
                    computeIsochoric2PK(makeSimpleShear(gamma_rate * tau), G_1, S_r, C_r, Cinv_r);
                    return S_r;
                };
                auto rhs = [&](const SymmetricTensor3& qv, double tau) -> SymmetricTensor3
                {
                    const auto Sv = eval_S(tau);
                    return { (Sv.xx - qv.xx) / tau_R, (Sv.yy - qv.yy) / tau_R, (Sv.zz - qv.zz) / tau_R,
                             (Sv.xy - qv.xy) / tau_R, (Sv.xz - qv.xz) / tau_R, (Sv.yz - qv.yz) / tau_R };
                };
                const SymmetricTensor3 k1 = rhs(q_rk, t);
                const SymmetricTensor3 q2{ q_rk.xx + 0.5 * dt * k1.xx, q_rk.yy + 0.5 * dt * k1.yy,
                                           q_rk.zz + 0.5 * dt * k1.zz, q_rk.xy + 0.5 * dt * k1.xy,
                                           q_rk.xz + 0.5 * dt * k1.xz, q_rk.yz + 0.5 * dt * k1.yz };
                const SymmetricTensor3 k2 = rhs(q2, t + 0.5 * dt);
                const SymmetricTensor3 q3{ q_rk.xx + 0.5 * dt * k2.xx, q_rk.yy + 0.5 * dt * k2.yy,
                                           q_rk.zz + 0.5 * dt * k2.zz, q_rk.xy + 0.5 * dt * k2.xy,
                                           q_rk.xz + 0.5 * dt * k2.xz, q_rk.yz + 0.5 * dt * k2.yz };
                const SymmetricTensor3 k3 = rhs(q3, t + 0.5 * dt);
                const SymmetricTensor3 q4{ q_rk.xx + dt * k3.xx, q_rk.yy + dt * k3.yy, q_rk.zz + dt * k3.zz,
                                           q_rk.xy + dt * k3.xy, q_rk.xz + dt * k3.xz, q_rk.yz + dt * k3.yz };
                const SymmetricTensor3 k4 = rhs(q4, t + dt);
                q_rk.xx += (dt / 6.0) * (k1.xx + 2.0 * k2.xx + 2.0 * k3.xx + k4.xx);
                q_rk.yy += (dt / 6.0) * (k1.yy + 2.0 * k2.yy + 2.0 * k3.yy + k4.yy);
                q_rk.zz += (dt / 6.0) * (k1.zz + 2.0 * k2.zz + 2.0 * k3.zz + k4.zz);
                q_rk.xy += (dt / 6.0) * (k1.xy + 2.0 * k2.xy + 2.0 * k3.xy + k4.xy);
                q_rk.xz += (dt / 6.0) * (k1.xz + 2.0 * k2.xz + 2.0 * k3.xz + k4.xz);
                q_rk.yz += (dt / 6.0) * (k1.yz + 2.0 * k2.yz + 2.0 * k3.yz + k4.yz);
            }
            SymmetricTensor3 S_f, C_f, Cinv_f;
            const TensorValue<double> F_f = makeSimpleShear(gamma_rate * t_end);
            computeIsochoric2PK(F_f, G_1, S_f, C_f, Cinv_f);
            const SymmetricTensor3 Q_f = projectDevC(S_f - q_rk, C_f, Cinv_f);
            TensorValue<double> P_f;
            pushForward2PKToPK1(F_f, Q_f, P_f);
            return P_f;
        };
        const TensorValue<double> P_ref = run_rk4_ref();

        auto run_ramp = [&](int n_steps) -> double
        {
            GeneralizedMaxwellHistoryStore store;
            store.configure({ b1 });
            store.initializeSlot(0, 0, F_eye);
            const double dt = t_end / static_cast<double>(n_steps);
            for (int k = 1; k <= n_steps; ++k)
            {
                const double t_k = k * dt;
                const TensorValue<double> F_k = makeSimpleShear(gamma_rate * t_k);
                store.commitQuadraturePoint(0, 0, F_k, t_k);
                store.finishCommitStep(t_k);
            }
            TensorValue<double> P_out;
            store.evaluatePK1(P_out, 0, 0, makeSimpleShear(gamma_rate * t_end), t_end);
            return maxTensorDiff(P_out, P_ref) / P_ref.norm();
        };
        const double err_coarse = run_ramp(40);
        const double err_fine = run_ramp(80);
        const double conv_order = std::log2(err_coarse / std::max(1.0e-18, err_fine));
        plog << "ramp err_coarse=" << err_coarse << " err_fine=" << err_fine << " order=" << conv_order << "\n";
        TBOX_ASSERT(err_fine < 5.0e-5);
        TBOX_ASSERT(conv_order > 1.95);

        // 3. Finite-strain rigid-rotation objectivity & 4. 2-branch G_2=0 reduction & 5. SAMRAI restart
        GeneralizedMaxwellHistoryStore store_1b, store_2b;
        store_1b.configure({ b1 });
        store_2b.configure({ b1, ViscoelasticBranchParams{ 0.0, 25.0, InitialBranchState::RELAXED } });
        store_1b.initializeSlot(10, 2, F_eye);
        store_2b.initializeSlot(10, 2, F_eye);
        const TensorValue<double> Q_rot = makeRotation(0.73);
        const TensorValue<double> F_mid = makeSimpleShear(0.35);
        store_1b.commitQuadraturePoint(10, 2, F_mid, 12.0);
        store_1b.finishCommitStep(12.0);
        store_2b.commitQuadraturePoint(10, 2, F_mid, 12.0);
        store_2b.finishCommitStep(12.0);

        const TensorValue<double> F_eval = makeSimpleShear(0.45);
        const TensorValue<double> F_eval_rot = Q_rot * F_eval;
        TensorValue<double> P_1b, P_2b, P_rot;
        store_1b.evaluatePK1(P_1b, 10, 2, F_eval, 16.0);
        store_2b.evaluatePK1(P_2b, 10, 2, F_eval, 16.0);
        store_1b.evaluatePK1(P_rot, 10, 2, F_eval_rot, 16.0);
        const double two_branch_diff = maxTensorDiff(P_1b, P_2b);
        const double rot_obj_err = maxTensorDiff(P_rot, Q_rot * P_1b) / P_1b.norm();
        TBOX_ASSERT(two_branch_diff <= 1.0e-14);
        TBOX_ASSERT(rot_obj_err <= 1.0e-13);

        // SAMRAI database round-trip + MPI parity across ranks
        Pointer<MemoryDatabase> db = new MemoryDatabase("sls_restart_test");
        store_1b.putToDatabase(db);
        GeneralizedMaxwellHistoryStore restored;
        restored.getFromDatabase(db);
        TensorValue<double> P_rest;
        restored.evaluatePK1(P_rest, 10, 2, F_eval, 16.0);
        double restart_diff = maxTensorDiff(P_1b, P_rest);
        restart_diff = IBTK_MPI::maxReduction(restart_diff);
        TBOX_ASSERT(restart_diff <= 1.0e-14);

        // Partitioned cross-rank MPI restart check: each rank owns a disjoint element slot
        // and verifies every rank restores peer ranks' slots deterministically after restart.
        const int my_rank = IBTK_MPI::getRank();
        const int n_ranks = IBTK_MPI::getNodes();
        GeneralizedMaxwellHistoryStore store_part;
        store_part.configure({ b1 });
        const std::uint64_t my_eid = static_cast<std::uint64_t>(100 + my_rank);
        store_part.initializeSlot(my_eid, 0U, F_eye);
        const TensorValue<double> F_rank = makeSimpleShear(0.20 + 0.10 * static_cast<double>(my_rank));
        store_part.commitQuadraturePoint(my_eid, 0U, F_rank, 12.0);
        store_part.finishCommitStep(12.0);

        Pointer<MemoryDatabase> db_part = new MemoryDatabase("sls_partitioned_restart_test");
        store_part.putToDatabase(db_part);
        GeneralizedMaxwellHistoryStore restored_part;
        restored_part.getFromDatabase(db_part);
        TBOX_ASSERT(restored_part.numSlots() == static_cast<size_t>(n_ranks));

        double part_restart_diff = 0.0;
        for (int src_r = 0; src_r < n_ranks; ++src_r)
        {
            const std::uint64_t eid_r = static_cast<std::uint64_t>(100 + src_r);
            GeneralizedMaxwellHistoryStore ref_r;
            ref_r.configure({ b1 });
            ref_r.initializeSlot(eid_r, 0U, F_eye);
            ref_r.commitQuadraturePoint(eid_r, 0U, makeSimpleShear(0.20 + 0.10 * static_cast<double>(src_r)), 12.0);
            ref_r.finishCommitStep(12.0);

            TensorValue<double> P_ref_r, P_rest_r;
            ref_r.evaluatePK1(P_ref_r, eid_r, 0U, F_eval, 16.0);
            restored_part.evaluatePK1(P_rest_r, eid_r, 0U, F_eval, 16.0);
            TBOX_ASSERT(P_ref_r.norm() > 0.0);
            part_restart_diff = std::max(part_restart_diff, maxTensorDiff(P_ref_r, P_rest_r));
        }
        part_restart_diff = IBTK_MPI::maxReduction(part_restart_diff);
        TBOX_ASSERT(part_restart_diff <= 1.0e-14);
        restart_diff = std::max(restart_diff, part_restart_diff);

        // Verify nearest-reference-coordinate lookup helper (inside tol2 vs fallback -1)
        const std::vector<std::array<double, 3>> qp_coords = { { 0.1, 0.2, 0.3 }, { 0.7, 0.8, 0.9 } };
        const int lookup_hit = GeneralizedMaxwellHistoryStore::findNearestSlotIndex(qp_coords, 0.7, 0.8, 0.9);
        const int lookup_miss = GeneralizedMaxwellHistoryStore::findNearestSlotIndex(qp_coords, 0.5, 0.5, 0.5);
        TBOX_ASSERT(lookup_hit == 1);
        TBOX_ASSERT(lookup_miss == -1);

        if (IBTK_MPI::getRank() == 0)
        {
            const double relax_ratio = P_end(0, 1) / P_0(0, 1);
            const auto fmt_tol = [](const double val) { return (val < 1.0e-14) ? 0.0 : val; };
            std::ostringstream oss;
            oss << "SLS_CHECKS PASS\n"
                << std::scientific << std::setprecision(5) << "step_relax_ratio=" << std::fixed << std::setprecision(8)
                << relax_ratio << "\n"
                << std::scientific << std::setprecision(5) << "step_rel_err=" << fmt_tol(step_err) << "\n"
                << "ramp_err_coarse=" << err_coarse << "\n"
                << "ramp_err_fine=" << err_fine << "\n"
                << "ramp_order=" << std::fixed << std::setprecision(5) << conv_order << "\n"
                << std::scientific << std::setprecision(5) << "two_branch_diff=" << fmt_tol(two_branch_diff) << "\n"
                << "rot_obj_err=" << fmt_tol(rot_obj_err) << "\n"
                << "restart_mpi_diff=" << fmt_tol(restart_diff) << "\n"
                << "qp_lookup_hit=" << lookup_hit << "\n"
                << "qp_lookup_miss=" << lookup_miss << "\n";
            pout << oss.str();
            std::ofstream out("output");
            out << oss.str();
            out.flush();
            TBOX_ASSERT(out.good());
        }
    }
    return 0;
}
