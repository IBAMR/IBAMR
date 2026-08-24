# Multiplicative CAV current-master merge readiness

Date: 2026-08-24

This handoff covers only the multiplicative CAV recovery through `cav-6g4`. It does not implement, validate, or certify CAV-RAS.

## Outcome

The multiplicative CAV implementation is locally merge-ready on the current-master restack:

- immutable independent-audit authority: `f0d32afbd27d031748f884a8f8ebe6647d97e7cf`, suite `m5-f0d32afb-20260824T001159Z`, all six lanes and synthesis `PASS`, no blocking finding;
- target master: `54055638ea2e8a080e3617baba707a4506b86f4d`;
- patch-equivalent rehearsal: `82e8fd6c362e3fe2730e12827aeab7b7e88afdfb`, tree `035e54fb4452abc77db1cd80b389e0b56732c62a`, 94/94 range-diff entries `=`;
- documentation-only merge-preparation tip: local branch `codex/implicit-ib-cav-master-merge-preparation`, based on the rehearsal and changing no file outside `plans/implicit-ib/`;
- fresh Debug 2D/3D build: passed;
- focused feature `attest`: 39/39 twice;
- broad `attest -R IB/`: 24/24;
- all 41 master-based review refs in Section 6 of `cav-implementation-plan.md`: exact object IDs verified 41/41 and verified as rehearsal ancestors.

The independent audit and the current-master implementation validation are separate evidence. No new independent-audit claim is made for the rebased SHA.

## Feature/test boundary check

Every CAV review unit that changes production code also changes its focused native test source or fixtures in the same unit. The inventory below is based on the exact diff between consecutive master-based review tips. The tests use the normal IBAMR `attest` workflow; external paired replay is additional evidence and never substitutes for native coverage.

| Unit | One review claim | Same-unit native test seam |
|---|---|---|
| `cav-1a` | Expose a narrow nonowning view of live Stokes operator state. | Existing shell-semantics executable; `live_operator_state` fixture. |
| `cav-1c` | Make matrix IB kernels follow the live Fortran kernel formulas/order. | Shared 2D/3D matrix/live kernel-consistency executable and fixtures; affected legacy expected outputs. |
| `cav-2` | Add deterministic pressure-seeded RELAXED/STRICT patch construction. | Existing CAV-reference and transfer-parity executables, including 3D construction coverage. |
| `cav-3a` | Add the BLAS/LAPACK LU local-solve backend. | Existing shell-semantics executable; LU, additive, lifecycle, and singular expected-failure fixtures. |
| `cav-3b` | Add retained robust BLAS/LAPACK solve modes and failure policy. | Existing shell-semantics executable; robust, SVD, QR expected-failure, and Cholesky-policy fixtures. |
| `cav-5` | Add the pressure-seeded multiplicative CAV smoother to the production level/FAC path. | Existing implicit-Stokes-IB component executable plus shared CAV shell-semantics executable; RELAXED/STRICT fixtures. |
| `cav-6a` | Add passive FAC-cycle stage observation without changing the FAC algorithm. | Existing implicit-Stokes-IB component executable and pressure-CAV fixture. |
| `cav-6b` | Keep the outer Stokes-IB Jacobian matrix-free and cover regrid/FAC production behavior. | Existing implicit integrator and component executables; CAV regrid and pressure-CAV fixtures. |
| `cav-6c` | Centralize correction composition across concrete local backends. | Existing shell-semantics executable; exact stagewise fixture exercises all attached backends. |
| `cav-6d` | Replace per-patch full-domain residual actions with affected-row updates. | Existing shell-semantics executable; exact stagewise result and sparsity-pattern invalidation expected failure. |
| `cav-6e` | Remove the redundant fixed-LU alias/API surface. | Existing LU, additive-LU, and singular-LU fixtures redirected through the retained generic backend. |
| `cav-6f2d` | Fill RT0 fine ghosts before zero-presmoothing restriction. | Existing implicit-Stokes-IB component executable and legal pressure-CAV hierarchy fixture. |
| `cav-6g1` | Borrow creator-owned supplied Galerkin matrices without copying or taking lifecycle ownership. | Existing shell and implicit-FAC executables; handle-identity/lifetime and Galerkin borrowing fixtures. |
| `cav-6g2` | Build the legacy velocity-seed pair map only for its STRICT consumer. | Existing 2D/3D CAV construction executable; absence/count/order expected outputs. |
| `cav-6g3` | Represent unrestricted pressure-CAV ownership partitions as absent. | Existing CAV shell and implicit-FAC executables; representation and backend-lifecycle coverage. |
| `cav-6g4` | Cache hierarchy-dependent FAC residual PETSc work vectors. | Existing implicit-FAC executable; rediscretized/Galerkin lifecycle counts, repeated residuals, and reinitialization. |

The following review units introduce no production functionality and therefore do not create a production-code-now/tests-later boundary:

| Unit | Review role |
|---|---|
| `cav-0` | Authoritative plan, committed ledger, and independent-audit instructions. |
| `cav-1b` | Raw export/comparator contract and mutation-sensitive native test support. |
| `cav-4` | Backend-neutral multiplicative composition semantics regression. |
| `cav-6f1` | Native low-stiffness physical-residual diagnostics and fixtures. |
| `cav-6f2a` | Test-only live construction/export schema and comparator coverage. |
| `cav-6f2b` | Test-only local/FAC/V-cycle/FGMRES trace export. |
| `cav-6f2c` | External N0/N1 paired common-arithmetic replay tooling. |
| `cav-6f2e` | External N2-N4 paired replay and backend cross-check tooling. |
| `cav-6f3` | Native 3D construction/solve lifecycle coverage using the existing dimensional source. |

## Exact validation commands

Run from the configured build directory, with the checkout's actual `attest` path and explicit environment tool paths:

```bash
/private/tmp/ibamr-cav-master-restack-rehearsal/attest -j8 \
  --mpi-executable /opt/homebrew/bin/mpiexec \
  --numdiff-executable /opt/homebrew/bin/numdiff \
  --verbose \
  -R '^(IBTK/sc_interp_matrix_live_kernel_consistency|IB/implicit_stokes_ib_solver_components_01|navier_stokes/stokes_petsc_level_solver_(cav_shell_semantics|shell_multiplicative_semantics))'
```

This scope passed 39/39 twice. The separate broad scope was:

```bash
/private/tmp/ibamr-cav-master-restack-rehearsal/attest -j8 \
  --mpi-executable /opt/homebrew/bin/mpiexec \
  --numdiff-executable /opt/homebrew/bin/numdiff \
  --verbose -R IB/
```

It passed 24/24. The initial managed-sandbox attempt is excluded because macOS denied Open MPI socket setup and PETSc hostname access before solver initialization; the identical commands passed with the required host access. Test duration is CI operational metadata only and is not performance evidence.

## Local submission order

1. Submit the foundation `*-master` review branches in the exact Section 6 order through foundation 6.4.
2. Submit `cav-0-recovery-plan-master`, then the CAV `*-master` branches in order from `cav-1a` through the `cav-6g4` production boundary.
3. Carry the cumulative audit/finalization and current plan records from `codex/implicit-ib-cav-master-merge-preparation` only as clearly labeled documentation/evidence commits; do not backdate them into the audited SHA.
4. Rebase or resolve review feedback forward on new candidate commits. Do not rewrite any cleaned-foundation, historical, audited, or frozen-result ref.

No branch has been pushed and no PR has been opened.

## Retained nonblocking items

- S3 `N-001`: Release expected-output portability remains nonblocking.
- The improvised Lane C `16x16 -> 32x32 -> 64x64`, `K=1e4` carrier did not converge. Do not advertise that configuration until it is separately diagnosed.
- CAV-RAS remains outside this merge scope and begins only after the user decides the multiplicative-CAV submission has reached the desired upstream state.
