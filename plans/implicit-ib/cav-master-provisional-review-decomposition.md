# Multiplicative CAV provisional review decomposition

Date: 2026-08-24

This is a local decomposition inventory, not a PR submission manifest. The boundaries below reproduce the current history exactly, but they require further consolidation/splitting analysis before they become a reviewable collection of branches. Do not push the 41 provisional `*-master` refs or open PRs from them. The separately published `origin/codex/implicit-ib-cav-multiplicative-current-master` branch is cumulative preservation/collaboration history, not an endorsed PR topology. Exact boundary object IDs and immutable provenance are in Section 6 of `cav-implementation-plan.md`; cumulative implementation validation is in `cav-master-integration-status.md`. This document does not cover CAV-RAS.

## Why this is not review-ready

- The literal reconstruction yields 15 prerequisite foundation boundaries and 25 multiplicative-CAV boundaries before documentation finalization. That is provenance, not a practical upstream stack.
- Several `cav-6f*` units exist to support independent reproduction and audit observability. We must decide which parts belong in upstream native tests, which belong in review-only evidence, and which should not be proposed for IBAMR master.
- The final branches must preserve same-unit feature tests while minimizing dependency depth and avoiding production refactors bundled with unrelated audit tooling.
- The cleaned foundation boundaries predate the prospective CAV feature-test rule and include the ledger-approved support-API test deferral. The final upstream decomposition must address that scope deliberately.
- No current `*-master` boundary should be treated as an approved PR head merely because its patch replay and cumulative integration tests are green.

## Evidence convention

Each PR body should distinguish:

- **same-unit native evidence:** the focused `attest` cases introduced or materially updated by that review unit;
- **cumulative current-master evidence:** fresh Debug 2D/3D build, focused multiplicative-CAV scope 39/39 twice, and separately labeled broad `attest -R IB/` 24/24 at rehearsal `82e8fd6c...`;
- **independent authority:** frozen M5 suite at exact `f0d32afb...`, all six lanes and synthesis `PASS`.

Do not report cumulative or independent evidence as though it were rerun at an intermediate PR tip. Before submission, use the focused selector listed below from that exact head and record the result in the PR body. Use the broad `-R IB/` scope required by the plan; a focused case does not substitute for it.

## Provisional foundation boundaries

| Order | Base branch | Head branch | Proposed title | One claim | Focused native evidence |
|---:|---|---|---|---|---|
| 1 | `master` at `54055638...` | `codex/implicit-ib-0-infrastructure-master` | Add PETSc options helper infrastructure | Add the common options helper required by the implicit-IB solver stack. | Existing options-helper consumers and dimensional library builds; this historical infrastructure boundary adds no standalone feature executable. |
| 2 | `codex/implicit-ib-0-infrastructure-master` | `codex/implicit-ib-1.1-support-apis-master` | Add implicit IB support APIs | Add the support surfaces used by later Stokes-IB operator components. | Historical ledger Finding 1 intentionally defers direct support-API coverage; do not claim a new focused case or reopen that disposition unless this API delta is changed. |
| 3 | `codex/implicit-ib-1.1-support-apis-master` | `codex/implicit-ib-1.2-stokes-ib-operator-master` | Add the Stokes-IB nonlinear operator | Add the nonlinear operator and shared time/operator plumbing. | `IB/implicit_stokes_ib_solver_components_01.max_levels=1`. |
| 4 | `codex/implicit-ib-1.2-stokes-ib-operator-master` | `codex/implicit-ib-1.3-stokes-ib-jacobian-operator-master` | Add the Stokes-IB Jacobian operator | Add the live Jacobian using the shared schedule and plumbing. | Existing one-level `implicit_stokes_ib_solver_components_01` Jacobian coverage. |
| 5 | `codex/implicit-ib-1.3-stokes-ib-jacobian-operator-master` | `codex/implicit-ib-1.4-stokes-ib-fac-level-operator-master` | Add Stokes-IB FAC level-operator support | Add the FAC level operator while retaining centralized residual-vector handling. | Existing focused solver-components test. |
| 6 | `codex/implicit-ib-1.4-stokes-ib-fac-level-operator-master` | `codex/implicit-ib-1.5-stokes-ib-jacobian-fac-preconditioner-master` | Add the Stokes-IB Jacobian FAC preconditioner | Add multilevel FAC orchestration for the live Jacobian. | `implicit_stokes_ib_solver_components_01` at one, two, and three levels. |
| 7 | `codex/implicit-ib-1.5-stokes-ib-jacobian-fac-preconditioner-master` | `codex/implicit-ib-2-integrator-operator-path-master` | Use shared Stokes-IB solver components in the integrator | Route the implicit integrator through the shared operator/time-schedule path. | Existing implicit-IB integration and solver-component cases. |
| 8 | `codex/implicit-ib-2-integrator-operator-path-master` | `codex/implicit-ib-3.1-cell-linear-prolongation-master` | Add linear cell-pressure prolongation | Add the cell-centered pressure prolongation operator and parity cases. | Linear pressure-prolongation parity fixtures. |
| 9 | `codex/implicit-ib-3.1-cell-linear-prolongation-master` | `codex/implicit-ib-3.2-cell-linear-restriction-master` | Add linear cell-pressure restriction | Add the matching cell-centered pressure restriction operator and parity cases. | Linear pressure-restriction parity fixtures. |
| 10 | `codex/implicit-ib-3.2-cell-linear-restriction-master` | `codex/implicit-ib-4-transfer-parity-tests-master` | Add Stokes transfer-operator parity coverage | Add transfer parity and RT0 depth/preallocation regression coverage without new solver behavior. | Transfer parity executable/fixtures. |
| 11 | `codex/implicit-ib-4-transfer-parity-tests-master` | `codex/implicit-ib-5-coupling-aware-vanka-subdomains-master` | Add coupling-aware ASM subdomain support | Add deterministic coupling-aware subdomain construction and hide mutable cache internals. | Coupling-aware ASM construction/application tests. |
| 12 | `codex/implicit-ib-5-coupling-aware-vanka-subdomains-master` | `codex/implicit-ib-6.1-shell-backend-framework-petsc-master` | Add the PETSc shell-smoother backend framework | Add parsed shell-backend selection with the PETSc implementation. | Shared shell-semantics executable and PETSc fixtures. |
| 13 | `codex/implicit-ib-6.1-shell-backend-framework-petsc-master` | `codex/implicit-ib-6.2-eigen-reference-shell-backend-master` | Add the Eigen reference shell backend | Add a serial reference local-solve backend using shared sweep logic. | Eigen-reference shell fixtures. |
| 14 | `codex/implicit-ib-6.2-eigen-reference-shell-backend-master` | `codex/implicit-ib-6.3-eigen-optimized-shell-backends-master` | Add optimized Eigen shell backends | Add retained optimized Eigen dense-solve modes behind the same composition seam. | Optimized Eigen backend fixtures. |
| 15 | `codex/implicit-ib-6.3-eigen-optimized-shell-backends-master` | `codex/implicit-ib-6.4-stokes-eigen-schur-shell-backend-master` | Add the Stokes Eigen-Schur shell backend | Add the Stokes-specialized Schur local solve while reusing shared sweep logic. | Eigen-Schur backend fixtures. |

## Provisional multiplicative-CAV boundaries

| Order | Base branch | Head branch | Proposed title | One claim | Focused native evidence |
|---:|---|---|---|---|---|
| 16 | `codex/implicit-ib-6.4-stokes-eigen-schur-shell-backend-master` | `codex/implicit-ib-cav-0-recovery-plan-master` | Add the multiplicative CAV recovery and audit plan | Commit the authoritative plan, historical ledger, and independent-audit protocol only. | Documentation cross-check; no production test claim. |
| 17 | `codex/implicit-ib-cav-0-recovery-plan-master` | `codex/implicit-ib-cav-1a-live-operator-audit-schema-master` | Expose live Stokes operator state | Add a narrow nonowning live-state view without serialization or solver behavior changes. | `navier_stokes/...shell_multiplicative_semantics_2d.live_operator_state`. |
| 18 | `codex/implicit-ib-cav-1a-live-operator-audit-schema-master` | `codex/implicit-ib-cav-1b-raw-export-comparator-contract-master` | Add the CAV raw export comparator contract | Add test-only raw mapping/comparison and mutation-sensitive controls. | Raw-export contract plus comparator control fixtures. |
| 19 | `codex/implicit-ib-cav-1b-raw-export-comparator-contract-master` | `codex/implicit-ib-cav-1c-matrix-live-kernel-consistency-master` | Align matrix IB kernels with live Fortran operators | Make the matrix builder follow production Fortran coordinate, stencil, kernel, tensor-product, and traversal formulas. | `IBTK/sc_interp_matrix_live_kernel_consistency` in 2D/3D plus affected legacy outputs. |
| 20 | `codex/implicit-ib-cav-1c-matrix-live-kernel-consistency-master` | `codex/implicit-ib-cav-2-patch-construction-audit-master` | Add pressure-seeded CAV patch construction | Add deterministic pressure-seeded RELAXED/STRICT patch construction without application semantics. | 2D/3D CAV reference-parity and transfer-parity construction cases. |
| 21 | `codex/implicit-ib-cav-2-patch-construction-audit-master` | `codex/implicit-ib-cav-3a-blas-lapack-lu-backend-master` | Add the BLAS/LAPACK LU shell backend | Add creator-owned local factors/workspaces, LU solve, lifecycle, and failure propagation. | LU, additive-LU, and singular expected-failure shell fixtures. |
| 22 | `codex/implicit-ib-cav-3a-blas-lapack-lu-backend-master` | `codex/implicit-ib-cav-3b-blas-lapack-robust-modes-master` | Add robust BLAS/LAPACK shell modes | Add the retained robust solve modes and explicit unsupported/failure policy. | Robust, SVD, QR expected-failure, and Cholesky-policy fixtures. |
| 23 | `codex/implicit-ib-cav-3b-blas-lapack-robust-modes-master` | `codex/implicit-ib-cav-4-multiplicative-composition-master` | Pin multiplicative correction-composition semantics | Add exact per-patch state/regression coverage without production behavior changes. | BLAS/LAPACK stagewise fixture. |
| 24 | `codex/implicit-ib-cav-4-multiplicative-composition-master` | `codex/implicit-ib-cav-5-smoother-parity-master` | Add the pressure-seeded multiplicative CAV smoother | Integrate pressure-seeded RELAXED/STRICT CAV into the level/FAC path. | CAV shell RELAXED/STRICT plus live pressure-CAV component fixture. |
| 25 | `codex/implicit-ib-cav-5-smoother-parity-master` | `codex/implicit-ib-cav-6a-fac-stage-observability-master` | Add passive FAC stage observation | Add narrow read-only FAC-cycle observation without changing cycle semantics. | Existing live pressure-CAV component fixture. |
| 26 | `codex/implicit-ib-cav-6a-fac-stage-observability-master` | `codex/implicit-ib-cav-6b-fac-vcycle-fgmres-parity-master` | Keep the outer Stokes-IB Jacobian matrix-free | Preserve the production matrix-free outer operator and cover FAC/regrid behavior. | `IB/implicit_stokes_ib_01.cav_regrid` and live pressure-CAV component cases. |
| 27 | `codex/implicit-ib-cav-6b-fac-vcycle-fgmres-parity-master` | `codex/implicit-ib-cav-6c-shared-correction-composition-master` | Centralize shell correction composition | Give all concrete local backends one backend-neutral traversal/update seam. | Exact stagewise shell fixture across attached backends. |
| 28 | `codex/implicit-ib-cav-6c-shared-correction-composition-master` | `codex/implicit-ib-cav-6d-affected-row-residual-update-master` | Update multiplicative residuals on affected rows | Remove per-patch full-domain residual actions while preserving exact original-residual semantics. | Stagewise shell result and affected-row pattern-invalidation expected failure. |
| 29 | `codex/implicit-ib-cav-6d-affected-row-residual-update-master` | `codex/implicit-ib-cav-6e-blas-lapack-lu-alias-cleanup-master` | Remove the duplicate fixed-LU backend alias | Keep one BLAS/LAPACK backend surface with explicit local-solver selection. | Redirected LU, additive-LU, and singular-LU fixtures. |
| 30 | `codex/implicit-ib-cav-6e-blas-lapack-lu-alias-cleanup-master` | `codex/implicit-ib-cav-6f1-live-numerical-diagnostics-master` | Add live low-stiffness CAV diagnostics | Add native original momentum/divergence/IB residual and consistency reporting. | Pressure-CAV `K=1`, `1e2`, and `1e4` fixtures. |
| 31 | `codex/implicit-ib-cav-6f1-live-numerical-diagnostics-master` | `codex/implicit-ib-cav-6f2a-live-export-schema-master` | Add the live CAV export schema | Add test-only live construction/provenance export and round-trip checks. | Live export and raw comparator fixtures. |
| 32 | `codex/implicit-ib-cav-6f2a-live-export-schema-master` | `codex/implicit-ib-cav-6f2b-live-dynamic-traces-master` | Add live multiplicative CAV traces | Add test-only local solve, FAC, V-cycle, FGMRES, and physical-summary traces. | Existing live pressure-CAV and shell fixtures with enabled/disabled trace paths. |
| 33 | `codex/implicit-ib-cav-6f2b-live-dynamic-traces-master` | `codex/implicit-ib-cav-6f2c-paired-common-arithmetic-replay-master` | Add paired CAV construction and smoother replay | Add external N0/N1 common-arithmetic tooling and mutation controls only. | Committed replay self-checks; external evidence does not replace native `attest`. |
| 34 | `codex/implicit-ib-cav-6f2c-paired-common-arithmetic-replay-master` | `codex/implicit-ib-cav-6f2d-rt0-restriction-ghost-invariant-master` | Fill RT0 restriction ghosts in zero-pre FAC cycles | Satisfy the existing RT0 fine-ghost invariant without changing transfer formulas. | Legal pressure-CAV hierarchy fixture through normal FAC/FGMRES. |
| 35 | `codex/implicit-ib-cav-6f2d-rt0-restriction-ghost-invariant-master` | `codex/implicit-ib-cav-6f2e-multilevel-common-arithmetic-replay-master` | Complete the multiplicative CAV replay matrix | Extend external replay through transfers, V-cycle, FGMRES, and backend cross-checks only. | Committed replay self-checks plus retained native pressure-CAV fixture. |
| 36 | `codex/implicit-ib-cav-6f2e-multilevel-common-arithmetic-replay-master` | `codex/implicit-ib-cav-6f3-production-robustness-master` | Cover the 3D CAV solver lifecycle | Reuse the existing dimensional source for one legal 3D init/deinit/reinit solve case. | `...cav_shell_semantics_3d.pressure_relaxed_blas_lu`. |
| 37 | `codex/implicit-ib-cav-6f3-production-robustness-master` | `codex/implicit-ib-cav-6g1-borrowed-galerkin-operator-master` | Borrow supplied Galerkin operators | Use the standard creator-owned nonowning `Mat` lifetime on the nonaugmented FAC path. | Supplied-operator lifetime and Galerkin-borrowing fixtures. |
| 38 | `codex/implicit-ib-cav-6g1-borrowed-galerkin-operator-master` | `codex/implicit-ib-cav-6g2-lazy-velocity-seed-map-master` | Build velocity-seed pairs on demand | Retain the legacy velocity STRICT map only where it is consumed. | 2D/3D construction absence/count/order outputs. |
| 39 | `codex/implicit-ib-cav-6g2-lazy-velocity-seed-map-master` | `codex/implicit-ib-cav-6g3-absent-unrestricted-partition-master` | Omit unrestricted ownership partitions | Represent pressure-CAV unrestricted ownership as absent across supported backends. | CAV shell representation/lifecycle and live pressure-CAV fixtures. |
| 40 | `codex/implicit-ib-cav-6g3-absent-unrestricted-partition-master` | `codex/implicit-ib-cav-6g4-fac-residual-vector-cache-production-master` | Cache FAC residual work vectors | Allocate hierarchy-dependent PETSc residual work vectors at setup and reuse them during apply. | Rediscretized and Galerkin lifecycle-count/repeated-residual fixtures. |

## Post-production documentation tips

The local `codex/implicit-ib-cav-6g4-fac-residual-vector-cache-master` tip adds the candidate-finalization/audit-handoff commits through the patch-equivalent cumulative rehearsal. The local `codex/implicit-ib-cav-master-merge-preparation` branch carries later frozen M5, M6, branch-mapping, and merge-readiness records. These tips change no production source beyond order 40. Submit or retain them as explicitly documentation/evidence-only changes; do not combine them with a new semantic fix.

## Future PR test block

Adapt the focused selector to the exact unit above and record the exact checkout path and discovered cases:

```text
Focused native attest
- Command: <checkout>/attest --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R <focused-selector>
- Result: <passed>/<discovered>
- Deterministic rerun: <result>

Broad IB regression
- Command: <checkout>/attest --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R IB/
- Result: <passed>/<discovered>
- Exclusions: <none or explicit record>
```

No timing, speedup, throughput, or benchmark claim belongs in these PR bodies. Test duration may be recorded only as CI operational metadata.
