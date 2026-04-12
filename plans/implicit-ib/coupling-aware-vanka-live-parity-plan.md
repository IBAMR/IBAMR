# Coupling-Aware Vanka Live Parity Plan

Date: 2026-04-12

## Objective
Identify the first concrete divergence between live IBAMR (`examples/IB/implicit/ex0`) and live MATLAB (`/Users/boyceg/code/implicit_ib_coupling_aware_vanka`) for the coupling-aware Vanka multigrid preconditioner, for both `RELAXED` and `STRICT` closure policies.

## Audit Snapshot (Repro Baseline)
- Time: 2026-04-12.
- Driver: `examples/IB/implicit/ex0/main2d` (serial, `np=1`).
- Finest-grid requirement: `N_FINE >= 32` (enforced by the audit runner).
- SAJ path for retained evidence: matrix-free (`USE_MATRIX_BASED_SAJ=FALSE`, i.e. default unless explicitly overridden).
- Stage-D comparison mode for retained sweep evidence: pressure-gauge-projected (`--stage-d-pressure-gauge-projected=true`).
- Comparator tolerances used in retained evidence: `abs_tol=1e-12`, `rel_tol=1e-10`.

## Artifact Contract
- Bundle format: MatrixMarket + JSON.
- Primary sweep artifact root: `/tmp/ibamr_cav_k_sweep`.
- Comparator report roots:
  - `/tmp/ibamr_cav_k_sweep/k*/reports`
  - `/tmp/ibamr_cav_live_parity/reports`

## Index And Sign Conventions
- State block order is unchanged in both implementations: `[u; v; p]`.
- MATLAB parity DOF metadata uses `dof = j + N*i` within each field block (`reshape(..., ny, nx)`, column-major, `y` fastest).
- IBAMR DOF matching uses metadata (`kind`, `axis`, `x`, `y`) instead of raw integer ranges.
- Known Stokes sign difference (Stage B matrix comparisons):
  - IBAMR lower-left block uses `+A01'`.
  - MATLAB lower-left block uses `-A01'`.
  - Comparator normalizes this with MATLAB pressure-row sign flip in Stage B.
- Stage D smoother/vector comparisons:
  - always apply MATLAB->IBAMR DOF permutation,
  - no pressure-row sign flip,
  - no marker-force `ds` scaling,
  - pass criterion uses pressure-gauge-projected metrics by default (`--stage-d-pressure-gauge-projected=true`), while raw metrics are still reported.
- Coupling-aware seed traversal naming:
  - 2D: `I_J`, `J_I`.
  - 3D: `I_J_K`, `J_K_I`, `K_I_J`.
  - `ex0` parity runs use 2D `I_J` to match MATLAB traversal.

## Stage Checklist (Current)

### RELAXED
- [x] Stage A: Marker count/positions/forces parity.
- [x] Stage B: Subdomain order/membership + level/subdomain matrix parity.
- [x] Stage C: Restriction/prolongation operator parity.
- [ ] Stage D: Finest-level preconditioned apply parity (`y = M^{-1}A x`) for all tested `K`.
- [x] Earliest mismatch recorded.

### STRICT
- [x] Stage A: Marker count/positions/forces parity.
- [x] Stage B: Subdomain order/membership + level/subdomain matrix parity.
- [x] Stage C: Restriction/prolongation operator parity.
- [ ] Stage D: Finest-level preconditioned apply parity.
- [x] Earliest mismatch recorded.

## Commands
- Build:
  - `CCACHE_DIR=/tmp/ccache CCACHE_TEMPDIR=/tmp/ccache/tmp cmake --build /Users/boyceg/code/ibamr-objs-opt --target IB-implicit-ex0 -j8`
- Full parity audit (both policies, one `K`):
  - `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_audit.sh`
- K sweep (gauge-projected Stage D):
  - `OUTPUT_ROOT=/tmp/ibamr_cav_k_sweep CASE_ID=ex0_k<value> ELASTIC_K=<value> STAGE_D_PRESSURE_GAUGE_PROJECTED=true STAGE_D_COND_SAFETY_FACTOR=1.0 STAGE_D_COND_THRESHOLD_FLOOR=1.0e-12 STAGE_D_COND2_BACKEND=matlab /Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_audit.sh`
- Comparator direct run:
  - `python3 /Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_compare.py --ibamr-dir <IB_BUNDLE> --matlab-dir <MAT_BUNDLE> --max-stage D --rel-tol 1e-10 --abs-tol 1e-12 --normalize-pressure-row-sign true --stage-d-pressure-gauge-projected true --stage-d-cond-safety-factor 1.0`
  - Audit runner control for conditioning-aware threshold factor:
    - `STAGE_D_COND_SAFETY_FACTOR=<C> /Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_audit.sh`

## Retained Evidence (Current)

### 2026-04-12: Stage A/B/C parity is established
- `RELAXED` and `STRICT` pass Stages A/B/C under the current comparator.
- Marker-force parity uses MATLAB force scaling by `ds` (`F_matlab * ds` vs IBAMR marker force units).

### 2026-04-12: First-sweep ASM state matches Stage-B exports
- Live first-sweep `getASMSubdomains()` and extracted local blocks match Stage-B exported subdomains/blocks exactly on smoother levels.
- Artifact root:
  - `/tmp/ibamr_first_sweep_vs_stageb/ibamr/first_sweep_check/policy=RELAXED`

### 2026-04-12: Backend/subsolver sensitivity does not explain Stage D
- Changing shell backend / subdomain solver options does not materially change IBAMR smoother outputs (machine-precision internal deltas).
- Artifact root:
  - `/tmp/ibamr_smoother_backend_impact2/ibamr/backend_subsolver_impact`

### 2026-04-12: Gauge-projected Stage-D K sweep (historical pre-floor / pre-MATLAB-cond2)
- Sweep values: `K = 1e0, 1e2, 1e4, 1e6, 1e8`.
- `RELAXED` (`rel_tol=1e-10`):
  - `K=1e0`: PASS
  - `K=1e2`: PASS
  - `K=1e4`: PASS
  - `K=1e6`: FAIL (`level 2 pre output`, `rel_linf=1.731e-10`, `rel_l2=1.401e-10`) near threshold.
  - `K=1e8`: FAIL (`level 1 pre output`, `rel_linf=2.064e-08`, `rel_l2=2.023e-08`).
- `STRICT`: this sweep snapshot is stale after strict-policy Stage-D path alignment updates; do not treat this row as current STRICT status.
- Artifact root:
  - `/tmp/ibamr_cav_k_sweep`

### 2026-04-12: Full K sweep with MATLAB cond2 backend + threshold floor
- Configuration:
  - `STAGE_D_COND2_BACKEND=matlab`
  - `STAGE_D_COND_THRESHOLD_FLOOR=1.0e-12`
  - `STAGE_D_COND_SAFETY_FACTOR=1.0`
  - `STAGE_D_PRESSURE_GAUGE_PROJECTED=true`
- Artifact root:
  - `/tmp/ibamr_cav_k_sweep_matlab_cond2_floor1e12_20260412`
- Summary by `K` and policy (fixed/conditioning outcomes are Stage-D output-vector metrics; overall Stage-D status can still fail earlier on smoother-level checks):
  - `K=1e0`
    - `RELAXED`: Stage D PASS; fixed PASS; conditioning PASS; `cond2(A00)=3.051551e+00`, threshold `1.0e-12` (floor active).
    - `STRICT`: Stage D PASS; fixed PASS; conditioning PASS; `cond2(A00)=3.051551e+00`, threshold `1.0e-12` (floor active).
  - `K=1e2`
    - `RELAXED`: Stage D PASS; fixed PASS; conditioning PASS; `cond2(A00)=1.736043e+02`, threshold `1.0e-12` (floor active).
    - `STRICT`: Stage D PASS; fixed PASS; conditioning PASS; `cond2(A00)=1.736043e+02`, threshold `1.0e-12` (floor active).
  - `K=1e4`
    - `RELAXED`: Stage D PASS; fixed PASS; conditioning PASS; `cond2(A00)=1.722976e+04`, threshold `3.825775e-12`.
    - `STRICT`: Stage D PASS; fixed PASS; conditioning FAIL; `cond2(A00)=1.722976e+04`, threshold `3.825775e-12`.
  - `K=1e6`
    - `RELAXED`: Stage D FAIL due smoother gate (`level 2 pre output`, `rel_linf=1.731e-10`, `rel_l2=1.401e-10`); output-vector fixed PASS, conditioning PASS; `cond2(A00)=1.722845e+06`, threshold `3.825485e-10`.
    - `STRICT`: Stage D FAIL due smoother gate (`level 1 pre output`, `rel_linf=2.479e-10`, `rel_l2=2.231e-10`); output-vector fixed FAIL, conditioning PASS; `cond2(A00)=1.722845e+06`, threshold `3.825485e-10`.
  - `K=1e8`
    - `RELAXED`: Stage D FAIL due smoother gate (`level 1 pre output`, `rel_linf=2.064e-08`, `rel_l2=2.023e-08`); output-vector fixed FAIL, conditioning FAIL; `cond2(A00)=1.722844e+08`, threshold `3.825482e-08`.
    - `STRICT`: Stage D FAIL due smoother gate (`level 1 pre output`, `rel_linf=4.886e-09`, `rel_l2=6.565e-09`); output-vector fixed FAIL, conditioning FAIL; `cond2(A00)=1.722844e+08`, threshold `3.825482e-08`.

### 2026-04-12: IBAMR-only regression guard tests added (not MATLAB parity claims)
- Purpose: lock down live operator/preconditioner/preconditioned-operator residual behavior for representative stiffness values in attest.
- Test family:
  - `tests/IB/implicit_stokes_ib_operator_preconditioner_residual_sweep_01.cpp`
  - `tests/IB/implicit_stokes_ib_operator_preconditioner_residual_sweep_01.K=1e2.input`
  - `tests/IB/implicit_stokes_ib_operator_preconditioner_residual_sweep_01.K=1e4.input`
  - `tests/IB/implicit_stokes_ib_operator_preconditioner_residual_sweep_01.K=1e6.input`
  - `tests/IB/implicit_stokes_ib_operator_preconditioner_residual_sweep_01.K=1e8.input`
- Current test configuration:
  - depth fixed by `MAX_LEVELS` (no separate `TEST_DEPTH`),
  - scalar `K` per input file,
  - coupling-aware multiplicative shell smoother with `shell_pc_type = "multiplicative-eigen-pseudoinverse"`.
- Verification command:
  - `/Users/boyceg/code/IBAMR/attest --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R 'IB/implicit_stokes_ib_operator_preconditioner_residual_sweep_01.K='`
- Status: all four K-labeled tests pass.

### 2026-04-12: Conditioning interpretation for Stage-D limits
- Expected forward-error floor for stable solves scales as `O(cond2(A00) * eps_machine)` (not `eps_machine - 1/cond2(A00)`).
- Measured `cond2(A00)` trend on finest level (RELAXED):
  - `K=1e0`: `cond2≈3.051551e+00`, `cond2*eps≈6.775805e-16`
  - `K=1e2`: `cond2≈1.736043e+02`, `cond2*eps≈3.854790e-14`
  - `K=1e4`: `cond2≈1.722976e+04`, `cond2*eps≈3.825775e-12`
  - `K=1e6`: `cond2≈1.722845e+06`, `cond2*eps≈3.825485e-10`
  - `K=1e8`: `cond2≈1.722844e+08`, `cond2*eps≈3.825482e-08`
- Observed Stage-D mismatch scaling is consistent with this conditioning trend.

### 2026-04-12: Recommended Stage-D acceptance criterion (preconditioner apply parity)
- For Stage-D preconditioner-apply comparisons, use a conditioning-aware relative threshold instead of a fixed global `1e-10` cutoff.
- Recommended criterion (gauge-projected metric):
  - `rel_err <= C * cond2(A00) * eps_machine`
- Where:
  - `rel_err` is the reported gauge-projected relative mismatch (`rel_linf` or `rel_l2`),
  - `cond2(A00)` is the finest-level velocity-block condition number for the tested `K`,
  - `eps_machine` is IEEE double precision epsilon,
  - `C` is a modest safety factor for implementation details and accumulation (start with `C=1`, increase only if justified by reproducible evidence).
- This criterion is intended for RELAXED-path parity decisions and should be recorded together with the exact `cond2(A00)` value used for each audited case.

### 2026-04-12: Comparator/report update for dual Stage-D outcomes
- Updated `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_compare.py` so Stage D now emits both:
  - fixed-threshold result (`rel_tol`-based), and
  - conditioning-aware result (`rel_err <= C * cond2(A00) * eps_machine`).
- `cond2(A00)` is now computed live from the exported finest-level `A_level*.mtx` and `dof_map_level*.json` for the run case (velocity block extraction), with selectable backends:
  - `matlab` (default, uses sparse `svds`),
  - `python` (fallback estimator).
- Stage-D report JSON now includes:
  - `stage_d_outcome_fixed_threshold`,
  - `stage_d_outcome_conditioning_aware`,
  - detailed `cond2` diagnostics (`a00_shape`, `a00_nnz`, backend metadata, singular-value/estimator metadata).
- Console/text report Stage-D line now prints both outcomes and the live `cond2(A00)` threshold used.
- Conditioning-aware threshold now applies a floor:
  - `threshold = max(C * cond2(A00) * eps_machine, 1e-12)`.

### 2026-04-12: Post-STRICT-path-fix live check (`K=1e2`)
- Command:
  - `OUTPUT_ROOT=/tmp/ibamr_cav_cond_report_check_escalated CASE_ID=ex0_cond_report_check_escalated ELASTIC_K=1.0e2 /Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_audit.sh`
- Artifact root:
  - `/tmp/ibamr_cav_cond_report_check_escalated`
- RELAXED Stage D:
  - fixed-threshold: PASS (`rel_linf_ref=3.2996489153587185e-14`, `rel_l2_ref=7.620725508442902e-15`, `rel_tol=1e-10`)
  - conditioning-aware: PASS (`cond2(A00)=173.58078578977316`, `threshold=3.854267700326667e-14`, `C=1.0`)
- STRICT Stage D:
  - fixed-threshold: PASS (`rel_linf_ref=3.156336826397238e-14`, `rel_l2_ref=6.558693025970851e-15`, `rel_tol=1e-10`)
  - conditioning-aware: PASS (`cond2(A00)=173.58078578977316`, `threshold=3.854267700326667e-14`, `C=1.0`)

### 2026-04-12: Stage-D reporting semantics
- Comparator top-level pass/fail status remains fixed-threshold-driven for compatibility.
- Conditioning-aware status is emitted alongside fixed-threshold status in Stage-D details and should be used for conditioning-aware acceptance decisions.
- Stage-D can fail before output-vector acceptance due smoother-level diagnostics (`preconditioned_apply_{pre,post}_smooth_*`); in those cases, output-vector fixed/conditioning outcomes are still emitted in details when available.

## Next Action
- Keep Stage A/B/C locked.
- Use the emitted dual Stage-D outcomes to classify each run under both criteria (fixed and conditioning-aware) before deeper mismatch localization.
