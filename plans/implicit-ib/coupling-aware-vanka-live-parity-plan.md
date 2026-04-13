# Coupling-Aware Vanka Live Parity Plan

Date: 2026-04-13

## Objective
Identify the first concrete divergence between live IBAMR (`examples/IB/implicit/ex0`) and live MATLAB (`/Users/boyceg/code/implicit_ib_coupling_aware_vanka`) for the coupling-aware Vanka multigrid preconditioner, for both `RELAXED` and `STRICT` closure policies.

## Non-Negotiable Audit Requirements
These requirements were missing or inconsistently enforced in earlier parity attempts and are now treated as mandatory:
1. Use live objects/solvers only (no surrogate operators or offline stand-ins).
2. Require exact Stage-B subdomain semantics:
   - subdomain sequence/order must match exactly (seed/order correspondence),
   - overlap/nonoverlap DOF membership for each corresponding subdomain is set-equality (intra-subdomain DOF order is non-semantic),
   - resulting local block matrices must match exactly after the established DOF/sign normalization.
3. Require exact transfer-operator parity (`P`, `R`) across levels.
4. Require coarse-solver semantics parity before interpreting Stage-D discrepancies.
5. Apply explicit MATLAB->IBAMR DOF reindexing from metadata (never assume index ranges align).
6. Handle pressure gauge correctly (gauge-projected comparisons where pressure nullspace exists).
7. Enforce known Stokes sign-convention handling for coupled blocks (including the `A10` mismatch: IBAMR `A10 = A10' = -Div` vs MATLAB `A10 = -A10' = Div`) and interpret Stage-D relative differences with `cond2(A00)`-aware expectations.

## Audit Snapshot (Repro Baseline)
- Time: 2026-04-12.
- Driver: `examples/IB/implicit/ex0/main2d` (serial, `np=1`).
- Finest-grid requirement: `N_FINE >= 32` (enforced by the audit runner).
- SAJ path for retained evidence: matrix-free (`USE_MATRIX_BASED_SAJ=FALSE`, i.e. default unless explicitly overridden).
- Stage-D comparison mode for retained sweep evidence: pressure-gauge-projected (`--stage-d-pressure-gauge-projected=true`).
- Comparator tolerance used in retained evidence: `rel_tol=1e-10`.

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
  - Stage-A metadata parity is required (`case_id`, `closure_policy`, `finest_level`, `num_curve_points`, `dt`, `rho`, `mu`, `marker_spacing_ds`, `level_global_dof_counts`),
  - Stage-D pass/fail uses gauge-projected final output metrics when pressure DOFs exist,
  - Stage-D coarse semantics check is required (`--stage-d-require-coarse-semantics=true`):
    - coarse RHS compares `preconditioned_apply_coarse_rhs_level0.mtx` after MATLAB pressure-row sign normalization,
    - coarse correction compares `preconditioned_apply_coarse_correction_level0.mtx` after pressure-gauge projection,
    - coarse semantic comparisons gate Stage-D pass/fail together with final-output acceptance,
  - Stage-D smoother-level vector parity is required in strict profile (`--stage-d-require-smoother-semantics=true`),
  - Stage-D strict MATLAB `cond2(A00)` backend can be enforced (`--stage-d-cond2-strict-backend=true`, fail instead of python fallback),
  - Stage B coarsest-level report note is `skipped subdomain checks (ibamr seeds=..., matlab seeds=...); coarsest-level solve is direct semantics and is gated in Stage D coarse semantic vectors`,
  - Stage-D smoother-level vectors (`preconditioned_apply_{pre,post}_smooth_*`) can be diagnostic-only in non-strict profile (`--stage-d-require-smoother-semantics=false`),
  - `--stage-d-pressure-gauge-projected` controls smoother-diagnostic metric presentation.
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
  - `OUTPUT_ROOT=/tmp/ibamr_cav_k_sweep CASE_ID=ex0_k<value> ELASTIC_K=<value> STAGE_D_PRESSURE_GAUGE_PROJECTED=true STAGE_D_REQUIRE_COARSE_SEMANTICS=true STAGE_D_REQUIRE_SMOOTHER_SEMANTICS=true STAGE_D_COND_SAFETY_FACTOR=5.0 STAGE_D_COND_THRESHOLD_FLOOR=1.0e-12 STAGE_D_FIXED_THRESHOLD=1.0e-10 STAGE_D_COND2_BACKEND=matlab STAGE_D_COND2_STRICT_BACKEND=true /Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_audit.sh`
- Comparator direct run:
  - `python3 /Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_compare.py --ibamr-dir <IB_BUNDLE> --matlab-dir <MAT_BUNDLE> --max-stage D --rel-tol 1e-10 --normalize-pressure-row-sign true --stage-d-pressure-gauge-projected true --stage-d-require-coarse-semantics true --stage-d-require-smoother-semantics true --stage-d-cond-safety-factor 5.0 --stage-d-cond-threshold-floor 1e-12 --stage-d-fixed-threshold 1e-10 --stage-d-cond2-backend matlab --stage-d-cond2-strict-backend true`
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

### 2026-04-12: Superseded Stage-D runs archived
- Earlier Stage-D sweeps using older acceptance semantics (`fixed only`, then `fixed + conditioning` with hard smoother gating) are superseded by the composite acceptance policy below.
- Superseded artifact roots kept for traceability:
  - `/tmp/ibamr_cav_k_sweep`
  - `/tmp/ibamr_cav_k_sweep_matlab_cond2_floor1e12_20260412`

### 2026-04-12: Full K sweep with updated composite Stage-D acceptance
- Acceptance rule applied in comparator:
  - Gauge-projected final output metric only (`preconditioned_apply_output_level_fine`),
  - `threshold = max(1e-10, C * cond2(A00) * eps_machine, 1e-12)`,
  - `C = 5`,
  - smoother-level diagnostics are warnings (non-gating).
- Configuration:
  - `STAGE_D_COND2_BACKEND=matlab`
  - `STAGE_D_COND_SAFETY_FACTOR=5.0`
  - `STAGE_D_FIXED_THRESHOLD=1.0e-10`
  - `STAGE_D_COND_THRESHOLD_FLOOR=1.0e-12`
  - `STAGE_D_PRESSURE_GAUGE_PROJECTED=true`
- Artifact root:
  - `/tmp/ibamr_cav_k_sweep_acceptance_relaxed_criteria_20260412`
- Results by `K` and policy:
  - `K=1e0`
    - `RELAXED`: Stage D PASS (`composite_threshold=1.0e-10`, `rel_linf=1.509e-15`, `rel_l2=9.514e-16`).
    - `STRICT`: Stage D PASS (`composite_threshold=1.0e-10`, `rel_linf=3.137e-15`, `rel_l2=2.131e-15`).
  - `K=1e2`
    - `RELAXED`: Stage D PASS (`composite_threshold=1.0e-10`, `rel_linf=3.300e-14`, `rel_l2=7.621e-15`).
    - `STRICT`: Stage D PASS (`composite_threshold=1.0e-10`, `rel_linf=3.156e-14`, `rel_l2=6.559e-15`).
  - `K=1e4`
    - `RELAXED`: Stage D PASS (`composite_threshold=1.0e-10`, `rel_linf=2.048e-12`, `rel_l2=2.108e-12`).
    - `STRICT`: Stage D PASS (`composite_threshold=1.0e-10`, `rel_linf=1.289e-11`, `rel_l2=1.379e-11`).
  - `K=1e6`
    - `RELAXED`: Stage D PASS (`composite_threshold=1.913e-09`, `rel_linf=5.461e-11`, `rel_l2=3.509e-11`); smoother-level mismatch retained as diagnostic (`level 2 pre output`).
    - `STRICT`: Stage D PASS (`composite_threshold=1.913e-09`, `rel_linf=1.154e-10`, `rel_l2=1.153e-10`); smoother-level mismatch retained as diagnostic (`level 1 pre output`).
  - `K=1e8`
    - `RELAXED`: Stage D FAIL (`composite_threshold=1.913e-07`, `rel_linf=2.696e-06`, `rel_l2=8.303e-07`).
    - `STRICT`: Stage D FAIL (`composite_threshold=1.913e-07`, `rel_linf=1.318e+00`, `rel_l2=4.000e-01`).

### 2026-04-12: Comparator semantics tightened for Requirement #4 (coarse solver semantics)
- Comparator update (`run_ex0_cav_live_parity_compare.py`):
  - Stage A/B/C/smoother/coarse-semantic checks use relative tolerances only,
  - Stage D now enforces coarse-level semantic checks (default-on via `--stage-d-require-coarse-semantics=true`),
  - coarse RHS check compares `preconditioned_apply_coarse_rhs_level0.mtx` after MATLAB pressure-row sign normalization,
  - coarse correction check compares `preconditioned_apply_coarse_correction_level0.mtx` after pressure-gauge projection.
- Audit runner update (`run_ex0_cav_live_parity_audit.sh`):
  - forwards `STAGE_D_REQUIRE_COARSE_SEMANTICS` to comparator (default `true`).
- Replay on retained sweep artifacts (`/tmp/ibamr_cav_k_sweep_acceptance_relaxed_criteria_20260412`, using `--stage-d-cond2-backend python` for quick confirmation):
  - `K=1e0,1e2,1e4,1e6`: `RELAXED` and `STRICT` remain PASS.
  - `K=1e8`: `RELAXED` and `STRICT` remain FAIL.

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
- For Stage-D preconditioner-apply comparisons, use final-output gauge-projected parity with a composite threshold and require coarse semantic vector parity.
- Recommended criterion (gauge-projected metric):
  - `rel_err <= max(1e-10, C * cond2(A00) * eps_machine, 1e-12)`
- Where:
  - `rel_err` is the reported gauge-projected relative mismatch (`rel_linf` or `rel_l2`),
  - `cond2(A00)` is the finest-level velocity-block condition number for the tested `K`,
  - `eps_machine` is IEEE double precision epsilon,
  - `C` is a modest safety factor for implementation details and accumulation (`C=5` in current retained runs).
- This criterion gates:
  - final preconditioned-apply output vectors (pressure-gauge-projected when pressure DOFs exist),
  - coarse RHS semantics (`preconditioned_apply_coarse_rhs_level0`) after pressure-row sign normalization,
  - coarse correction semantics (`preconditioned_apply_coarse_correction_level0`) after pressure-gauge projection.
- In non-strict profile, smoother-level diagnostics remain informative but non-gating.
- In strict semantic profile (`--stage-d-require-smoother-semantics=true`), smoother-level diagnostics are promoted to gating semantics.

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
  - `stage_d_outcome_composite_acceptance`,
  - `stage_d_outcome_with_coarse_semantics`,
  - `coarse_level_semantics` (raw + normalized coarse RHS/coarse correction diagnostics),
  - detailed `cond2` diagnostics (`a00_shape`, `a00_nnz`, backend metadata, singular-value/estimator metadata).
- Console/text report Stage-D line now prints both outcomes and the live `cond2(A00)` threshold used.
- Conditioning-aware threshold now applies a floor:
  - `threshold = max(C * cond2(A00) * eps_machine, 1e-12)`.
- Composite Stage-D acceptance threshold uses:
  - `max(1e-10, C * cond2(A00) * eps_machine, 1e-12)`.

### 2026-04-12: Stage-D reporting semantics (historical non-strict profile)
- Historical non-strict comparator profile used:
  - composite acceptance on final output-vector gauge-projected metrics,
  - coarse semantic vector acceptance when `--stage-d-require-coarse-semantics=true`,
  - smoother-level diagnostics as non-gating warnings.
- This profile is retained for traceability only and is superseded for acceptance decisions by the strict semantic profile documented below.

### 2026-04-13: Parity-Audit Semantics Hardening Implementation
- Comparator updates (`run_ex0_cav_live_parity_compare.py`):
  - Stage-B pass message now states explicit semantics:
    - subdomain sequence/order is enforced,
    - overlap/nonoverlap DOFs are compared as sets per subdomain.
  - Stage-B now includes default-on IBAMR first-sweep consistency validation for smoother levels when first-sweep artifacts are present:
    - `subdomains_first_sweep_level*` vs `subdomains_level*`,
    - `A_subdomain_first_sweep_level*_k*` vs `A_subdomain_level*_k*`.
  - Missing-file/parse failures at stage boundaries are converted into structured Stage-level `FAIL` messages with explicit file paths (no Python traceback in normal report flow).
  - Report metadata now records:
    - IBAMR git SHA,
    - MATLAB reference repo SHA,
    - MATLAB `v_cycle.m` override usage/path/SHA256.
- Runner updates (`run_ex0_cav_live_parity_audit.sh`):
  - Verifies required bundle artifacts after each IBAMR/MATLAB export pair and before comparator invocation.
  - On missing artifacts, emits clear policy/stage failure with explicit log pointers (`run_ibamr.stdout`, `run_matlab.stdout`) and writes structured report files.
  - Forwards reproducibility metadata fields to comparator reports.
- Re-run outcomes with hardened scripts:
  - `K=1e6`: `RELAXED=PASS`, `STRICT=PASS` (Stage-D composite threshold `~1.913e-09`).
  - `K=1e8`: `RELAXED=FAIL`, `STRICT=FAIL` (Stage-D composite threshold `~1.913e-07`).
  - Stage-B reports retain the coarsest-level skip note plus explicit coarse-semantics gating note.

### 2026-04-13: Skeptical Re-Audit Follow-Up Hardening
- Additional comparator hardening (`run_ex0_cav_live_parity_compare.py`):
  - Stage-D audits now **require** first-sweep consistency artifacts on all smoother levels (instead of treating missing first-sweep artifacts as a skip).
  - Stage-D overall pass/fail now requires preconditioned-apply **input vector parity** (`preconditioned_apply_input_level_fine.mtx`) in addition to output/composite/coarse-semantic checks.
  - Stage-D failure message now reports explicit mismatch clauses (`input`, `output`, `coarse_semantics`, `smoother_semantics`) to reduce ambiguity in root-cause interpretation.
  - Report metadata now records whether Stage-B first-sweep consistency was required (`stage_b_require_first_sweep_consistency`).
- Additional runner hardening (`run_ex0_cav_live_parity_audit.sh`):
  - Stage-D artifact validation now requires IBAMR first-sweep files on smoother levels:
    - `subdomains_first_sweep_level*`,
    - `A_subdomain_first_sweep_level*_k*`.
- MATLAB wrapper reproducibility hardening (`run_ex0_cav_live_parity_export_matlab.m`):
  - uses per-run unique temporary scratch directory (no shared `tempdir` filename collisions),
  - resets MATLAB path via `restoredefaultpath()` before adding only:
    - parity scratch override directory (`v_cycle.m`, `-begin`),
    - MATLAB reference repository root (`-end`),
  - restores original path on exit and cleans temporary directory.
- Fresh live re-runs (outside sandbox) used for this follow-up:
  - `/tmp/ibamr_cav_live_parity_reaudit_20260413_k1e6`
  - `/tmp/ibamr_cav_live_parity_reaudit_20260413_k1e8`
  - `/tmp/ibamr_cav_live_parity_reaudit_20260413_k1e6_hardened`
  - `/tmp/ibamr_cav_live_parity_reaudit_20260413_k1e8_hardened`
- Observed outcomes in this non-strict profile:
  - `K=1e6`: `RELAXED=PASS`, `STRICT=PASS`.
  - `K=1e8`: `RELAXED=FAIL`, `STRICT=FAIL`.

### 2026-04-13: Strict semantic profile hardening and re-run
- Comparator updates (`run_ex0_cav_live_parity_compare.py`):
  - Stage-A now gates on core metadata parity (case/policy, level counts, `dt`, `rho`, `mu`, `ds`) before marker checks.
  - Added strict smoother semantics gate (`--stage-d-require-smoother-semantics=true`).
  - Added strict MATLAB cond2 backend mode (`--stage-d-cond2-strict-backend=true`) to fail instead of falling back to python.
  - Stage-D overall outcome now requires all enabled semantic clauses: input, output composite acceptance, coarse semantics, smoother semantics.
- Runner updates (`run_ex0_cav_live_parity_audit.sh`):
  - forwards `STAGE_D_REQUIRE_SMOOTHER_SEMANTICS` (default `true`),
  - forwards `STAGE_D_COND2_STRICT_BACKEND` (default `true`).
- Fresh strict-profile full sweep artifact root:
  - `/tmp/ibamr_cav_live_parity_reaudit_codex_20260413_strictprofile`
- Strict-profile outcomes:
  - `K=1e0`: `RELAXED=PASS`, `STRICT=PASS`.
  - `K=1e2`: `RELAXED=PASS`, `STRICT=PASS`.
  - `K=1e4`: `RELAXED=PASS`, `STRICT=PASS`.
  - `K=1e6`: `RELAXED=FAIL`, `STRICT=FAIL` (smoother semantics mismatch only; output composite and coarse semantics still pass).
  - `K=1e8`: `RELAXED=FAIL`, `STRICT=FAIL` (output composite, coarse semantics, and smoother semantics all fail).

### 2026-04-13: Rebuild + runner stale-artifact hardening + live rerun
- Rebuilt target:
  - `CCACHE_DIR=/tmp/ccache CCACHE_TEMPDIR=/tmp/ccache/tmp cmake --build /Users/boyceg/code/ibamr-objs-opt --target IB-implicit-ex0 -j8`
- Comparator diagnostics hardening (`run_ex0_cav_live_parity_compare.py`):
  - when MATLAB `cond2(A00)` backend fails, Stage-D failure now includes `stdout_tail`, `stderr_tail`, and any partial `cond2_out.txt` tail.
- Audit runner correctness hardening (`run_ex0_cav_live_parity_audit.sh`):
  - clear per-policy IBAMR/MATLAB bundle directories before Stage C and Stage D exports,
  - fail explicitly with structured report output if IBAMR export command fails,
  - fail explicitly with structured report output if MATLAB export command fails,
  - mirror Stage-C export failures into the policy-level report files (so stale prior reports are not reused).
- Root-cause found for ambiguous in-sandbox reruns:
  - PETSc/MPI initialization fails under sandbox (`getdomainname(): Operation not permitted`), so sandbox runs are not valid for live parity evidence.
  - Fresh evidence must be collected outside sandbox.
- Fresh post-rebuild live artifact roots (outside sandbox, strict profile defaults):
  - `/tmp/ibamr_cav_live_parity_rebuild_20260413_k1e6`
  - `/tmp/ibamr_cav_live_parity_rebuild_20260413_k1e8`
- Post-rebuild outcomes:
  - `K=1e6`:
    - `RELAXED=FAIL`, `STRICT=FAIL`.
    - Failure clause is smoother semantics only (`smoother_semantics_ok=false`), while Stage-D output composite acceptance and coarse semantics are satisfied.
    - Example strict mismatch: `level 1 pre output rel_ref(linf=2.479e-10, l2=2.231e-10)`.
  - `K=1e8`:
    - `RELAXED=FAIL`, `STRICT=FAIL`.
    - Failure clauses include output composite acceptance, coarse semantics, and smoother semantics.
    - `cond2(A00)≈1.722844e+08`, Stage-D composite threshold `≈1.913e-07`.

## Next Action
- Keep Stage A/B/C and metadata parity locked.
- Debug and resolve Stage-D smoother semantic divergence at `K=1e6` (first strict-profile semantic failure).
- For `K=1e6` root-cause:
  - compare level-1 and level-2 pre-smoother vectors (`preconditioned_apply_pre_smooth_output_level{1,2}.mtx`) against first-sweep local-subdomain solve outputs to localize whether divergence starts in subdomain solve, ordering, or accumulation,
  - replay with same bundles and stricter diagnostics around per-subdomain contribution accumulation order,
  - keep using live object exports only (no surrogate operators).
- After `K=1e6` smoother parity is fixed, re-evaluate whether `K=1e8` should be treated as an in-scope parity target or high-conditioning out-of-scope regime.
