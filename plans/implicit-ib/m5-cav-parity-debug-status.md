# MS5 CAV Parity Debug Status

Date: 2026-03-14
Scope: 2D, serial, CAV parity against MATLAB/Octave and Julia reference workflows.

Current audit expansion plan:
- [m5-expanded-julia-audit-matrix.md](/Users/boyceg/code/IBAMR/plans/implicit-ib/m5-expanded-julia-audit-matrix.md)

## Test Drivers and Commands

### 1) IBAMR export + bridge parity orchestrator
Script:
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/run_ibamr_cav_subdomain_bridge_parity.sh`

Typical command form:
- `tests/external-reference/cav-stokes-ib/reference-tools/run_ibamr_cav_subdomain_bridge_parity.sh /Users/boyceg/code/ibamr-objs-dbg <work_dir> /Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.multilevel_l2.input <K> <mode> <N>`

Modes used:
- `one_level`
- `two_level`

K values used in this pass:
- `1e2`
- `1e4`
- `1e6`

N values used in this pass:
- `16`
- `8`

### 2) Octave comparator (direct)
Script:
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_ibamr_exported_subdomains_octave.m`

Typical command form:
- `/opt/homebrew/bin/octave --quiet --eval "addpath('/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools'); compare_ibamr_exported_subdomains_octave('<export_dir>', '<matlab_out_dir>', <tol>);"`

### 3) Julia comparator (direct)
Script:
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_ibamr_exported_subdomains_julia.jl`

Typical command form:
- `JULIA_DEPOT_PATH=/tmp/ibamr-julia-depot julia /Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_ibamr_exported_subdomains_julia.jl <export_dir> <julia_out_dir> <matlab_out_dir> <tol>`

### 4) Core IBAMR parity executable
Executable:
- `/Users/boyceg/code/ibamr-objs-dbg/tests/IB/implicit_stokes_ib_operator_chain_parity_01`

Input used via orchestrator:
- `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.multilevel_l2.input`

## What We Are Checking and How

### A) Subdomain index-set parity
Check:
- Relaxed and strict CAV overlap sets exactly match reference construction.

How:
- Recompute reference sets in Octave/Julia from exported maps and row patterns.
- Compare sorted set equality subdomain-by-subdomain.

### B) Subdomain matrix parity
Check:
- Overlap block matrices match for both:
  - A00 velocity block diagnostics
  - full-saddle local block (`L(IS_k, IS_k)`) used by smoother semantics

How:
- Export canonical triplets per overlap block.
- Compare sparsity and values with tolerance.

### C) Subdomain processing order parity
Check:
- Same sweep metadata/order.

How:
- Compare sweep trace metadata `(sweep, subdomain, stage)` frame-by-frame.

### D) Subdomain RHS parity
Check:
- Same residual restriction vector on each subdomain at each frame.

How:
- Reconstruct `r = b - L*y_prev` in reference path.
- Compare restricted `r(IS_k)` against IBAMR-trace-derived state.

### E) One-level smoother action parity
Check:
- Same post-sweep state and trace.

How:
- Run one sweep with matching subdomains and pseudoinverse local solve semantics.
- Compare intermediate trace vectors and final smoother output vector.

### F) Transfer and V-cycle parity (two-level)
Check:
- Restriction (`R`), prolongation (`P`), and one V-cycle action parity.

How:
- Compare exported vectors and traces for:
  - `restriction_output_coarse`
  - `prolongation_output_fine`
  - V-cycle trace and output

## Current Status

## Working
1. Subdomain selection parity: PASS.
2. Subdomain block matrix parity: PASS.
3. Subdomain processing order parity: PASS.
4. Subdomain RHS parity: PASS (machine-precision agreement at first mismatch frame).
5. IBAMR export/import pipeline for structural checks and one-level action checks: PASS.
6. External Octave/Julia two-level bridge implementations now preserve the coarse-corrected iterate into post-smoothing (previous restart-from-zero bug removed).

## Oracle coverage
1. Octave bridge is currently the MATLAB/Octave-semantics baseline oracle for one-level smoother debugging:
- overlap-set parity
- overlap submatrix parity
- sweep metadata/order parity
- frame-by-frame one-level sweep trace parity
- final one-level smoother output parity

2. Julia bridge currently covers:
- overlap-set parity
- overlap submatrix parity
- sweep metadata/order parity
- frame-by-frame one-level sweep trace parity
- final one-level smoother output parity
- transfer/V-cycle vector and trace parity checks once two-level exports complete

3. Julia bridge is now a faster IBAMR-trace oracle for one-level drift localization, but Octave remains the baseline reference path when the question is specifically "does IBAMR match MATLAB/Octave semantics?".

## Partially working
1. One-level smoother parity:
- Not exact at `1e-12` combined abs/rel criterion.
- Mismatch is small numerical drift.
- After switching harness local solve from Eigen-SVD to direct LAPACK SVD, mismatch improved substantially.

Quantified improvement for `K=1e4` one-level:
- Before (Eigen-based local SVD in harness):
  - `max_abs ~= 6.054e-05`, `max_rel ~= 1.653e-09`
  - failing entries at tol=1e-12: 759
- After (LAPACK-based local SVD in harness):
  - `max_abs ~= 5.458e-07`, `max_rel ~= 1.642e-10`
  - failing entries at tol=1e-12: 72

## Broken / blocked
1. Two-level gate is unstable:
- `implicit_stokes_ib_operator_chain_parity_01` intermittently hangs in two-level mode for tested cases (`K=1e4`, `N=8/16`).
- This still blocks consistent transfer/V-cycle parity completion even though the external Octave/Julia post-smoothing bridge bug has now been fixed.

2. One-level strict parity threshold:
- `1e-12` criterion still fails for smoother trace/output due to tiny accumulated numerical differences.
- Practical agreement is much tighter now, but not bitwise/1e-12 exact.

## Convergence-Focused Interpretation

The main goal is not one-level bitwise parity; it is to explain why IBAMR needs more outer iterations, and sometimes diverges, at high stiffness compared to the MATLAB reference.

Current read:
1. One-level relaxed CAV sweep semantics are no longer the leading suspect:
- PETSc-delta replay reproduces the exported IBAMR one-level relaxed sweep trace exactly on the reproduced `K=1e4`, `N=8` case.
- This effectively validates subdomain order, local update application, and pressure projection for the one-level sweep path.

2. The strongest robustness clue now comes from the audited MATLAB reference itself:
- the MATLAB audit removed the `K=1e4` breakdown on the reference side only after switching the active multigrid path to deterministic pseudoinverse solves for both
  - local relaxed-CAV patch solves, and
  - the coarsest full-level solve;
- the audit also documented near-singular coarse-level relaxed patches, especially on levels coarser than the fine grid.
- with that audited reference contract, the `N=16` MATLAB/Octave smoke/example solve stays in single-digit GMRES counts across stiffness decades:
  - `K=1e2 -> 4`
  - `K=1e3 -> 5`
  - `K=1e4 -> 7`
  - `K=1e5 -> 8`
  - `K=1e6 -> 9`
- the important qualitative point is that MATLAB iteration growth with stiffness is mild and nearly linear, whereas the current IBAMR `ASM` path grows much faster and eventually diverges on the tested stiff `ex0` case.

3. That shifts the likely IBAMR robustness gap away from one-level sweep ordering and toward multilevel/coarse semantics and local-solve contract:
- coarse-grid solve/nullspace handling,
- transfer/V-cycle semantics,
- and whether the active IBAMR local patch solves are actually using MATLAB-style pseudoinverse semantics on stiff cases.

4. There is also a configuration caveat that can make the comparison apples-to-oranges:
- the shipped `examples/IB/implicit/ex0/input2d` default still uses `pc_type = "asm"` with coupling-aware subdomains, not the shell `coupling_aware_vanka_matlab` path;
- the alternate `petsc_options.coupling_aware_vanka_matlab.dat` file switches to the shell path, but it also explicitly sets `-stokes_ib_pc_level_sub_pc_type lu`;
- in `ibtk/src/solvers/impls/PETScLevelSolver.cpp`, the shell `coupling_aware_vanka_matlab` path defaults those local subsolves to `PCSVD` before `KSPSetFromOptions()`, so that options file can override the intended pseudoinverse-style default back to LU.

Working hypothesis:
- if a stiff IBAMR run is using default ASM, then worse iteration counts than MATLAB are expected because the preconditioner family is different;
- if a stiff IBAMR run is using the shell CAV path but its local subsolves are still being overridden to LU-like solves, then the remaining robustness gap is plausibly the same issue the MATLAB audit already exposed;
- if stiff IBAMR runs already use shell CAV with true SVD local solves, then the remaining gap is most likely in the multilevel/coarse path, not the validated one-level sweep itself.

## Convergence Experiment: ASM vs Shell CAV on Stiff `ex0`

On `2026-03-14`, a focused stiffness comparison was run against:
- optimized example executable: `/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/main2d`
- base input: `/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d.running`
- fixed outer solver contract: right-preconditioned `FGMRES`
- varying only the effective level-preconditioner family on the active stiff example

Persistent repro helper:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh`

Command used:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d.running`

Observed outer-iteration results:

| Mode | `K=1e4` | `K=1e6` | `K=1e7` |
| --- | --- | --- | --- |
| coupling-aware `ASM` | converged in `8` | converged in `30` | diverged at `200` (`DIVERGED_ITS`) |
| shell CAV (`coupling_aware_vanka_matlab`) | converged in `4` | converged in `10` | converged in `12` |

Implications:
1. The dominant convergence gap on this stiff example is already explained by preconditioner family:
- coupling-aware `ASM` degrades much faster with stiffness than the shell CAV path;
- by `K=1e7`, `ASM` diverges while the shell CAV path still converges rapidly.

2. This strongly supports the idea that "IBAMR is less robust than MATLAB at high stiffness" is not primarily a one-level sweep-order bug:
- the current shell CAV path is much closer to the robust MATLAB reference behavior than the default `ASM` path is.

3. A follow-up attempt to separate shell-local `LU` versus shell-local `SVD` via PETSc options files did **not** succeed:
- the shell subsolver options in those files were reported as unused on the multilevel `ex0` path;
- therefore the shell rows above should be interpreted as the built-in shell `coupling_aware_vanka_matlab` contract actually active in code, not as a validated LU-vs-SVD shell comparison.

4. Additional configuration caveat:
- the repository's PETSc options file `examples/IB/implicit/ex0/petsc_options.coupling_aware_vanka_matlab.dat` uses option names like `-stokes_ib_pc_level_shell_pc_type` and `-stokes_ib_pc_level_sub_pc_type`;
- on the tested multilevel `ex0` path, those shell-specific options were reported as unused, so that file should not currently be treated as a reliable oracle for shell-local solve policy.

## Convergence Experiment: Standard `ex0` `N=16` Sweep

To get closer to the audited MATLAB `N=16` smoke/example scale, the same helper was then run against:
- base input: `/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d`
- same optimized executable: `/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/main2d`
- same outer solver contract: right-preconditioned `FGMRES`
- with the helper forcing `ASM` versus shell `coupling_aware_vanka_matlab` directly in the input database instead of relying on shell-specific PETSc option names

Command used:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d`

Observed outer-iteration results:

| Mode | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| coupling-aware `ASM` | diverged at `200` | diverged at `200` | diverged at `200` | diverged at `200` | diverged at `200` |
| shell CAV (`coupling_aware_vanka_matlab`) | `6` | `8` | `10` | `12` | `16` |

Interpretation:
1. On this standard `ex0` hierarchy, the current `ASM` path is not a viable robustness match for the MATLAB reference at all.
2. The shell CAV path is dramatically better than `ASM`, so moving from `ASM` to MATLAB-style shell CAV is not an optimization detail; it is a qualitative solver-behavior change.
3. However, shell CAV still does **not** match the audited MATLAB robustness trend:
- audited MATLAB `N=16` smoke/example counts were `4, 5, 7, 8, 9` for `K=1e2, 1e3, 1e4, 1e5, 1e6`;
- IBAMR shell CAV on this `ex0` case gave `6, 8, 10, 12, 16`;
- so escaping `ASM` explains a large part of the gap, but not all of it.
4. This makes the next most likely culprit the multilevel/coarse path under shell CAV rather than the already-mostly-validated one-level relaxed sweep semantics.

## Outer Krylov Side Clarification

On `2026-03-14`, the outer-Krylov contract was rechecked against both the local PETSc source tree and the MATLAB/Octave driver behavior:
- IBAMR's implicit solver setup explicitly selects `KSPFGMRES` in `src/IB/IBImplicitStaggeredHierarchyIntegrator.cpp`;
- PETSc `3.23.3` documents that left preconditioning is the default for most Krylov methods **except** `KSPFGMRES`, which only supports right preconditioning;
- the PETSc `FGMRES` implementation advertises supported norms only for `PC_RIGHT`, so the current IBAMR `-ib_ksp_type fgmres` runs are already on the right-preconditioned path even without an explicit `-ib_ksp_pc_side right`;
- MATLAB/Octave `gmres(..., M)` is a left-preconditioned built-in driver, so it should continue to be treated as a smoke/example script rather than as the canonical outer iteration contract for IBAMR parity.

Conclusion:
- no IBAMR code change to left-preconditioned `GMRES` is needed, and making that change would actually move the comparison away from the audited Chapter 4 / shell-CAV contract;
- adding `-ib_ksp_pc_side right` to PETSc option files may still be worthwhile as an explicit documentation/guardrail choice, but it should not change solver semantics on the current `FGMRES` path.

## Convergence Experiment: `test.m`-Compatible Outer Driver

On `2026-03-14`, the `ex0` stiffness helper was extended with a third shell-CAV mode that approximates the MATLAB `test.m` outer driver:
- outer solver: `GMRES` instead of `FGMRES`
- preconditioning side: left
- norm type: preconditioned residual
- tolerance / cap: `rtol = 1e-8`, `max_it = 150`
- restart: `150`, so there is effectively no restart before the MATLAB-style iteration cap

Command used:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d`

Observed standard-`ex0` `N=16` outer-iteration results:

| Mode | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| shell CAV, current IBAMR-style outer driver | `6` | `8` | `10` | `12` | `16` |
| shell CAV, `test.m`-compatible outer driver | `6` | `8` | `10` | `11` | `14` |
| audited MATLAB smoke/example reference | `4` | `5` | `7` | `8` | `9` |

Interpretation:
1. Matching the MATLAB smoke driver's outer-GMRES contract does help slightly at the upper end:
- `K=1e5`: `12 -> 11`
- `K=1e6`: `16 -> 14`
2. However, the effect is modest relative to the full remaining gap to MATLAB:
- the outer-driver change does **not** close the `K=1e2..1e4` gap at all on this case;
- even after aligning to a `test.m`-style left-preconditioned `GMRES`, IBAMR remains substantially weaker than the audited MATLAB counts.
3. So the outer-driver mismatch is real enough to be worth recording, but it is not the main explanation for the robustness gap.

## Next Diagnostic Direction: Preconditioner-Only FAC Apply

On `2026-03-14`, the direct-V-cycle question was rechecked against the IBAMR FAC implementation:
- `IBTK::FACPreconditioner::solveSystem()` zeroes the initial guess and applies exactly one FAC cycle, so it is the closest in-process IBAMR analogue to MATLAB's `M(x) = v_cycle(L, x, zeros(...), ...)`;
- the current parity harness exports one-level sweep traces and a two-level V-cycle trace, but that two-level path is still a reconstructed cycle assembled from exported `L`, `R`, `P`, smoother overlap sets, and the audited pseudoinverse contract;
- therefore the current smooth-probe V-cycle diagnostics are already useful, but they do not yet prove that the *production* FAC scheduling path diverges at the same place.

Practical implication:
- to localize where IBAMR and MATLAB start to separate substantially, the next cleaner comparison target is a direct single application of `FACPreconditioner::solveSystem()` on a fixed RHS, not the full outer Krylov solve;
- ideally that RHS should include both a smooth probe and a real stiff-run Krylov residual snapshot, so we can distinguish "synthetic smooth-mode drift" from the residuals that actually drive the observed robustness gap.

## Convergence Check: `MAX_LEVELS=1` Collapses to the Coarsest Solve

To separate one-level patch action from multilevel hierarchy effects, the same helper was rerun against the same standard `ex0` input while forcing:
- `MAX_LEVELS = 1`

Command used:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh --max-levels 1 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d`

Observed outer-iteration results:

| Mode | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| coupling-aware `ASM` | `1` | `1` | `1` | `1` | `1` |
| shell CAV (`coupling_aware_vanka_matlab`) | `1` | `1` | `1` | `1` | `1` |

Interpretation:
1. This does **not** show that the one-level `ASM` or shell CAV smoother is exact in general.
2. With the current `ex0` input, forcing `MAX_LEVELS=1` makes the only remaining level the FAC coarsest level, and the FAC driver calls the coarsest solver path instead of the usual level smoother path.
3. On this input, the coarsest solver is configured as a whole-level `PETSC_LEVEL_SOLVER` with `ksp_type = "preonly"` and `pc_type = "svd"`, so the `MAX_LEVELS=1` run is effectively testing a coarsest pseudoinverse solve, not the ordinary one-level shell/ASM relaxation behavior.
4. The main remaining suspect is therefore still the hierarchy/V-cycle path itself:
- coarse solve semantics,
- coarse pressure-nullspace handling,
- interlevel transfer/correction,
- or some other multilevel control-path difference.
5. This also means that comparing the `MAX_LEVELS=1` `ex0` counts directly to the audited MATLAB `N=16` counts is not the right parity target:
- the important result is not that IBAMR matches MATLAB at one level,
- it is that the hierarchy configuration matters enough to completely change the outer-solver behavior.
6. Follow-up:
- forcing `coarse_solver_type = LEVEL_SMOOTHER` initially exposed a real coarsest-level null-dereference bug in `StaggeredStokesIBLevelRelaxationFACOperator`;
- after fixing that bug, the corrected `MAX_LEVELS=1` experiment still produced `1` outer iteration for both `ASM` and shell CAV across `K=1e2..1e6`;
- therefore this particular one-level `ex0` benchmark is not discriminating enough to explain the stiffness-robustness gap, even once the exact-coarse-solve artifact is removed.

## Convergence Experiment: Galerkin-like Stokes Toggle

To test whether the remaining shell-CAV gap was mainly caused by rediscretized Stokes operators on coarse levels, the standard multilevel `N=16` `ex0` sweep was rerun with:
- `rediscretize_stokes = FALSE`
- `res_rediscretized_stokes = FALSE`

Command used:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh --rediscretize-stokes FALSE --res-rediscretized-stokes FALSE /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d`

Observed outer-iteration results:

| Mode | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| coupling-aware `ASM` | diverged at `200` | diverged at `200` | diverged at `200` | diverged at `200` | diverged at `200` |
| shell CAV (`coupling_aware_vanka_matlab`) | `6` | `8` | `10` | `12` | `17` |

Interpretation:
1. Turning both rediscretization flags off does **not** materially close the shell-CAV gap to the audited MATLAB reference.
2. On this case it is slightly worse than the default shell-CAV trend at the stiffest point (`17` instead of `16` at `K=1e6`).
3. So the remaining robustness gap is not well explained by a simple `rediscretize_stokes` versus Galerkin coarse-operator switch by itself.

## Convergence Experiment: RT0 Velocity Transfer on FAC Corrections

The standard multilevel `N=16` `ex0` sweep was then rerun with FAC correction-transfer methods chosen to better match the MATLAB reference velocity transfer semantics:
- `U_prolongation_method = RT0_REFINE`
- `U_restriction_method = RT0_COARSEN`

Command used:
- `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh --u-prolongation-method RT0_REFINE --u-restriction-method RT0_COARSEN /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0 /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d`

Observed outer-iteration results:

| Mode | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| coupling-aware `ASM` | diverged at `200` | diverged at `200` | diverged at `200` | diverged at `200` | diverged at `200` |
| shell CAV (`coupling_aware_vanka_matlab`) | `6` | `8` | `9` | `11` | `15` |

Interpretation:
1. Aligning the FAC correction velocity transfers with RT0 semantics gives a modest but consistent shell-CAV improvement at moderate and high stiffness.
2. Relative to the default shell-CAV sweep (`6, 8, 10, 12, 16`), this reduces the counts to `6, 8, 9, 11, 15`.
3. This is still not enough to match the audited MATLAB trend (`4, 5, 7, 8, 9`), but it is the first multilevel configuration change in IBAMR on this case that clearly improves the shell-CAV stiffness trend without changing the one-level local solve policy.
4. Current read: transfer semantics are a real contributor to the remaining robustness gap, even if they are not the whole story.

## Diagnostic Experiment: Smooth Two-Level V-cycle Probe

To make the multilevel comparison less sensitive to arbitrary probe content, the two-level bridge workflow was rerun on the reproduced `K=1e4`, `N=8` case with:
- `vcycle_probe_type = "smooth_eulerian_mode"`
- `vcycle_probe_mode_x = 1`
- `vcycle_probe_mode_y = 1`

Command used:
- `VCYCLE_PROBE_TYPE=smooth_eulerian_mode VCYCLE_PROBE_MODE_X=1 VCYCLE_PROBE_MODE_Y=1 /Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/run_ibamr_cav_subdomain_bridge_parity.sh /Users/boyceg/code/ibamr-objs-dbg /tmp/ibamr-smooth-vcycle-probe /Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.multilevel_l2.input 1.0e4 two_level 8`

Observed bridge outcome:
1. The smooth probe export completed successfully and now writes `vcycle_input_state.level1.txt` in addition to the usual `vcycle_input_rhs.level1.txt`.
2. The Julia bridge was extended to emit per-stage norm summaries for the two-level trace:
- `julia_vcycle_norm_report_relaxed.level1.txt`
- `julia_vcycle_norm_report_strict.level1.txt`
3. On the relaxed smooth probe, the first two-level trace mismatch occurs already at:
- `stage = y_pre`
- not first at `coarse_rhs`, `coarse_sol`, `fine_correction`, or `y_post`
4. The relaxed norm report shows no dramatic new jump isolated to the coarse solve:
- relative stage-norm differences stay roughly `1.4e-12` to `1.5e-12` from `y_pre` through `y_after_correction`
- `y_post` remains the same order, about `1.1e-12`
5. The strict smooth-probe norm report shows the same qualitative pattern:
- relative stage-norm differences are about `2.6e-12` through `y_after_correction`
- they grow only mildly by `y_post`, to about `4.9e-12`

Interpretation:
1. This smooth-probe experiment does **not** support the idea that the coarse solve is the first point where the two-level relaxed V-cycle starts to separate.
2. On this case, the multilevel mismatch is already present in the fine-level pre-smoothing output `y_pre`, and the later residual/coarse/correction stages mostly propagate that difference instead of introducing a qualitatively new norm jump.
3. That does not rule out multilevel transfer/coarse semantics as contributors to outer-iteration robustness, but it does mean the smooth-probe norm diagnostic points back to the fine-level smoother semantics as the earliest visible source of V-cycle divergence on this reproduced export.

## Key Implementation Changes In This Pass

In solver/fac code:
- `/Users/boyceg/code/IBAMR/src/IB/StaggeredStokesIBLevelRelaxationFACOperator.cpp`

Changed:
- fixed the `coarse_solver_type = "LEVEL_SMOOTHER"` path so the coarsest level now initializes a regular level solver instead of leaving `d_level_solvers[d_coarsest_ln]` null and then crashing when `solveCoarsestLevel()` routes through `smoothError()`.

Intent:
- make corrected one-level/coarsest-smoother experiments possible without PETSc segfaults.

Outcome:
- removes a real null-dereference bug from the Stokes+IB FAC operator;
- confirms the corrected `MAX_LEVELS=1` experiment is still non-discriminating on this `ex0` case even after the bug is fixed.

In parity/export tooling:
- `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp`
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_ibamr_exported_subdomains_julia.jl`
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/run_ibamr_cav_subdomain_bridge_parity.sh`

Changed:
- added `vcycle_probe_type = "smooth_eulerian_mode"` plus `vcycle_probe_mode_x`, `vcycle_probe_mode_y` to generate a low-frequency V-cycle probe state from the Eulerian DOF map;
- export now writes `vcycle_input_state.level*.txt` alongside the RHS for V-cycle bridge cases;
- Julia bridge now writes per-stage two-level norm reports for relaxed and strict V-cycles;
- Julia bridge no longer throws away higher-level/two-level diagnostics just because a lower-level one-level mismatch was already detected;
- the bridge runner can now inject the smooth V-cycle probe through environment variables.

Intent:
- support step-by-step multilevel diagnosis on a common smooth probe instead of only on the original index-based synthetic vector.

Outcome:
- makes it possible to see that, on the reproduced `K=1e4`, `N=8` smooth-probe case, the relaxed two-level trace first diverges at `y_pre`, not first at the coarse solve stage.

In test harness:
- `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp`

Changed:
- `solve_pseudoinverse()` switched from Eigen `JacobiSVD` to direct LAPACK `gesvd` call via PETSc BLAS/LAPACK headers.

Intent:
- Reduce solver-backend mismatch between IBAMR harness and reference pipelines.

Outcome:
- Significant reduction in smoother parity drift, but not full 1e-12 equality yet.

In external bridge comparators:
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_ibamr_exported_subdomains_octave.m`
- `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_ibamr_exported_subdomains_julia.jl`

Changed:
- Two-level V-cycle reference path now carries `y_after_correction` into the post-smoothing sweep instead of restarting that sweep from zero.
- Julia bridge now compares one-level sweep traces frame-by-frame against exported IBAMR traces and reports the first mismatching frame/dof.
- Julia bridge `solve_pseudoinverse()` now uses explicit LAPACK SVD with the same `max(m,n) * eps * sigma_max` cutoff rule used by the C++ parity harness.
- The C++ parity harness now supports targeted local-solve dumps of:
  - `rhs_full`
  - `residual_full`
  - selected full-row `MatGetRow` diagnostics with exact `Float64` bit payloads for the row entries and state values
- The C++ parity harness now also exports one-level PETSc-backed sidecar traces for:
  - full residual-before-local-solve (`smoother_residual_trace_*`)
  - applied full-space local correction (`smoother_delta_trace_*`)
- The Julia bridge now supports two extra one-level replay modes when those sidecars are present:
  - PETSc-residual replay: use exported harness residuals, but still solve the local pseudoinverse in Julia
  - PETSc-delta replay: use exported harness local corrections directly and replay only the update/projection path
- Float-valued export writers in the harness now emit scientific `max_digits10` formatting instead of the previous default-format `setprecision(17)`.

Intent:
- Make external Octave/Julia two-level oracle semantics match the IBAMR harness V-cycle structure.

Outcome:
- Removes a known oracle bug from the multilevel bridge path and upgrades Julia from a final-action-only checker to a trace-aware one-level oracle.
- On the reproduced `K=1e4`, `N=8` export path, the Julia relaxed first-mismatch frame moved from `29` to `43` after aligning the local pseudoinverse backend with the C++ harness, much closer to the Octave relaxed first mismatch at frame `45`.
- Direct frame-43 local-solve comparison now shows:
  - exported local block matches the harness-local dumped block exactly,
  - Julia and the harness-local C++ pseudoinverse agree to about `1e-11` on the exact same dumped local block/RHS,
  - the harness-local delta matches the IBAMR trace delta to about `1e-12`,
  - the remaining larger frame-43 discrepancy comes from the incoming local RHS/state drift, not from the pseudoinverse backend itself.
- Frame-growth analysis on the reproduced relaxed level-0 sweep now shows:
  - the first parity-threshold crossing occurs at stage `0` after subdomain `10` (not at the pressure-projection stage),
  - by frame `10`, the harness and exported `y_before` states still agree to about `7e-13`,
  - at that same frame, the dominant difference is already in the full residual/matvec path before local restriction: same-run comparisons show about `2.9e-11` full-residual mismatch from the export state and about `5.5e-12` even when recomputing `b - L*y` from the exact harness `y_before`,
  - the locally restricted RHS at frame `10` then differs by about `2e-11` from the export path and about `2e-12` from the exact harness-state recomputation,
  - by frame `36`, both effects matter: the incoming full state has drifted to about `1.6e-11`, and the local RHS reconstructed from the harness state still differs from the harness local RHS by about `2.3e-10`.
- Repeating these comparisons against same-run exports (not the original reproduced export directory) gives the same numbers, so the residual mismatch is not an artifact of cross-run variation.
- Targeted row-level dumps at frame `10` / row `120` and frame `36` / row `41` now show:
  - on the exact in-process harness data, PETSc `MatMult` residual and an explicit `MatGetRow` row-dot reconstruction agree to machine precision,
  - the frame-36 row-41 mismatch persists even when replaying the dumped row using the exact serialized `Float64` bit patterns for the coefficients and `y_before` entries,
  - switching the export writers to scientific `max_digits10` formatting does not remove that row-41 replay gap,
  - therefore the remaining same-state residual mismatch is no longer evidence of a wrong exported operator/RHS/state; it is consistent with low-level in-process accumulation or excess-precision behavior that the offline Julia/Octave replay does not reproduce bit-for-bit.
- On a fresh reproduced `K=1e4`, `N=8`, two-level export with the new sidecar traces:
  - plain Julia semantic replay still first mismatches at relaxed frame `43`,
  - PETSc-residual replay pushes the first mismatch later, to relaxed frame `59`, with much smaller error,
  - PETSc-delta replay matches the exported IBAMR relaxed one-level trace exactly.
- Targeted exact local-solve inspection at the PETSc-residual first-mismatch frame (`59`) now shows:
  - the dumped harness local block matches the export-side block exactly,
  - the dumped harness local RHS matches the exported PETSc residual sidecar exactly,
  - the dumped harness local delta matches the exported PETSc delta sidecar exactly,
  - Julia on that exact frame-59 harness block/RHS differs from the harness local delta by only about `4.5e-13`,
  - in PETSc-residual replay, the frame-59 `y_before` state is already off by about `7.3e-11`, while the frame-59 local delta differs by only about `1.8e-12`.
- Interpretation:
  - the sweep ordering, local update application, and pressure-projection path are now effectively confirmed,
  - the residual path is the dominant source of early offline replay drift,
  - after replacing the residual path with the exported PETSc residuals, the remaining drift is confined to the residual-to-delta local solve/backend step,
  - by the PETSc-residual first-mismatch frame, most of that remaining drift is accumulated backend noise from earlier local solves rather than a large new per-frame error at frame `59`,
  - once the exported PETSc local deltas themselves are replayed, the Julia bridge can reproduce the IBAMR one-level trace exactly on this reproduced case.
- These changes do not by themselves resolve the separate intermittent two-level hang in the IBAMR parity executable.

## New In-Process FAC/Smoother Trace Hooks

Intent:
- Stop inferring the relaxed CAV sweep only from exported/reconstructed bridge data and instead trace the real in-process IBAMR smoother path during a FAC application.

Implementation:
- Added opt-in FAC stage tracing in `IBTK::FACPreconditioner` via:
  - `IBAMR_FAC_TRACE_FILE`
  - `IBAMR_FAC_TRACE_TARGET_CALL`
- Added opt-in per-subdomain shell-CAV tracing in `IBTK::PETScLevelSolver::PCApply_CouplingAwareVankaMatlab()` via:
  - `IBAMR_PETSC_LEVEL_TRACE_FILE`
  - `IBAMR_PETSC_LEVEL_TRACE_TARGET_CALL`
  - `IBAMR_PETSC_LEVEL_TRACE_TARGET_LEVEL`
- The level-solver trace writes one row per stage:
  - `start`
  - `restricted_rhs`
  - `local_delta`
  - `post_prolong`
  - `postprocess`
- Each row records level/solve-call/subdomain plus `L2` and `L_inf` norms of:
  - the shell RHS,
  - the accumulated global correction,
  - the full shell residual,
  - the restricted local RHS,
  - the local subdomain correction.

Smoke check:
- Rebuilt `IB-implicit-ex0` successfully with the new trace hooks.
- Ran `/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/main2d input2d.running` with:
  - `IBAMR_PETSC_LEVEL_TRACE_FILE=/tmp/ibamr_petsc_level_trace.txt`
  - `IBAMR_PETSC_LEVEL_TRACE_TARGET_CALL=1`
  - `IBAMR_PETSC_LEVEL_TRACE_TARGET_LEVEL=1`
- The run converged in `10` outer iterations and produced `/tmp/ibamr_petsc_level_trace.txt`.
- That trace confirms we can now observe the true in-process relaxed sweep step-by-step at the level-solver seam, instead of only reconstructing it offline.

Interpretation:
- There is not a distinct `StaggeredStokesIBFACLevelSmootherOps` class in this tree, but the actual hook point the smoother uses is the combination of:
  - `StaggeredStokesIBLevelRelaxationFACOperator::smoothError()` for level residual/error handoff, and
  - `PETScLevelSolver::PCApply_CouplingAwareVankaMatlab()` for the relaxed coupling-aware Vanka sweep itself.
- So this new PETSc-level trace is the right place to compare a common smooth vector through both IBAMR and MATLAB-style V-cycle logic before outer Krylov effects mix things together.

## Focus Shift: Julia Reference Package

Intent:
- Stop treating `compare_ibamr_exported_subdomains_julia.jl` as a second independent Julia reference implementation and instead source canonical relaxed CAV semantics from `/Users/boyceg/code/implicit_ib_coupling_aware_vanka/julia-reference`.

Outcome:
- The exported-subdomain Julia bridge now loads the audited Julia reference package directly from:
  - `/Users/boyceg/code/implicit_ib_coupling_aware_vanka/julia-reference`
  - overridable with `IBAMR_JULIA_REFERENCE_ROOT`
- The canonical relaxed overlap construction in the bridge now comes from `ImplicitIBCAVReference.extract_coupled_dofs(...)` instead of the older local `relaxed_set(...)` copy.
- The canonical local/coarse pseudoinverse in the bridge now comes from `ImplicitIBCAVReference.solve_reference_pseudoinverse(...)` instead of the older local `solve_pseudoinverse(...)` copy.
- The bridge still retains IBAMR-specific diagnostics that are not part of the canonical package:
  - strict-overlap reconstruction,
  - exported algebraic `R/P` replay,
  - PETSc-residual replay,
  - PETSc-delta replay.

Validation:
- Parse/load check:
  - `julia --startup-file=no -e 'include(".../compare_ibamr_exported_subdomains_julia.jl")'`
- Julia reference smoke solve check:
  - direct package call to `solve_smoke_example(N=8, depth=3, K=1e2)` returned `iterations=2`, `converged=true`
- Existing-export bridge check:
  - running the updated bridge on `/tmp/ibamr-smooth-vcycle-probe/two_level/ibamr-export` now fails immediately on relaxed-overlap structure, before the older floating-point drift diagnostics:
    - level `0`, subdomain `10`, seed `20`: IBAMR exported relaxed patch size `83`, canonical Julia-reference patch size `67`
    - level `1`, subdomain `37`, seed `74`: IBAMR exported relaxed patch size `127`, canonical Julia-reference patch size `118`
  - in both spot-checked mismatches, the canonical Julia-reference patch is a strict subset of the exported IBAMR patch; there were no `only_ref` DOFs.
- Current implication:
  - the relaxed path used by the Julia bridge is now tied to the audited canonical Julia reference package rather than to a duplicate copy embedded in the bridge script,
  - strict-path and PETSc-replay logic remain useful diagnostics, but they should no longer be confused with the canonical reference definition.
  - more importantly, the first current IBAMR-vs-canonical difference may now be structural patch construction, not only later residual/local-solve drift.

## Problem Setup Equivalence Check

Intent:
- Double check that the IBAMR cases under discussion are actually using the same geometry, IB point count, and stiffness convention as the canonical Julia reference package.

Canonical Julia reference contract:
- `build_reference_geometry()` in `/Users/boyceg/code/implicit_ib_coupling_aware_vanka/julia-reference/src/problem.jl` uses:
  - `dx = 1 / N`
  - `ds0 = dx / r_cyl` with `r_cyl = 0.25`
  - `approx = round(2*pi / ds0)`
  - `ds = 2*pi / approx`
  - ellipse parameters `alpha = 0.23`, `beta = r_cyl^2 / alpha`
- `apply_elastic_laplacian()` in both the Julia and MATLAB references applies the elastic operator as `K / ds^2` times the periodic second difference.

Concrete setup comparison:

| Configuration | Finest Eulerian grid | IB geometry | IB points | `ds` | Input stiffness meaning |
| --- | --- | --- | --- | --- | --- |
| Julia/MATLAB reference, `N=8` | `8 x 8` | ellipse `alpha=0.23`, `beta=0.25^2/0.23` | `13` | `0.483321946706122` | `K` is the reference elastic coefficient; the operator uses `K / ds^2` |
| Julia/MATLAB reference, `N=16` | `16 x 16` | same ellipse | `25` | `0.251327412287183` | same |
| Julia/MATLAB reference, `N=32` | `32 x 32` | same ellipse | `50` | `0.125663706143592` | same |
| `implicit_stokes_ib_operator_chain_parity_01` | fixed `32 x 32` in code | same ellipse | `50` | `0.125663706143592` | input `SPRING_STIFFNESS = K`; the test then applies `spring_k = K / ds^2` internally |
| `examples/IB/implicit/ex0/input2d` | finest `32 x 32` for the shipped `N=16`, `MAX_LEVELS=2` case | circle `r = 0.2` | `40` | input `DS = 0.03125` | input `SPRING_STIFFNESS = K / ds`; for `K = 1e4`, this is `3.2e5` |

Update on `2026-03-15`:
- `implicit_stokes_ib_operator_chain_parity_01.cpp` was updated so the parity test no longer hard-wires the finest-grid contract.
- It now derives:
  - finest Eulerian `N` from the input `N`, `REF_RATIO`, and `MAX_LEVELS`,
  - `dx = L / N_finest`,
  - `dt = 0.5 * dx`,
  - `n_ib = round(2*pi / (dx / r_cyl))`
  from the requested problem setup before the IB geometry and linearized operator are built.
- The bridge export path now also writes `problem_setup.txt` so the effective parity-test setup can be checked directly from the generated artifacts.
- Direct validation:
  - a rebuilt one-level `N=16` export wrote `problem_setup.txt` with `finest_grid_cells = 16`, `dx = 0.0625`, `dt = 0.03125`, `num_lag_nodes = 25`, `lag_ds = 0.251327412287183`;
  - the existing intermittent two-level hang still prevented a matching direct artifact check on the base-`N=16`, `MAX_LEVELS=2` case during this pass, but the new setup logic should make that case use `finest_grid_cells = 32`, `dt = 0.015625`, and `num_lag_nodes = 50`.
- Follow-up comparison against the canonical Julia bridge on the corrected one-level `N=16` export:
  - the export-side setup now truly matches the Julia `N=16` geometry/stiffness contract;
  - a targeted debug hook in `StaggeredStokesPETScMatUtilities.cpp` showed that the IBAMR relaxed construction was building its initial velocity set from the PETSc row pattern of `A00`, not from the numerically nonzero entries of the row:
    - the debug `initial_velocity_dofs` matched `A00.level0.row_pattern.txt`,
    - but it was strictly larger than the canonical Julia seed-row set extracted from the actual matrix values.
  - root cause:
    - `build_initial_velocity_set_from_seed_components()` was inserting every column returned by `MatGetRow()` without checking whether the corresponding stored value was numerically zero.
    - this widened the relaxed patch through structural zeros before the pressure/closure expansion.
  - fix:
    - `build_initial_velocity_set_from_seed_components()` now filters out zero-valued entries from `MatGetRow()`.
  - after the fix, the same targeted one-level `N=16` debug case gave:
    - `initial_velocity_dofs` exactly matching the Julia seed-row velocity set,
    - `involved_cell_dofs` exactly matching the Julia pressure set,
    - `closure_dofs` / `overlap_dofs` exactly matching the canonical Julia relaxed patch.
  - comparator result on the corrected export:
    - the relaxed mismatch disappeared,
    - the Julia bridge now fails only on the non-canonical `STRICT` path (`strict IS mismatch at level 0 subdomain 37 seed 74`).

Current read:
- The dedicated parity test now matches the canonical Julia reference setup for both:
  - one-level `N=16` comparisons, which should use `25` IB vertices and `dt = 0.03125`;
  - two-level base-`N=16` comparisons, whose finest grid is `32 x 32`, matching the Julia `N=32` setup with `50` IB vertices and `dt = 0.015625`.
- `ex0` is still a different benchmark problem:
  - different geometry (`r = 0.2` circle instead of the audited ellipse),
  - different IB point count (`40` instead of `50` on the shipped case),
  - different input stiffness units (`K / ds` at input rather than `K` passed through an internal `1 / ds^2` scaling).

Implication:
- When comparing against the canonical Julia reference package, we should treat:
  - the parity test as the correct setup-match vehicle for the `N=32` fine-grid contract;
  - `ex0` as a useful robustness benchmark, but not a like-for-like Julia reference setup.
- The one-level bridge gate is now a valid strict Julia-problem parity case once the rebuilt executable is used.
- The relaxed overlap-construction mismatch on the one-level `N=16` case is now explained by a real IBAMR bug in the `A00` seed-row expansion: structural zeros were being treated as active couplings.
- Since this utility is also used by the actual PETSc level solver, it is now a plausible direct contributor to the shell-CAV robustness gap, not just a parity-harness artifact.

## Structural-Zero Relaxed-CAV Bug

Targeted diagnosis:
- A new targeted debug hook was added to `StaggeredStokesPETScMatUtilities.cpp` so a single relaxed or strict coupling-aware subdomain can dump:
  - `seed_velocity_dofs`
  - `initial_seed_components`
  - `initial_velocity_dofs`
  - `involved_cell_dofs`
  - `closure_dofs`
  - `overlap_dofs`
- On the corrected one-level `N=16` setup-matched parity case, the first targeted dump for relaxed level `0`, subdomain `37`, seed `74` showed:
  - `initial_velocity_dofs` matched the exported `A00.level0.row_pattern.txt` row, not the actual numerically nonzero `A00` row,
  - the canonical Julia reference `extract_coupled_dofs()` uses the actual numeric row entries of the velocity block,
  - therefore IBAMR was treating structurally present zero entries as active velocity couplings when constructing relaxed CAV patches.

Root cause:
- `build_initial_velocity_set_from_seed_components()` called `MatGetRow()` and inserted every returned column index into `initial_velocity_dofs` without checking whether the corresponding stored value was numerically zero.

Fix:
- `build_initial_velocity_set_from_seed_components()` now requests both columns and values from `MatGetRow()` and skips zero-valued entries.

Post-fix one-level parity result:
- The same one-level `N=16` canonical Julia comparison now no longer reports any relaxed mismatch.
- The Julia bridge on the corrected export now fails only on the non-canonical strict path:
  - `strict IS mismatch at level 0 subdomain 37 seed 74`

Interpretation:
- For the canonical relaxed CAV path, the one-level structural mismatch is now explained and fixed.
- Since this code path is shared with the real PETSc level solver, this is a credible solver-quality fix, not merely a bridge-export cleanup.

## Post-Fix Stiffness Sweep

Validation run:
- Rebuilt `IB-implicit-ex0` in `/Users/boyceg/code/ibamr-objs-opt`.
- Reran `/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_convergence_compare.sh /Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0`

New iteration table:

| Mode | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| `asm` | `4` | `5` | `8` | `12` | `87` |
| `shell_matlab` | `2` | `4` | `5` | `5` | `9` |
| `shell_testm` | `2` | `4` | `5` | `5` | `5` |

Comparison to the pre-fix shell-CAV baseline on this `ex0` case:
- previous `shell_matlab` counts were `6, 8, 10, 12, 16`
- post-fix `shell_matlab` counts are `2, 4, 5, 5, 9`

Current implication:
- The structural-zero relaxed-subdomain bug was not a cosmetic parity issue; it had a large effect on shell-CAV convergence robustness.
- After the fix, the shell `coupling_aware_vanka_matlab` path on this benchmark is now much closer to the audited MATLAB reference trend, including matching the old `K=1e6 -> 9` reference endpoint on the `N=16` table.

## Direct Count Comparison Against Current Julia Reference

Current Julia-package sweep:
- Ran the current canonical Julia package entry point directly:
  - `solve_smoke_example(N=16, depth=4, K=...)`
- Observed iteration counts:
  - `K=1e2 -> 4`
  - `K=1e3 -> 5`
  - `K=1e4 -> 7`
  - `K=1e5 -> 8`
  - `K=1e6 -> 7`

Side-by-side with the post-fix IBAMR `ex0` shell path:

| Case | `K=1e2` | `K=1e3` | `K=1e4` | `K=1e5` | `K=1e6` |
| --- | --- | --- | --- | --- | --- |
| Julia reference `solve_smoke_example(N=16, depth=4)` | `4` | `5` | `7` | `8` | `7` |
| IBAMR `ex0` `shell_matlab` | `2` | `4` | `5` | `5` | `9` |
| IBAMR `ex0` `shell_testm` | `2` | `4` | `5` | `5` | `5` |

Interpretation:
- The structural-zero fix clearly removed the large shell-CAV robustness deficit that was present before.
- The remaining difference is now much smaller and mixed in sign:
  - IBAMR `ex0` is still slightly worse than the current Julia package at `K=1e6`,
  - but it is now better on the milder stiffness points of this particular benchmark.
- Since `ex0` is still not the same geometry/unit contract as the Julia smoke problem, the remaining differences in this table should not be over-interpreted as solver bugs by themselves.
- There is also a new documentation follow-up:
  - the current Julia package gives `K=1e6 -> 7` at `N=16, depth=4`, which differs from the older audited MATLAB table entry `9`.
  - That means the package-backed Julia reference should now be treated as the live reference for future comparisons, and the older MATLAB audit table may need an update or an explicit “historical branch result” label.

## Matched Outer-GMRES Smoke Driver

New capability:
- Added an opt-in Julia-style outer-solve mode to `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp`.
- This mode:
  - uses the aligned ellipse setup already used by the parity executable,
  - builds the full explicit IBAMR multilevel eliminated Eulerian operator hierarchy,
  - forms the smoke RHS from the same linearized elastic force construction (`S * (A * X)`),
  - runs left-preconditioned `GMRES` with `KSP_NORM_PRECONDITIONED`,
  - uses the parity test’s reconstructed multilevel CAV V-cycle as a shell preconditioner,
  - and writes a compact summary file.

Implementation details:
- Added recursive multilevel V-cycle application and the shell-PC / summary plumbing in:
  - `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp`
- Key new entry points are around:
  - `MatrixVcycleShellContext`
  - `apply_multilevel_vcycle(...)`
  - `run_matrix_gmres_smoke(...)`

Follow-up correction:
- The first “matched depth-2” smoke comparison was not actually matched. The IBAMR input still had base `N=16` with `MAX_LEVELS=2`, so its finest grid was `32`, while Julia’s `solve_smoke_example(N=16, depth=2)` uses finest `16`.
- After fixing the hierarchy contract to use base `N=8` so the finest grid is truly `16`, the parity executable exported:
  - `dx = 0.0625`
  - `dt = 0.03125`
  - `num_lag_nodes = 25`
  - `lag_ds = 2*pi/25`
- That matches the canonical Julia package geometry/stiffness contract for the finest-`16` ellipse problem.

Root-cause fix for the smoke mismatch:
- The parity executable intentionally uses Julia-style springs with pair stiffness `K / ds^2`, but the smoke matrix path was still forming the Eulerian force terms with the unscaled transpose interpolation operator.
- Concretely, it was using:
  - `SAJ = J' * A * J`
  - `rhs = J' * (A * X)`
- for a Julia-normalized force Jacobian `A` that still needs the discrete spread factor.
- The corrected contract is:
  - `SAJ = (-dt * ds / h^2) * (J' * A * J)` in 2D
  - `rhs = (ds / h^2) * (J' * (A * X))`
- After patching `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp` to apply that scale:
  - the matched `N=16`, depth-2 smoke RHS agrees with the Julia reference to about `6.8e-12`,
  - and the big operator mismatch collapses to the expected pressure-row sign convention (`-div u = 0` in IBAMR vs `+div u = 0` in the Julia reference).

Corrected matched `N=16`, depth-2 comparison against the canonical Julia package:
- Julia reference (`solve_smoke_example(N=16, depth=2, K=...)`):
  - `K=1e2 -> 4`
  - `K=1e4 -> 6`
  - `K=1e6 -> 7`
- IBAMR aligned parity smoke driver after the spread-scale fix (`matrix_krylov_policy = relaxed`):
  - `K=1e2 -> 5`
  - `K=1e4 -> 6`
  - `K=1e6 -> 8`

Updated interpretation:
- The earlier “IBAMR matched smoke driver is completely wrong” result was mostly a parity-harness contract bug, not evidence that the reconstructed multilevel CAV V-cycle itself was fundamentally off.
- After fixing the spread scale, the true matched smoke driver is now within one iteration of Julia on the tested cases:
  - exact match at `K=1e4`,
  - one extra iteration at `K=1e2`,
  - one extra iteration at `K=1e6`.
- The remaining matrix difference on the matched problem is now dominated by the known pressure-block sign convention (`-div u = 0` in IBAMR to keep the saddle form symmetric).
- This substantially lowers the priority of further smoke-driver contract debugging and shifts attention back to the actual robustness problem in the real solver/preconditioner path.

Live shell/FAC contract check:
- The real FAC `SAJ` matrix was then checked directly in `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_solver_components_01.cpp` on a Julia-aligned finest-`16` ellipse case with:
  - `X_RADIUS = 0.23`
  - `Y_RADIUS = 0.27173913043478265`
  - `25` IB points
  - pair spring stiffness `K / ds_theta^2 = 1.5831434944115276e3`
- Result:
  - `fd_relative_error = 1.88727e-06`
  - `saj_cell_scaled_relative_error = 0`
  - `saj_theta_scaled_relative_error = 0.748673`
- Interpretation:
  - the live FAC `SAJ` path matches the native IBAMR nodal-force scaling
    - `SAJ = (-dt / h^2) * J' A J`
  - and does **not** match the Julia-density scaling used by the corrected explicit smoke harness
    - `SAJ = (-dt * ds_theta / h^2) * J' A J`
- So there are currently two equivalent-but-different normalization contracts in play:
  - explicit parity/matrix smoke:
    - Julia-density normalization
    - pair stiffness `K / ds_theta^2`
    - spread scale `ds_theta / h^2`
  - live shell/FAC path:
    - native IBAMR nodal-force normalization
    - pair stiffness `K / ds_theta`
    - spread scale `1 / h^2`
- This explains why the corrected explicit smoke driver can now match Julia closely while the real shell/FAC path still needs an explicit normalization choice to run the *same* problem.

Follow-up implementation:
- `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp` now supports:
  - `SPRING_NORMALIZATION = JULIA_DENSITY`
  - `SPRING_NORMALIZATION = IBAMR_NODAL`
- The export metadata now records:
  - `spring_normalization`
  - the actual `spring_pair_k`
  - the resulting `matrix_spread_scale`
- This should make future parity and live-shell comparisons much less ambiguous.

Live FAC smooth-probe mode:
- Added a new `run_live_fac_probe` mode to `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp`.
- This mode does **not** try to reuse the condensed matrix-smoke RHS, since that RHS lives in the finest-level condensed DOF space while the real FAC preconditioner acts on a composite hierarchy vector.
- Instead it applies the actual `StaggeredStokesIBJacobianFACPreconditioner::solveSystem()` to a common smooth hierarchy-space probe built level-by-level from the same analytic mode used by the reconstructed V-cycle probe:
  - `vcycle_probe_type = "smooth_eulerian_mode"`
  - `vcycle_probe_mode_x = 1`
  - `vcycle_probe_mode_y = 1`
- The mode can optionally export the hierarchy-space input/output vectors:
  - `live_fac_probe_rhs.level*.txt`
  - `live_fac_probe_sol.level*.txt`
- It also works with the existing in-process trace hooks:
  - `IBAMR_FAC_TRACE_FILE`
  - `IBAMR_FAC_TRACE_TARGET_CALL`
  - `IBAMR_PETSC_LEVEL_TRACE_FILE`
  - `IBAMR_PETSC_LEVEL_TRACE_TARGET_CALL`
  - `IBAMR_PETSC_LEVEL_TRACE_TARGET_LEVEL`

First matched live probe run:
- Ran the optimized parity executable on the matched finest-`16`, depth-2 ellipse setup with:
  - `SPRING_NORMALIZATION = IBAMR_NODAL`
  - shell level solver contract (`pc_type = shell`, `shell_pc_type = coupling_aware_vanka_matlab`)
  - smooth probe export directory: `/tmp/ibamr-live-fac-probe-true-n16-export`
  - FAC trace: `/tmp/ibamr-live-fac-probe-true-n16-fac-trace.txt`
  - level trace: `/tmp/ibamr-live-fac-probe-true-n16-level-trace.txt`
- Summary:
  - `live_fac_probe_success = 1`
  - `live_fac_probe_rhs_l2_norm.level0 = 5.8275209137333865`
  - `live_fac_probe_sol_l2_norm.level0 = 0.62157099509375002`
  - `live_fac_probe_rhs_l2_norm.level1 = 11.655041827466771`
  - `live_fac_probe_sol_l2_norm.level1 = 2.4770650110670793`
- The top-level FAC trace looked suggestive at first:
  - `rhs_input(level 1) l2 = 7.28440114216673762e-01`
  - `y_pre(level 1) l2 = 9.73618221974046211e-02`
  - `residual_fine(level 1) l2 = 3.69476404172269457e-01`
  - `coarse_rhs(level 0) l2 = 0`
  - `coarse_sol(level 0) l2 = 0`
  - `y_post(level 1) l2 = 1.54816563191692458e-01`
- Follow-up diagnosis:
  - this is a diagnostics artifact, not evidence that the coarse correction is actually zero;
  - `FACPreconditioner::recordDebugStage()` builds a one-level `SAMRAIVectorReal` and calls weighted `L1/L2/max` norms using the original control-volume indices in `ibtk/src/solvers/impls/FACPreconditioner.cpp`;
  - in the live probe, those control-volume weights come from `HierarchyMathOps`, and level 1 covers the full coarse domain on this hierarchy;
  - as a result, the traced coarse-level norms are exactly zero even though the underlying coarse patch data are nonzero.
- Evidence:
  - the exported raw coarse solution vector remained nonzero:
    - `live_fac_probe_sol_l2_norm.level0 = 0.62157099509375002`
    - `/tmp/ibamr-live-fac-probe-true-n16-export/live_fac_probe_sol.level0.txt`
  - repeating the same live probe with explicit RT0 velocity transfer still gave `coarse_rhs(level 0) l2 = 0` in the FAC trace while producing a nonzero exported coarse solution:
    - FAC trace: `/tmp/ibamr-live-fac-probe-true-n16-rt0-fac-trace.txt`
    - level trace: `/tmp/ibamr-live-fac-probe-true-n16-rt0-level-trace.txt`
- Updated interpretation:
  - the zero `coarse_rhs`/`coarse_sol` entries in the FAC trace should not be used as evidence that the live FAC V-cycle skips the coarse solve on fully covered hierarchies;
  - this is a trace-norm bug or limitation, so future coarse-stage debugging should rely on raw exported level vectors or an unweighted trace path instead of the current weighted per-level FAC norms.

FAC trace fix:
- Patched `ibtk/src/solvers/impls/FACPreconditioner.cpp` so the per-level debug trace uses unweighted norms for its one-level trace views, i.e. it no longer reuses masked control-volume weights when reporting `coarse_rhs` and `coarse_sol`.
- Rebuilt the optimized parity target and reran the same matched live probe.
- The updated FAC trace now reports nonzero coarse stages:
  - `coarse_rhs(level 0) l2 = 2.77934465099973416`
  - `coarse_sol(level 0) l2 = 6.21681725294810139e-01`
- That coarse-solution trace now matches the raw exported coarse vector very closely:
  - trace `coarse_sol(level 0) l2 = 6.21681725294810139e-01`
  - raw export `live_fac_probe_sol.level0.txt` L2 norm `= 6.21570995093750023e-01`
- New trace artifact:
  - `/tmp/ibamr-live-fac-probe-true-n16-fac-trace.unweighted.txt`
- Updated interpretation:
  - the live FAC probe does produce a real nonzero coarse RHS and nonzero coarse solve on this matched smooth mode;
  - so coarse-path comparisons against the Julia reference are back on the table, and the old “coarse RHS is zero” conclusion should be considered retired.

Canonical Julia live-probe comparison:
- Added [compare_live_fac_probe_julia.jl](/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_live_fac_probe_julia.jl) to compare the repaired unweighted FAC trace against the canonical Julia reference V-cycle on the same matched smooth probe.
- Also updated the live probe export in [implicit_stokes_ib_operator_chain_parity_01.cpp](/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_operator_chain_parity_01.cpp) to write `problem_setup.txt` alongside the exported `live_fac_probe_rhs.level*.txt` / `live_fac_probe_sol.level*.txt` files.
- On the first pass, this comparison looked bad very early:
  - `y_pre(level 1)` IBAMR `1.55783973629987682` vs Julia `1.19195030548358893`
  - `residual_fine(level 1)` IBAMR `5.91163758435973552` vs Julia `1.20310321110469332`
- But that turned out to be a live-probe configuration mistake: the temporary input only set `coupling_aware_asm_closure_policy = RELAXED`, so the shell level solver silently stayed on the default `asm_subdomain_construction_mode = GEOMETRICAL` instead of `COUPLING_AWARE`.
- Evidence:
  - `StaggeredStokesPETScLevelSolver` defaults `d_asm_subdomain_construction_mode = GEOMETRICAL` in [StaggeredStokesPETScLevelSolver.h](/Users/boyceg/code/IBAMR/include/ibamr/StaggeredStokesPETScLevelSolver.h)
  - the incorrect live probe traced only `64` restricted subdomains on the fine level;
  - after adding
    - `asm_subdomain_construction_mode = "COUPLING_AWARE"`
    - `coupling_aware_asm_seed_axis = 0`
    - `coupling_aware_asm_seed_stride = 1`
    the corrected level trace showed `256` restricted subdomains, matching the Julia-style seed count.
- Corrected coupling-aware live probe result:
  - `rhs_input(level 1)` still matches to roundoff
  - `y_pre(level 1)` IBAMR `1.58431590811160872` vs Julia `1.19195030548358893`
  - `residual_fine(level 1)` IBAMR `2.95862987699915392` vs Julia `1.20310321110469332`
  - `coarse_rhs(level 0)` IBAMR `1.29458474257789513` vs Julia `5.28390467938935293e-01`
  - `coarse_sol(level 0)` IBAMR `5.10725612982514510e-01` vs Julia `2.08805847318780075e-01`
- So fixing the subdomain-construction mode removed a large artificial mismatch downstream of `y_pre`, but the real live shell-CAV path is still materially stronger than the canonical Julia reference on this smooth probe even before the coarse solve.
- Fine-level pre-smoother localization:
  - the corrected PETSc level trace contains three shell-apply blocks; the first block is the pre-smoother sweep, the second is an identical repeat, and the third is the post-smoother with a reduced RHS;
  - comparing that first block against a Julia replay of one canonical sweep, the first material difference is at:
    - `stage = local_delta`
    - `subdomain = 0`
    - `delta_rel ≈ 1.12420272054878975e-02`
  - At that point the incoming `restricted_rhs` still matches in norm, so the earliest visible difference is in the applied local patch solve/update, not in the restricted residual itself.
- Interpretation:
  - the shell level solver is now using the correct coupling-aware seed sweep, but the first real difference appears inside the very first local patch solve;
  - since [PETScLevelSolver.cpp](/Users/boyceg/code/IBAMR/ibtk/src/solvers/impls/PETScLevelSolver.cpp) already defaults `coupling_aware_vanka_matlab` local subsolves to `PCSVD`, the next suspect is not “LU versus SVD” by default, but some remaining difference in the local patch matrix/vector contract or in the exact PETSc `PCSVD` behavior relative to the canonical Julia pseudoinverse.
- One practical caveat:
  - the PETSc options file still reports several unused `ib_` and `stokes_ib_pc_level_sub_*` options in this probe mode, so the current probe should be read as “real shell CAV level smoother with in-process FAC tracing,” not a full outer-Krylov benchmark.

Live-probe comparison corrections and first real mismatch:
- The initial live-FAC local-patch comparison had three avoidable apples-to-oranges issues:
  - the canonical MATLAB/Julia docs have moved to the standalone solver repo:
    - `/Users/boyceg/code/implicit_ib_coupling_aware_vanka/docs/matlab-reference-algorithm-contract.md`
    - `/Users/boyceg/code/implicit_ib_coupling_aware_vanka/docs/matlab-reference-audit.md`
  - live probe vectors were being compared in raw IBAMR/PETSc ordering instead of canonical `(u, v, p)` ordering;
  - the live shell operator uses the IBAMR symmetric sign convention on the pressure row (`-div u = 0`), so local full-saddle blocks have to be compared after the corresponding pressure-row sign flip.
- Fixes made:
  - `live_fac_probe_export_dir` now also writes `eulerian_dof_map.level*.txt`;
  - [compare_live_fac_probe_julia.jl](/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_live_fac_probe_julia.jl) now:
    - remaps live probe vectors into canonical Julia ordering,
    - applies the IBAMR pressure-row sign convention to the reference hierarchy operators,
    - replays the fine-level sweep in the remapped IBAMR seed order instead of assuming `1,2,3,...,N^2`;
  - [compare_live_fac_local_dump_julia.jl](/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_live_fac_local_dump_julia.jl) now:
    - remaps targeted patch dumps via the exported DOF map,
    - compares the dumped local block against the pressure-row-aligned canonical block,
    - replays the exact dumped local matrix/RHS with the Julia pseudoinverse,
    - and, when available, compares the dumped `y_before` / full residual against the remapped reference residual path.
- What these corrected comparisons show:
  - Subdomain `0` on the first pre-smoother block is not a real mismatch.
    - After remap, the patch set matches exactly.
    - After pressure-row sign alignment, the local block matches to roundoff:
      - `same_dofs_local_A_pressure_rowsigned_rel = 1.41721840393849903e-16`
    - The local RHS matches exactly.
    - Replaying the exact dumped local matrix/RHS with the Julia pseudoinverse matches the dumped PETSc local delta to roundoff:
      - `dump_matrix_local_delta_replay_rel = 1.41727891251884138e-15`
    - The apparent subdomain-0 difference was therefore a comparison-contract bug, not a solver bug.
  - Subdomain `1` is also not a real mismatch.
    - The patch set matches as a set after remap.
    - The full remapped residual matches to roundoff:
      - `global_residual_full_rel = 4.43059490744281945e-17`
    - The pressure-row-aligned local block matches to roundoff:
      - `same_dofs_local_A_pressure_rowsigned_rel = 1.41721840393849903e-16`
    - The local RHS and local delta both match to roundoff:
      - `same_dofs_local_rhs_rel = 2.17379026897365338e-16`
      - `same_dofs_local_delta_pressure_rowsigned_rel = 2.16075942959057541e-15`
    - So the earlier “first real mismatch at subdomain 1” report was also a comparison-order artifact.
  - After the remap and sign fixes, the first material live fine-level mismatch moves out to:
    - `stage = local_delta`
    - `subdomain = 7`
    - from the corrected live comparator:
      - `y_rel = 6.74640927432948586e-04`
      - `delta_rel = 2.02515278898753731e-02`
  - The targeted subdomain-7 dump confirms that this is the first genuinely interesting patch:
    - patch set matches as a set after remap (`dump_only_count = 0`, `canonical_only_count = 0`);
    - exact PETSc local solve replay is still fine:
      - `dump_matrix_local_delta_replay_rel = 1.89644207031734428e-15`
    - but the incoming full residual is now genuinely different:
      - `global_residual_full_rel = 6.83469181831082295e-02`
      - `global_residual_full_max_abs = 4.29001899289448319e-01`
    - the restricted local RHS differs accordingly:
      - `same_dofs_local_rhs_rel = 2.49877012818236315e-02`
      - `same_dofs_local_rhs_max_abs = 3.33950137308061112e-02`
    - and even after pressure-row sign alignment the local block is no longer identical:
      - `same_dofs_local_A_pressure_rowsigned_rel = 1.51405343059465760e-02`
      - `same_dofs_local_A_pressure_rowsigned_max_abs = 3.12647583575426324e+00`
    - so the resulting local delta differs by a few percent:
      - `same_dofs_local_delta_pressure_rowsigned_rel = 3.85237306721021577e-02`
      - `same_dofs_local_delta_pressure_rowsigned_max_abs = 2.00350480830713967e-03`
    - entry-level localization of the pressure-row-aligned local-block difference shows only four differing entries, all in the `v-v` block:
      - `v(7,2), v(7,2)`: dump `46.416025524462704`, reference `43.289549688708441`, abs diff `3.126475835754263`
      - `v(6,2), v(6,2)`: dump `43.862953505266283`, reference `42.647892704740990`, abs diff `1.215060800525293`
      - `v(6,2), v(7,2)`: dump `-1.1859809559723606`, reference `-2.2146713492312244`, abs diff `1.0286903932588638`
      - `v(7,2), v(6,2)`: dump `-1.1859809559723609`, reference `-2.2146713492312244`, abs diff `1.0286903932588636`
    - block summary:
      - all `p-*`, `*-p`, `u-*`, and `*-u` blocks agree to roundoff;
      - the entire local-block discrepancy lives in the `v-v` block (`fro ≈ 3.6561772045258745`, `4` differing entries).
    - interpretation of that structure:
      - this does not look like a general pressure-nullspace or saddle coupling mismatch;
      - it also does not look like a pure Stokes discretization mismatch, since the rest of the local velocity-pressure stencil matches;
      - the remaining local-block discrepancy is concentrated in a tiny vertical-velocity subblock near the immersed structure and is therefore most consistent with an IB coupling / `SAJ` contribution mismatch on that patch.
- Updated interpretation:
  - the earliest live-FAC divergence is no longer “the first local patch solve” in general;
  - after correcting DOF order, seed order, and sign convention, the first several patches are consistent with the canonical Julia reference;
  - the first real mismatch on this smooth finest-`16`, depth-2 Julia-density probe is at subdomain `7`, and it includes both:
    - a residual-path mismatch before the local solve, and
    - a smaller but real local-block/operator mismatch on that same patch;
  - the local `PCSVD` backend itself is not the culprit there, since replaying the dumped exact matrix/RHS reproduces the dumped PETSc local delta to machine precision.

Normalization split and corrected live-shell comparison:
- Follow-up matrix forensics on the four subdomain-7 `v-v` entries showed that the explicit bridge-exported finest-level operators were already consistent with the canonical Julia reference:
  - `/tmp/ibamr-bridge-juliadensity-n16-export/L.level1.mtx` matched the canonical Julia `K = 1e2` operator to roundoff on all four suspect entries;
  - `/tmp/ibamr-bridge-juliadensity-n16-export/SAJ.level1.mtx` also matched the canonical Julia `SAJ` contribution to roundoff;
  - the pure Stokes remainder `L - SAJ` matched the `K = 0` Julia operator to roundoff.
- The live shell subdomain-7 dump, however, matched a different operator exactly:
  - the dumped patch matrix was *exactly* equal to
    - `Stokes + (1 / ds_theta) * SAJ`
    - with `1 / ds_theta = 3.9788735772973833`;
  - numerically:
    - `dump vs bridge L rel = 1.5140534305946611e-02`
    - `dump vs bridge (L - SAJ + (1/ds) * SAJ) rel = 0`
- Interpretation:
  - the earlier Julia-density live-probe mismatch was not exposing a local shell-matrix bug inside subdomain extraction;
  - it was exposing a normalization mismatch between:
    - the explicit parity/bridge export path, which can represent the canonical Julia-density contract directly, and
    - the live shell/FAC solver path, which uses IBAMR's native nodal-force scaling.

Comparator correction for native live-shell runs:
- Patched the live Julia comparators so they now recover the effective canonical Julia stiffness from `problem_setup.txt` using `spring_normalization`:
  - `JULIA_DENSITY`: `K_reference = spring_pair_k * ds_theta^2`
  - `IBAMR_NODAL`: `K_reference = spring_pair_k * ds_theta`
- Files updated:
  - `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_live_fac_probe_julia.jl`
  - `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_live_fac_local_dump_julia.jl`
- This makes the native live shell probe comparable to the canonical Julia reference on the same effective stiffness without pretending that the live solver path uses the explicit Julia-density matrix contract internally.

Refreshed IBAMR_NODAL live-shell result:
- Reran the optimized live FAC probe with the corrected native contract:
  - input: `/tmp/implicit_stokes_ib_operator_chain_parity_01.live_fac_probe.true_n16.coupling.input`
  - export: `/tmp/ibamr-live-fac-probe-true-n16-export`
  - FAC trace: `/tmp/ibamr-live-fac-probe-true-n16-coupling-fac-trace.refresh.txt`
  - level trace: `/tmp/ibamr-live-fac-probe-true-n16-coupling-level-trace.refresh.txt`
- The refreshed export now includes:
  - `eulerian_dof_map.level0.txt`
  - `eulerian_dof_map.level1.txt`
- Corrected stage comparison against Julia:
  - `rhs_input(level 1)` matches to roundoff
  - `y_pre(level 1)` relative norm difference `≈ 1.9172303326790371e-03`
  - `residual_fine(level 1)` relative norm difference `≈ 2.0702573423548882e-03`
  - `y_post(level 1)` relative norm difference `≈ 5.1334554672441369e-03`
  - first material fine-level difference moved out to:
    - `stage = local_delta`
    - `subdomain = 37`
    - `delta_rel ≈ 6.7882935991718366e-02`
- This is a major contraction of the live-shell gap relative to the mis-normalized Julia-density comparison.

Targeted native subdomain-37 dump:
- Dumped the first remaining nodal mismatch patch:
  - dump dir: `/tmp/ibamr-live-fac-local-dump-true-n16-sub37/level1_call1_apply1_subdomain37`
- The result is much cleaner than the old Julia-density subdomain-7 report:
  - on the shared DOFs, the pressure-row-aligned local block matches to roundoff:
    - `same_dofs_local_A_pressure_rowsigned_rel = 2.3693806002248539e-15`
  - the full remapped residual matches to roundoff:
    - `global_residual_full_rel = 1.6286051405832912e-15`
  - the local RHS matches to roundoff:
    - `same_dofs_local_rhs_rel = 3.0823630282182501e-15`
  - the local delta matches to roundoff:
    - `same_dofs_local_delta_pressure_rowsigned_rel = 5.9976426841372485e-15`
  - exact PETSc local solve replay still matches to roundoff:
    - `dump_matrix_local_delta_replay_rel = 5.4261330848101306e-15`
- The remaining mismatch at subdomain `37` is therefore no longer an operator or residual mismatch on the shared patch; it is a patch-set mismatch only:
  - `dump_patch_size = 127`
  - `canonical_patch_size = 118`
  - `dump_only_count = 9`
  - `canonical_only_count = 0`
  - `dump_only_sample = [114, 130, 146, 354, 370, 386, 610, 626, 642]`
- Inference from the canonical ordering:
  - those extra DOFs form a small `3 x 3` `(u, v, p)` strip near the seed location rather than a random set of indices.

Relaxed closure debug on subdomain 37:
- Used the built-in coupling-aware subdomain debug hook:
  - `IBAMR_CAV_SUBDOMAIN_DEBUG_FILE=/tmp/ibamr-cav-subdomain-debug-sub37.txt`
  - `IBAMR_CAV_SUBDOMAIN_DEBUG_LEVEL=1`
  - `IBAMR_CAV_SUBDOMAIN_DEBUG_SUBDOMAIN=37`
  - `IBAMR_CAV_SUBDOMAIN_DEBUG_POLICY=RELAXED`
- Important result:
  - for this subdomain, `overlap_dofs` and `closure_dofs` are identical in the debug dump;
  - the 9 extra DOFs are therefore **not** coming from the final RELAXED-mode reinsertion of `initial_velocity_dofs`;
  - they are already present in the closure generated from the involved-cell / cell-closure map.
- Updated interpretation:
  - the remaining native live-shell discrepancy on this aligned smooth probe is now localized to relaxed patch-construction semantics, not:
    - local `PCSVD`,
    - residual replay,
    - pressure sign convention,
    - or local matrix assembly on the shared patch.
  - The next code-level suspect is the cell-closure map / involved-cell selection path in
    `/Users/boyceg/code/IBAMR/src/navier_stokes/StaggeredStokesPETScMatUtilities.cpp`,
    especially the logic that maps one seed component to the relaxed closure around subdomain `37`.

Numerical-zero couplings in the seed stencil:
- Follow-up inspection of the live dumped seed row for the old first-mismatch nodal patch (`subdomain 37`, seed raw DOF `74`) found that the row was carrying many floating-point crumbs that were being treated as structural couplings:
  - stored seed-row nonzeros: `54`
  - entries with `|value| < 1e-10`: `23`
  - entries with `|value| < 1e-12`: `23`
  - smallest magnitude: `3.0512048810393941e-21`
- Representative tiny entries from the live row dump:
  - `v(4,6) = 3.0512048810393941e-21`
  - `v(5,2) = -2.7418310770260616e-19`
  - `v(7,6) = 3.3579970416151008e-19`
  - `v(6,2) = -6.6966313538529049e-18`
  - `v(7,2) = -9.8348454493508200e-18`
- The relaxed patch builder was still using an exact-zero test in
  `/Users/boyceg/code/IBAMR/src/navier_stokes/StaggeredStokesPETScMatUtilities.cpp`
  when expanding seed rows into `initial_velocity_dofs`, so all of these numerical crumbs were being treated as genuine one-step A00 neighbors.

Tolerance fix and immediate effect:
- Patched `build_initial_velocity_set_from_seed_components()` in
  `/Users/boyceg/code/IBAMR/src/navier_stokes/StaggeredStokesPETScMatUtilities.cpp`
  to ignore row entries smaller than a row-scaled numerical-zero tolerance:
  - `tol = max(ncols * eps * row_max_abs, 1e-14 * row_max_abs)`
- Rebuilt the optimized parity executable and reran the same aligned `IBAMR_NODAL` live probe.
- Immediate result:
  - the first material fine-level mismatch moved from `subdomain 37` to `subdomain 58`;
  - corrected stage differences on the fine level shrank again:
    - `y_pre(level 1)` relative difference `≈ 2.9996791067152731e-05`
    - `residual_fine(level 1)` relative difference `≈ 4.3864184467526229e-06`
    - `y_post(level 1)` relative difference `≈ 4.6864038677732027e-03`
- Updated interpretation:
  - yes, there were real floating-point near-zeros in the seed stencil that should be excluded from relaxed patch construction;
  - filtering those numerical crumbs materially improved live-shell parity and pushed the first remaining patch-level mismatch much later in the sweep.

Near-zero false alarm at the next reported subdomain:
- Dumped the new first reported local mismatch (`subdomain 58`) after the numerical-zero fix.
- Result:
  - on the shared DOFs, the pressure-row-aligned local block still matches to roundoff;
  - the remapped full residual still matches to roundoff;
  - the pressure-row-aligned local delta still matches to roundoff;
  - the apparently large relative local-delta difference there was only because both local updates were at machine scale:
    - dumped local-delta norm `≈ 4.07e-16`
    - canonical local-delta norm `≈ 4.69e-16`
- Conclusion:
  - subdomain `58` was not a meaningful remaining local-patch mismatch;
  - the live-probe comparator was still too eager to classify machine-epsilon local updates as “material”.

Comparator tightening for local trace diagnostics:
- Patched
  `/Users/boyceg/code/IBAMR/tests/external-reference/cav-stokes-ib/reference-tools/compare_live_fac_probe_julia.jl`
  so `first_material_level_difference` now requires:
  - a relative difference larger than `1e-2`, and
  - an absolute difference larger than `1e-12`
  before classifying a fine-level local-delta or local-state difference as material.
- After rerunning the aligned nodal live-probe comparison on the post-fix export:
  - `first_material_level_difference = none`
- Interpretation:
  - there is no longer a meaningful fine-level pre-smoother mismatch on this aligned smooth probe.

RT0 transfer check after the fine-level fixes:
- Built a corrected RT0 live-probe variant from the current coupling-aware nodal input by adding:
  - `U_prolongation_method = "RT0_REFINE"`
  - `U_restriction_method = "RT0_COARSEN"`
- Reran the same live FAC probe and compared it against the canonical Julia reference.
- Result:
  - the remaining coarse-stage gap nearly vanished:
    - `coarse_rhs(level 0)` relative difference dropped from `≈ 1.85e-01` to `≈ 5.82e-05`
    - `coarse_sol(level 0)` relative difference dropped from `≈ 3.04e-02` to `≈ 1.86e-04`
    - `y_after_correction(level 1)` relative difference dropped to `≈ 2.66e-05`
    - `y_post(level 1)` relative difference dropped to `≈ 1.97e-05`
  - the level trace still reports:
    - `first_material_level_difference = none`
- Updated interpretation:
  - after fixing numerical-zero patch growth, the fine-level relaxed sweep is effectively aligned with the canonical Julia reference on this smooth probe;
  - the remaining multilevel discrepancy was overwhelmingly a transfer-semantics issue;
  - explicit RT0 velocity transfer makes the live FAC V-cycle nearly match the Julia reference end-to-end on this test.

## Remaining Gaps

1. Verify the actual high-stiffness IBAMR runs under discussion are using the shell `coupling_aware_vanka_matlab` family when they are being compared against the audited MATLAB reference, rather than the much weaker default `ASM` path.

2. Diagnose/fix the two-level hang in `implicit_stokes_ib_operator_chain_parity_01`, since multilevel/V-cycle parity is now the most important open semantic check for convergence robustness.

3. Finish coarse-level / V-cycle parity on reproduced stiff shell-CAV cases:
- restriction/prolongation,
- FAC correction transfer method choices versus MATLAB-style RT0 velocity transfer,
- coarse solve semantics,
- coarse pressure-nullspace handling,
- post-correction post-smoothing.

4. Determine exactly what local and coarse solve contracts the shell multilevel path is actually using in code on stiff runs:
- true default SVD/pseudoinverse semantics,
- LU-like solves,
- or some mixed contract that differs from the audited MATLAB reference.

5. Treat the remaining one-level `1e-12` semantic-replay drift as secondary unless it can be shown to correlate with outer-iteration growth or divergence.
