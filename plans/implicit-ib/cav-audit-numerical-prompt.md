# Audit Lane N: Numerical Reproduction and Parity

**Prompt version:** `cav-audit-v4/N` (2026-08-20)

Use with `cav-independent-audit-prompt.md`. This lane is independent and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

## Question

Does the candidate numerically reproduce the declared live IBAMR problem and mathematically specified algorithm, with the first discrepancy localized rather than hidden by aggregate norms?

## Required evidence

### First-mismatch ladder

Walk these stages in order and stop promotion at the first unresolved mismatch:

1. configuration, geometry, coefficients, scales, boundary conditions, stiffness, dimension, and MPI size;
2. live original operator and named blocks;
3. DOF mapping, ownership, sign conventions, nullspace, and pressure gauge;
4. formal pressure-seed `E_h` provenance plus patch seed/order/membership/restricted sets, with legacy full-`A00` velocity-seed data kept separately labeled;
5. local matrices, RHS vectors, factorization status, solves, and backward errors;
6. per-patch original residual, correction, and state after mapping;
7. RAS macro geometry, owned/solve partitions, macro-local patch order/state, completed correction, restriction, and outer assembly;
8. transfers, synchronization, coarse RHS/solve/correction, and gauge;
9. smoother output;
10. V-cycle stages/output;
11. FGMRES on the original operator and original momentum/divergence/IB residuals.

Use provenance-bearing observable artifacts from live IBAMR objects and the separate detached oracle worktree.

For every matrix-builder-supported regularized-delta kernel, treat the original Fortran interaction implementations in `lagrangian_interaction2d.f.m4`, `lagrangian_interaction3d.f.m4`, and their included `lagrangian_delta.f.m4` definitions as the production source of truth. Require the focused 2D/3D matrix/live interpolation-and-spreading evidence before accepting an assembled IB operator. A C++/Fortran disagreement is not resolved by editing the matrix-free path unless a separate focused test independently reproduces a live-kernel bug.

### Mandatory matrix

All listed `K` values mean `K = 1, 10^2, 10^4`. `K=10^6` is an additional conditioning case; `K=10^8` remains diagnostic.

| ID | Grid/dimension/ranks | Construction/composition/traversal | Backend | Required endpoint |
|---|---|---|---|---|
| N0 | live periodic `4x4` 2D, 1 rank, `K=1` | RELAXED and STRICT; forward | Eigen reference | raw export, hand mapping, selected patch/local solve/comparator validation |
| N1 | single-level periodic `8x8` 2D, 1 rank, all listed `K` | RELAXED and STRICT; forward and reverse global multiplicative | Eigen reference | operator through every per-patch state and smoother |
| N2 | levels `8x8 -> 16x16` 2D, 1 rank, all listed `K` | RELAXED; forward pre-sweep, reverse post-sweep, and explicit palindrome | Eigen reference | transfers through V-cycle and original FGMRES residuals |
| N3 | levels `8x8 -> 16x16 -> 32x32` 2D, 1 rank, all listed `K` | RELAXED production multiplicative policy | production backend | M5/cav-6b V-cycle, FGMRES, physical residuals |
| N4 | single-level `8x8` 2D, 1 rank, `K=1,10^4` | RELAXED forward | Eigen reference, production Eigen/Schur, recovered BLAS/LAPACK LU and SVD | identical local matrices; solution/backward-error and sweep cross-check |
| R0 | single-level periodic `8x8` 2D, 1 rank, all listed `K` | pressure macro blocks `4x4`; standard ghost 0; IB ghost 0 and 2; RELAXED forward, plus STRICT at `K=1` | Eigen reference | every inner patch state, completed macro correction, restriction, outer assembly |
| R1 | levels `8x8 -> 16x16` 2D, 1 rank, all listed `K` | macro blocks `4x4`; standard ghost 0; IB ghost 2; forward and palindrome | Eigen reference | RAS smoother and V-cycle |
| R2 | levels `8x8 -> 16x16 -> 32x32` 2D, 1 rank, all listed `K` | macro blocks `8x8`; standard ghost 0; IB ghost 2; production traversal | production backend | M9/cav-9b V-cycle, FGMRES, physical residuals |
| P2 | single-level `16x16` 2D, 2 ranks, `K=10^4` | applicable production mode | production backend | native ownership, construction, lifecycle, transfer, and rank-consistency checks; sandbox parity only if advertised |
| P3 | single-level `8^3` 3D, 1 rank, `K=1` | applicable native mode | production backend | Debug/Release build and native smoke/lifecycle checks; no claim of 3D sandbox parity |

N0–N4 apply through M5 as their capabilities become available. R0–R2 apply through M9. P2/P3 are production-robustness checks, not substitutes for 2D rank-1 oracle parity. A candidate advertising 3D rank-2 or distributed RAS adds those configurations; otherwise unsupported requests must fail clearly.

### Two mandatory numerical tracks

Audit these tracks separately and report both. Neither can substitute for the other.

1. **Common-arithmetic algorithm parity:** map the independently exported candidate and oracle structural artifacts into canonical order. In one documented replay runtime, execute candidate semantics using candidate-exported ordered patch/macro/transfer artifacts and oracle semantics using oracle-exported artifacts. Run the paired replay on the mapped live-IBAMR operator/RHS and again on the independently exported oracle operator/RHS. Both replays in each pair must use the same local-solve routine, multiplication/traversal order, gauge operation, and floating-point runtime. Apply the fixed table below at every applicable layer in both common-input runs.
2. **Live production behavior:** rerun the actual IBAMR matrix-free/original operator, recovered local backend, multiplicative or RAS shell, FAC/V-cycle, and FGMRES path. Require fresh original-residual recomputation, local backward errors, incremental/fresh consistency, convergence reason, iteration history, and separate momentum/divergence/IB consistency metrics. For mandatory low-`K` rows, require the recomputed live original relative residual to be at most the declared `1e-10` solve tolerance; compare the scalar vector `(rho, ||r_m||_2, ||r_d||_2)` and iteration history with the oracle run under the fixed per-`K` table. Require the live-only coupling check to report both IB action norms and satisfy `rho_IB <= 1e-10`; do not compare it to a nonexistent oracle third residual row.

Retain the independently constructed raw candidate/oracle operators and raw entrywise live-versus-MATLAB corrections/residuals as diagnostics. Localize their first mismatch, but do not use vendor/runtime entrywise identity as the algorithm-parity definition once the common-arithmetic track is valid. Conversely, the common replay cannot establish construction correctness: exact mapped structure/order, Fortran-authoritative matrix/live checks, independently exported raw operators, and mapping/comparator controls remain mandatory.

### Exact metrics and fixed low-K targets

For mapped vectors `x` and `y`, define

`E_inf(x,y) = ||x-y||_inf / max(1.0, ||x||_inf, ||y||_inf)`.

The denominator floor is exactly `1.0`; therefore sub-unit data are checked by absolute max error. For mapped matrices with exactly equal dimensions and sparsity, vectorize entries in canonical row/column order and use the same formula. Report the unscaled absolute max error too.

For a local solve `A_i delta_i = r_i`, report

`eta_i = ||A_i delta_i-r_i||_inf / (||A_i||_inf ||delta_i||_inf + ||r_i||_inf + 1e-30)`.

For incremental/fresh residual consistency, use

`E_r = ||r_inc-(b-Ax)||_inf / max(1.0, ||r_inc||_inf, ||b-Ax||_inf)`.

For original physical residuals, the reduced Eulerian state is `x=(u,p)` and has no independent third IB row block. Report absolute norms and

`rho = ||b-Ax||_2 / max(||b||_2, 1e-30)`,

`rho_m = ||R_u(b-Ax)||_2 / max(||b_m||_2, sqrt(eps) ||b||_2, 1e-30)`,

`rho_d = ||R_p(b-Ax)||_2 / max(||b_d||_2, sqrt(eps) ||b||_2, 1e-30)`.

Independently evaluate the IB contribution from the original matrix-free live Jacobian and the FAC-owned live `d_SAJ_mat`:

`c_IB^mf = R_u[(A_live-A_Stokes)x]`, `c_IB^mat = R_u A_IB^FAC x`, `r_IB = c_IB^mf-c_IB^mat`,

`rho_IB = ||r_IB||_2 / max(||c_IB^mf||_2, ||c_IB^mat||_2, sqrt(eps)||b||_2, 1e-30)`.

Report both IB action norms and the absolute/relative coupling defect. Do not label the elasticity action itself as a third residual equation. Only report a Lagrangian kinematic/constraint residual when the candidate actually carries an independent corresponding unknown; do not fabricate one for the eliminated velocity path.

Apply these fixed targets at every applicable common-arithmetic parity layer, on both the live-operator common input and the oracle-operator common input. Apply the same table to the live scalar residual-summary vector and iteration-history comparison:

| K | Expected mapped agreement | Mandatory disposition |
|---:|---:|---|
| `1` | `E_inf <= 1e-13` | near-roundoff expectation; a larger discrepancy blocks PASS |
| `10^2` | `E_inf <= 1e-12` | fixed low-K criterion |
| `10^4` | `E_inf <= 2e-11` | fixed `1e-11`-scale low-K criterion; every final output also has hard ceiling `1e-10` |
| `10^6` | conditioning-aware only with measured condition estimate and backward error | never relax structure/gauge/residual reporting |
| `10^8` | diagnostic | no gate until a dated plan decision defines one |

The `1e-10` final ceiling is not permission to waive the per-K targets. A common-arithmetic or live scalar-summary target miss is **FAIL** when reproduced and **INCOMPLETE** when the cause/evidence cannot be established. A raw cross-runtime entrywise miss must be reported and localized; it is a defect only when it exposes incorrect structure/formula/semantics, a comparator/mapping failure, or failed live-production criteria.

### Export, mapping, and comparator validation

Complete every raw-export, tiny live-case, comparator control, and permutation/sign/gauge/omission/reordering mutation in the common protocol. The legal gauge mutation must change only pressure by a constant and must map to agreement; an equivalent velocity shift must fail. Preserve raw and mapped hashes.

### Production robustness

Own the numerical conclusions for Debug/Release 2D/3D build/smoke, supported MPI ranks, local-solver failure behavior, NaN/Inf checks, and boundary/periodic/tiny/no-IB cases. Use D/C evidence for lifecycle, cache, and resources. Unsupported scope must be rejected explicitly.

## Exclusions

- Do not accept old logs as current validation.
- Do not substitute surrogate production matrices/operators.
- Do not skip an earlier mismatch because a later aggregate norm passes.
- Do not use conditioning to excuse any `K <= 10^4` common-arithmetic or live scalar-summary discrepancy.
- Do not let either common-arithmetic replay rebuild candidate patch/macro/transfer artifacts with oracle construction code.
- Do not treat a common replay as evidence that the actual live backend converged or that the raw operator was constructed correctly.
- Do not judge architecture or speculative performance.
- Do not edit or fix the candidate, oracle, or comparator.

## Pass/fail

**PASS** requires every mandatory-matrix row applicable to the milestone, both two-source common-arithmetic runs, the separate live-production run, validated mapping/comparator controls, exact structural agreement, Fortran-authoritative matrix/live kernel evidence, per-K targets, physical residual components, supported-scope robustness, and activated sensitivity detection. A reproduced semantic/structural/formula mismatch or criterion violation is **FAIL**. Missing/invalid common replay, live run, matrix rows, raw artifacts, comparator controls, environment, or mutation activation is **INCOMPLETE**.

## Lane-specific report schema

After the common schema, provide: mandatory-matrix row status; configuration-by-layer artifact table; hashes for raw candidate/oracle exports, each mapped artifact, common-replay implementation, and comparator; separate results for live-input replay, oracle-input replay, and actual live execution; formulas and measured absolute/scaled errors; targets; local backward errors; earliest raw and common-arithmetic mismatches; original physical residual triplet; robustness results; sensitivity controls; highest layer reached; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
