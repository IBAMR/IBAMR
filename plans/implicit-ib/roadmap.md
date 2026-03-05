# Implicit IB Roadmap (High-Level Milestones)

## Summary
Focus on a staged roadmap to stabilize implicit IB, validate algorithmic correspondence, improve performance, and expand capability (parallel + kernel generalization).
Use deterministic CMake tests as the main confidence signal, but treat milestone checks as directional goals rather than strict gates.

## Roadmap Record
- Store and maintain roadmap at: `plans/implicit-ib/roadmap.md`
- Keep this file updated with:
1. Current reference configs
2. Milestones and status
3. Key comparison notes vs MATLAB
4. PR history and follow-up items

---

## Current State (2026-03-04)
1. Development is on branch `codex/implicit-ib-milestone-3`.
2. Milestone 1 remains under review in PR `#1887`.
3. Milestone 2 implementation work is complete on this branch:
   - PR-A utility and deterministic test work is complete.
   - PR-B solver wiring and integration test/smoke coverage is complete.
4. Reference implicit stabilization changes are applied in:
   - `src/IB/IBImplicitStaggeredHierarchyIntegrator.cpp`
5. A promoted finished reference implicit test is added in:
   - `tests/IB/implicit_stokes_ib_01.cpp`
   - `tests/IB/implicit_stokes_ib_01.baseline.input`
   - `tests/IB/implicit_stokes_ib_01.baseline.output`
   - `tests/CMakeLists.txt` (test registration)
6. Focused validation passing:
   - `attest -R 'IB/implicit_stokes_ib_01.baseline.input'`
7. Expanded validation on the milestone-1 branch is passing:
   - `attest -R '^IB/'` (all `IB` tests in current CMake build subset).
8. M2 solver-mode tests and transfer-invariant tests are passing in CMake and autotools builds.
9. M3 test work is complete for strict-CAV/classical-Vanka equivalence, RIDP sum-invariant checks, RT0 velocity transfer parity, and pressure transfer parity in 2D/3D.
10. End-to-end FAC runtime transfer-configuration validation is tracked in Milestone 4.

---

## Milestone 1 — Cleanup and stabilization
### Goal
Clean up implicit IB implementation and land a focused stabilization PR that restores reference implicit workflows and finished test coverage.

### Scope
1. Remove accidental/debug artifacts that interfere with reference implicit solves.
2. Stabilize the `IBImplicitStaggeredHierarchyIntegrator` execution path and required INS-side plumbing.
3. Promote core reference implicit tests from unfinished/development status into finished attest-managed tests.
4. Keep this milestone strictly reference-oriented (no new smoother feature development).

### Status
- In progress (PR open: `#1887`).
- Completed so far:
1. Reference stabilization edits applied for `IBImplicitStaggeredHierarchyIntegrator`.
2. Reference finished implicit attest test added and passing in focused runs.
3. Roadmap file moved to top-level `plans/implicit-ib/`.

### Explicitly Deferred from M1
1. Coupling-aware vs geometrical ASM behavior work.
2. New smoother modes/policies and tuning keys.
3. Transfer-method feature expansion (RT0/FAC/PETSc tuning work).
4. MATLAB parity/performance studies for coupling-aware smoothers.

---

## Milestone 2 — Coupling-Aware Vanka (CAV) smoother implementation
### Goal
Implement coupling-aware Vanka smoother construction for implicit Stokes-IB solves using A00-sparsity-driven subdomain construction.

### Scope
1. Construct CAV subdomains from operator sparsity in a deterministic way.
2. Keep geometric and CAV subdomain modes clearly separated (one or the other).
3. Preserve required ASM overlap/nonoverlap structure for PETSc usage.
4. Validate serial behavior first with focused deterministic tests and smoke runs.
5. Keep implementation general enough for later parallel extension.

### Status
- Complete on this branch.

### Delivery Plan (Stacked PRs)
1. PR-A (base: M1):
   - CAV subdomain construction utilities in `StaggeredStokesPETScMatUtilities`.
   - Deterministic serial tests for sparsity-to-cells, cells-to-closure, strict/relaxed closure, and overlap/nonoverlap invariants.
   - MATLAB parity audit document and mapping table for `extract_coupled_dofs(.m/.m2)` semantics.
2. PR-B (stacked on PR-A):
   - Solver wiring and input-driven mode selection for geometrical vs CAV construction.
   - Serial integration/smoke checks and logging consistency updates.
   - Status: complete on branch.

---

## Milestone 3 — Stokes infrastructure + RT0/RIDP transfer verification (without IB coupling)
### Goal
Verify that Stokes-only infrastructure corresponds to basic Vanka-style multigrid behavior and validate key RT0/RIDP transfer invariants.

### Scope
1. Define reference Stokes-only cases.
2. Compare IBAMR behavior with MATLAB reference implementations.
3. Document where strict CAV should match classical Vanka and where behavior intentionally differs.
4. Add focused tests/diagnostics for transfer invariants.
5. Check both FAC transfer settings and PETSc-side transfer settings.
6. Verify `R * Id * P` behavior:
   - row sums are `1` for velocity DOFs and `0` for pressure DOFs,
   - column sums are `1` for velocity DOFs and `0` for pressure DOFs.

### Status
- Complete on branch (`codex/implicit-ib-milestone-3`).
- Completed:
1. Strict-CAV/classical-Vanka equivalence tests added in 2D/3D.
2. RIDP sum-invariant tests added in 2D/3D.
3. RT0 velocity prolongation/restriction parity tests (matrix-based PETSc vs matrix-free SAMRAI) added in 2D/3D, with affine, piecewise-RT0, and nonlinear profiles.
4. Pressure transfer parity tests (`CONSTANT_REFINE`/`CONSERVATIVE_COARSEN`) added in 2D/3D, with affine and nonlinear profiles.
5. Transfer parity tests now include nontriviality checks and SAMRAI-vs-PETSc max-norm consistency checks for corresponding fields.

---

## Milestone 4 — MATLAB parity for Stokes + IB
### Goal
Build exact-comparison cases between IBAMR and MATLAB for Stokes+IB using programmatic structure generation, then identify and reduce discrepancies.

### Scope
1. Align geometry/parameters between codebases.
2. Compare subdomains, iteration traces, and residual behavior.
3. Track mismatches with likely causes and resolution plan.
4. Add an end-to-end FAC runtime transfer-configuration validation case that complements Milestone 3 operator-level parity tests.

### Status
- Not started.

---

## Milestone 5 — Bespoke Gauss-Seidel performance work
### Goal
Improve performance of the custom GS-like smoother while preserving intended solver semantics.

### Scope
1. Identify dominant runtime costs.
2. Apply targeted optimizations in small PRs.
3. Re-check behavior against Stokes-only and Stokes+IB reference cases.

### Status
- Not started.

---

## Milestone 6 — Parallel correctness and robustness
### Goal
Extend current serial-validated behavior to parallel runs with consistent subdomain and solver behavior.

### Scope
1. Validate coverage/ownership behavior across ranks.
2. Add representative MPI test scenarios.
3. Compare serial vs parallel outcomes on selected reference cases.

### Status
- Not started.

---

## Milestone 7 — Kernel generalization beyond `IB_4`
### Goal
Generalize implementation to support additional IB kernel functions while preserving coupling-aware logic and solver correctness.

### Scope
1. Parameterize kernel-dependent paths.
2. Add representative tests with non-`IB_4` kernels.
3. Document kernel-specific behavior differences where relevant.

### Status
- Not started.

---

## PR Strategy
Flexible batching, with clear milestone intent:
1. Each PR should have one primary milestone objective.
2. Keep scope explicit (what it does and does not do).
3. Keep roadmap file updated after each merged PR.

## Current Reference Configs
1. Build system: CMake + attest.
2. Validation focus (current):
   - `cmake --build <build_dir> -j10 --target tests-IB_implicit_stokes_ib_01`
   - `<build_dir>/attest --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R 'IB/implicit_stokes_ib_01.baseline.input'`
3. Reference implicit test intent:
   - run one or two implicit timesteps without crashes,
   - verify stable deterministic output under attest,
   - avoid coupling-aware/classical-Vanka feature coverage in Milestone 1.

## Key Comparison Notes vs MATLAB
- Deferred until Milestones 3 and 4.
- Current milestone does not treat MATLAB parity as a gate.

## PR History and Follow-up Items
- Add entries as PRs are opened/merged.
- Suggested entry format:
  - PR: `<number/link>`
  - Milestone: `<M#>`
  - Summary: `<what changed>`
  - Validation: `<tests/runs>`
  - Follow-up: `<next tasks>`

- Current active PR:
  - PR: `https://github.com/IBAMR/IBAMR/pull/1887`
  - Milestone: `M1`
  - Summary: stabilize implicit hierarchy integrator path and promote a finished implicit attest test
  - Validation: focused `IB/implicit_stokes_ib_01` runs and full `IB` subset attest runs
  - Follow-up: merge milestone-1 branch, then begin milestone-2 work on stacked branch
