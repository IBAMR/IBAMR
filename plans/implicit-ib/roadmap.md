# Implicit IB Roadmap (High-Level Milestones)

## Summary
Focus on a staged roadmap to stabilize implicit IB, validate algorithmic correspondence, improve performance, and expand capability (parallel + kernel generalization).
Use deterministic CMake tests as the main confidence signal, but treat milestone checks as directional goals rather than strict gates.

## Roadmap Record
- Store and maintain roadmap at: `/Users/boyceg/code/IBAMR/doc/plans/implicit-ib/roadmap.md`
- Keep this file updated with:
1. Current reference configs
2. Milestones and status
3. Key comparison notes vs MATLAB
4. PR history and follow-up items

---

## Current State (2026-03-01)
1. Development is on branch `codex/implicit-ib-milestone-1`.
2. Milestone 1 work is in-progress and currently under review (not yet merged).
3. Reference implicit stabilization changes are applied in:
   - `/Users/boyceg/code/IBAMR/src/IB/IBImplicitStaggeredHierarchyIntegrator.cpp`
   - `/Users/boyceg/code/IBAMR/src/navier_stokes/INSStaggeredHierarchyIntegrator.cpp`
   - `/Users/boyceg/code/IBAMR/include/ibamr/INSStaggeredHierarchyIntegrator.h`
4. A promoted finished reference implicit test is added in:
   - `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_01.cpp`
   - `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_01.baseline.input`
   - `/Users/boyceg/code/IBAMR/tests/IB/implicit_stokes_ib_01.baseline.output`
   - `/Users/boyceg/code/IBAMR/tests/CMakeLists.txt` (test registration)
5. Focused validation currently passing:
   - `attest -R 'IB/implicit_stokes_ib_01.baseline.input'`

---

## Milestone 1 — Cleanup and stabilization
### Goal
Clean up implicit IB implementation and land a focused stabilization PR that restores reference implicit workflows and finished test coverage.

### Scope
1. Remove accidental/debug artifacts that interfere with reference implicit solves.
2. Stabilize the `IBImplicitHierarchyIntegrator` execution path and required INS-side plumbing.
3. Promote core reference implicit tests from unfinished/development status into finished attest-managed tests.
4. Keep this milestone strictly reference-oriented (no new smoother feature development).

### Status
- In progress.
- Completed so far:
1. Reference stabilization edits applied for `IBImplicitStaggeredHierarchyIntegrator` and `INSStaggeredHierarchyIntegrator` synchronization path.
2. Reference finished implicit attest test added and passing in focused runs.
3. Roadmap file created under `doc/`.
- Remaining before closeout:
1. Final review and cleanup of staged diffs.
2. Land commit sequence and changelog entry for milestone 1.

### Explicitly Deferred from M1
1. Coupling-aware vs geometrical ASM behavior work.
2. New smoother modes/policies and tuning keys.
3. Transfer-method feature expansion (RT0/FAC/PETSc tuning work).
4. MATLAB parity/performance studies for coupling-aware smoothers.

---

## Milestone 2 — Stokes infrastructure verification (without IB coupling)
### Goal
Verify that Stokes-only infrastructure corresponds to basic Vanka-style multigrid behavior, including known special cases where strict CAV and classical Vanka should coincide.

### Scope
1. Define reference Stokes-only cases.
2. Compare IBAMR behavior with MATLAB reference implementations.
3. Document where strict CAV should match classical Vanka and where behavior intentionally differs.

### Status
- Not started.

---

## Milestone 3 — RT0 prolongation validation
### Goal
Add confidence that RT0 prolongation is behaving as expected, especially divergence-related behavior for velocity transfer.

### Scope
1. Add focused tests/diagnostics for transfer invariants.
2. Check both FAC transfer settings and PETSc-side transfer settings.
3. Record expected behavior and known caveats in roadmap notes.
4. Verify the RT0 implementation behavior directly.
5. Add a deterministic transfer test for `R * Id * P`:
   - row sums are `1` for velocity DOFs and `0` for pressure DOFs,
   - column sums are `1` for velocity DOFs and `0` for pressure DOFs.

### Status
- Not started.

---

## Milestone 4 — MATLAB parity for Stokes + IB
### Goal
Build exact-comparison cases between IBAMR and MATLAB for Stokes+IB using programmatic structure generation, then identify and reduce discrepancies.

### Scope
1. Align geometry/parameters between codebases.
2. Compare subdomains, iteration traces, and residual behavior.
3. Track mismatches with likely causes and resolution plan.

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
- Deferred until Milestones 2 and 4.
- Current milestone does not treat MATLAB parity as a gate.

## PR History and Follow-up Items
- Add entries as PRs are opened/merged.
- Suggested entry format:
  - PR: `<number/link>`
  - Milestone: `<M#>`
  - Summary: `<what changed>`
  - Validation: `<tests/runs>`
  - Follow-up: `<next tasks>`
