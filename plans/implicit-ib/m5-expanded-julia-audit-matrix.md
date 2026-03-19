# Milestone 5: Expanded Julia-First Audit Matrix

Date: 2026-03-16
Scope: `IB-implicit-ex0`, 2D, serial, shell CAV parity against the canonical Julia reference.

## Goal

Expand the current one-case smoke audit into a small but explicit matrix of Julia-reference cases, update the audit contract around those cases, and only then repeat the same matrix in IBAMR with an exactly matched solver contract.

This ordering is deliberate:
- first freeze the Julia/MATLAB-side contract over a broader case set,
- then document the exact expected behavior and artifacts,
- then compare IBAMR only on runs that are known to satisfy that contract.

## Completion Checklist

Status legend:
- `done`: completed and documented
- `in_progress`: active work, partially verified
- `blocked`: waiting on a bug fix or prerequisite
- `todo`: not started yet

### Stage A: Locked Two-Level `16 -> 8` Audit

- `done` Freeze the exact two-level Julia smoke contract in the canonical reference docs.
- `done` Freeze the matching IBAMR `ex0` contract for the same case.
- `done` Match outer Krylov contract.
- `done` Match inner FAC contract.
- `done` Match force convention for audit exports.
- `done` Match full-space transfer exports for the audit path.
- `done` Match fine/coarse Eulerian operators.
- `done` Match subdomain index sets.
- `done` Match local `A00` blocks.
- `done` Match local full-saddle blocks.
- `done` Verify Stage A stiffness cases `K = 1e0` through `1e7` on lightweight parity objects.
- `in_progress` Verify full `outer_sol` parity for Stage A high-stiffness cases `K = 1e5, 1e6, 1e7`.

### Stage B: Deep Hierarchy Audit

- `done` Decide the next canonical deep IBAMR audit case.
- `done` Document that `16 -> 8 -> 4` is the current next IBAMR-auditable deep case.
- `done` Document that `16 -> 8 -> 4 -> 2` is currently unsupported in `ex0` because of SAMRAI patch-size restrictions.
- `done` Confirm by experiment that the `2x2` coarsest hierarchy is rejected in both DEBUG and OPT.
- `done` Confirm by experiment that the `4x4` fallback hierarchy runs in OPT.
- `done` Fix the DEBUG-mode `16 -> 8 -> 4` crash.
- `done` Freeze the exact Stage B contract in the canonical reference docs.
- `done` Freeze the exact matching Stage B IBAMR contract in `ex0`.
- `in_progress` Re-run the full operator, transfer, force, subdomain, and outer-solve audit on the deep Stage B case across representative stiffness values.

### Stage C: Broader Grid Matrix

- `todo` Choose the Stage C grid set to certify on the Julia/MATLAB side.
- `todo` Freeze those Stage C cases in the canonical reference docs.
- `todo` Replicate those exact Stage C cases in `ex0`.
- `todo` Complete parity checks for operators, transfers, forces, subdomains, and outer solves on Stage C.

### Audit Exit Criteria

The audit is complete when all of the following are true:
- every Stage A checklist item is `done`
- every Stage B checklist item is `done`
- the intended Stage C scope has been explicitly chosen, and every chosen Stage C item is `done`
- the canonical reference docs and the IBAMR-local tracker agree on the certified case matrix
- every certified IBAMR comparison was run with a solver contract confirmed by `KSPView/PCView`

## Current Locked Baseline

The current audited baseline is:
- example: `examples/IB/implicit/ex0`
- geometry: Julia smoke geometry
- hierarchy: two levels, finest `16x16`, coarse `8x8`
- outer solver:
  - `gmres`
  - left preconditioning
  - preconditioned norm
  - `rtol = 1e-8`
  - `max_it = 150`
  - `restart = 150`
- inner FAC contract:
  - shell `coupling_aware_vanka_matlab`
  - local patch solve `preonly + svd`
  - coarse solve `preonly + svd`
  - `num_pre_sweeps = 1`
  - `num_post_sweeps = 1`
  - relaxation factor `1.0`
- audit-path conventions:
  - force export in Julia force-density units
  - full-space `P/R` export matched to Julia
  - pressure-row sign normalized in the comparator
  - deterministic numerical-zero filtering for subdomain growth and local block export

For this locked baseline, the audit now matches Julia on:
- marker positions
- marker forces
- fine/coarse Eulerian operators
- full-space transfer operators
- outer RHS consistency
- subdomain index sets
- local `A00` blocks
- local full-saddle blocks

The lightweight stiffness sweep for this baseline also matches through `K = 1e7`.

## Expanded Julia-First Audit Matrix

The audit should expand in stages so each new stage adds one major source of variation.

### Stage A: Same hierarchy, broader stiffness range

Purpose:
- verify that the current matched `16 -> 8` contract remains stable in the high-stiffness regime where discrepancies were previously reported.

Cases:
- finest `16x16`, depth `2`, coarse `8x8`
- `K = 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7`

Required checks:
- outer iteration count
- `marker_forces`
- `L.level0`, `L.level1`
- `outer_rhs_vs_Lseed`
- subdomain set parity
- local `A00` parity
- local full-saddle block parity
- full `outer_sol` parity for at least `K = 1e5, 1e6, 1e7`

Current status:
- completed for all listed `K` except the heavy full `outer_sol` comparison at `1e5+`
- no structural drift has been observed through `K = 1e7`

### Stage B: Same finest grid, deeper hierarchy

Purpose:
- determine whether parity holds when the coarse path becomes materially more important.

Target cases:
- canonical next IBAMR-auditable deep case:
  - finest `16x16`, depth `3`, coarsest `4x4`
- deeper canonical reference target retained on the Julia/MATLAB side:
  - finest `16x16`, depth `4`, coarsest `2x2`

Current IBAMR status:
- the `16 -> 8 -> 4 -> 2` `ex0` case is currently unsupported because SAMRAI rejects the `2x2` coarsest hierarchy with minimum patch size constraints
- the `16 -> 8 -> 4` case is therefore the canonical next deep audit case for IBAMR
- the earlier DEBUG-mode `SIGSEGV` on the `16 -> 8 -> 4` shell-CAV path is now fixed
- representative Stage B runs now complete in both DEBUG and OPT for `K = 1e2, 1e4, 1e6`
- lightweight deep parity objects remain matched through `K = 1e7` on:
  - marker forces
  - levelwise Eulerian operators
  - outer RHS consistency
  - per-level subdomain index sets
  - per-level local `A00` blocks
  - per-level local full-saddle blocks

Required checks:
- everything from Stage A
- explicit `P/R` and coarse-level operator parity on every level
- coarse solve contract confirmation via `KSPView/PCView`
- per-level subdomain parity on every relaxed level

Acceptance rule:
- no IBAMR comparison for a Stage B case is meaningful unless the Julia/MATLAB-side contract for that exact hierarchy has been documented first

Current representative Stage B iteration counts on the locked `16 -> 8 -> 4` audit path:
- `K = 1e2 -> 5`
- `K = 1e4 -> 9`
- `K = 1e6 -> 12`
- `K = 1e7 -> 17`

### Stage C: Broader finest-grid configurations

Purpose:
- check whether parity is specific to the `16 -> 8` smoke hierarchy or persists under grid refinement/coarsening changes.

Candidate cases:
- finest `8x8`, depth `2`
- finest `16x16`, depth `2`
- finest `32x32`, depth `2`
- possibly finest `32x32`, depth `3` if the Julia reference remains practical

Required checks:
- same contract items as Stage A
- explicit DOF-map-based comparison only; raw index ordering is not admissible

## Contract Update Checklist

Before each stage is used for IBAMR comparison, the Julia/MATLAB-side contract should be updated to record:
- exact geometry construction rule
- exact marker count and `ds` convention
- exact stiffness convention
- hierarchy depth and coarsest grid
- outer Krylov contract
- FAC/smoother contract
- transfer contract
- local and coarse solve contract
- force and pressure sign conventions
- expected audit artifacts and file names

The canonical home for that update remains the standalone solver repository:
- [matlab-reference-audit.md](/Users/boyceg/code/implicit_ib_coupling_aware_vanka/docs/matlab-reference-audit.md)
- [matlab-reference-algorithm-contract.md](/Users/boyceg/code/implicit_ib_coupling_aware_vanka/docs/matlab-reference-algorithm-contract.md)

This IBAMR-local document is the staging checklist for what needs to be reflected there before the next comparison phase is considered authoritative.

## IBAMR Audit Handoff Rule

For any new case, IBAMR comparison is only considered meaningful if all of the following are true:
- the case exists in the Julia/MATLAB audit matrix
- the case's contract is documented explicitly
- `ex0` is configured to match that contract exactly
- `KSPView/PCView` confirms the intended live solver contract
- the DOF-map-aware export bundle is complete

If any of those are false, the result should be labeled as exploratory, not parity-audited.

## Immediate Next Steps

1. Finish the Stage A full `outer_sol` checks for `K = 1e5, 1e6, 1e7`.
2. Finish the Stage B deep audit sweep on representative high-stiffness cases.
3. Add full deep `outer_sol` checks for the highest-priority Stage B stiffnesses.
4. Decide whether the next certified expansion should be deeper hierarchy coverage or broader Stage C grid coverage.
