# Canonical CAV Recovery Review Findings Ledger

**Canonicalization date:** 2026-08-18

**Original handoff source SHA-256:** `2eb0d1b842a503c7850946d88f17bebf2ee822f6feceb0ce9b7ac39b5258a6b6`

**Status:** authoritative recovery input; preserve every disposition and design intent below

This repository file canonically incorporates the historical handoff ledger that was originally supplied as an external attachment. Independent audit tasks must use this file from the exact candidate commit and must not depend on the attachment path. The attachment is provenance only; this file is the durable audit input.

## Finding 1: Add tests in branch 1.1 support APIs

- **Code area:** `ibtk/include/ibtk/PETScMatUtilities.h:139` and related branch-1.1 support APIs.
- **Original concern:** branch 1.1 adds public `PETScMatUtilities` delta/interpolation APIs, `IndexUtilities` support, `IBImplicitStrategy`/`IBMethod` hooks, and `StaggeredStokesPETScLevelSolver` behavior, but does not add narrow tests in the same branch.
- **Disposition:** intentionally deferred for this cleanup pass. Do not treat this as accidentally missed.
- **Implementation:** none.
- **Future action:** only add these tests if explicitly requested. If reopened, tests should land in the branch where the support APIs are introduced.

## Finding 2: Factor shared Stokes-IB operator plumbing

- **Code area:** `src/IB/StaggeredStokesIBOperator.cpp` and `src/IB/StaggeredStokesIBJacobianOperator.cpp`.
- **Original concern:** nonlinear and Jacobian operator paths duplicated Stokes operator time setup, velocity-state selection, Stokes apply setup, and boundary-condition forwarding.
- **Disposition:** implemented.
- **Implementation:**
  - `3fded45fc` — Factor Stokes-IB operator plumbing.
  - `d4a736372` — Reuse Stokes-IB operator plumbing in Jacobian.
- **Design intent:** keep public APIs unchanged; share only private/common plumbing so time and boundary semantics remain identical across nonlinear and Jacobian paths.

## Finding 3: Collapse duplicate FAC residual vector handling

- **Code area:** `src/IB/StaggeredStokesIBLevelRelaxationFACOperator.cpp`, `computeResidual()`.
- **Original concern:** rediscretized and non-rediscretized residual paths duplicated PETSc `Vec` creation, SAMRAI/PETSc copies, copyback, synchronization, ghost fill, and destruction.
- **Disposition:** implemented.
- **Implementation:** `d65dd80fe` — Factor Stokes-IB FAC residual vector handling.
- **Design intent:** centralize vector lifetime/copy handling while preserving the separate rediscretized and non-rediscretized residual formulas.

## Finding 4: Reuse transfer-entry construction instead of two full passes

- **Code area:** `ibtk/src/math/PETScMatUtilities.cpp`, linear cell and side prolongation builders.
- **Original concern:** PETSc preallocation and value-fill phases duplicated stencil/weight construction, side geometry, offsets, ownership checks, and `MatSetValues` preparation.
- **Disposition:** implemented.
- **Implementation:** `52dd7d9e4` — Share linear prolongation row construction.
- **Design intent:** share per-row candidate column/value construction between preallocation and insertion. Retain separate PETSc preallocation and value-fill phases. Do not introduce row-entry caching unless there is a separate, explicit reason.

## Finding 5: Deduplicate velocity-block embedding

- **Code area:** `src/navier_stokes/StaggeredStokesPETScLevelSolver.cpp`, `initializeSolverStateSpecialized()`.
- **Original concern:** `MATPREALLOCATOR` and final AIJ matrix paths duplicated `MatGetRow`, `AOApplicationToPetsc` column mapping, row mapping, `MatSetValues`, and `MatRestoreRow`.
- **Disposition:** implemented.
- **Implementation:** `efe5bc430` — Deduplicate velocity-block embedding.
- **Design intent:** use one local helper to copy velocity-block rows into a destination matrix; keep matrix assembly at the call sites.

## Finding 6: Hide mutable coupling-aware map internals

- **Code area:** `include/ibamr/StaggeredStokesPETScMatUtilities.h`, `PatchLevelCellClosureMapData`.
- **Original concern:** public maps, source indices, and built flags could be mutated independently, making cache correctness hard to reason about.
- **Disposition:** implemented.
- **Implementation:** `2e77def20` — Hide coupling-aware map cache internals.
- **Design intent:** preserve `PatchLevelCellClosureMapData` as the public cache type and keep `clear()` public, but make internals private. Tests that need introspection should use narrow private-access helpers.

## Finding 7: Share Eigen shell smoother sweep logic

- **Code area:** `ibtk/include/ibtk/private/PETScLevelSolverEigenLocalSolveShellBackend-inl.h` and the Stokes Schur shell backend.
- **Original concern:** generic Eigen local-solve and Stokes Schur backends duplicated additive/multiplicative traversal, correction application, residual update, and postprocessing logic.
- **Disposition:** implemented.
- **Implementation:**
  - `b222a691f` — Share Eigen local-solve sweep logic.
  - `75c612f55` — Reuse Eigen shell sweep logic in Schur backend.
- **Design intent:** shared Eigen shell infrastructure owns traversal/update mechanics; the Stokes Schur backend keeps only Schur-specific block algebra and local solve behavior.

## Finding 8: Consolidate transfer test helpers

- **Code area:** `tests/navier_stokes/stokes_petsc_mat_utilities_pressure_transfer.cpp` and side-transfer helper code.
- **Original concern:** pressure-transfer tests duplicated profile parsing and max-norm/nontrivial checks already present in side-transfer helper logic.
- **Disposition:** implemented.
- **Implementation:** `70f2fae95` — Share transfer test helper checks.
- **Design intent:** share parsing and common validation helpers while keeping cell-profile and side-profile filling separate where geometry differs.

## Recovery rule

The implementation task may preserve or strengthen these dispositions, but it must not silently reopen an intentionally deferred item, reintroduce resolved duplication or mutable internals, or change the design intent. Any proposed change to a disposition requires an explicit dated decision in `cav-implementation-plan.md` and independent review at the affected candidate SHA.
