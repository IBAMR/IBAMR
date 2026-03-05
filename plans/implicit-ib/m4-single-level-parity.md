# Milestone 4 Phase 1: Single-Level MATLAB/IBAMR Structural Parity

## Scope
Phase 1 adds a deterministic single-level parity gate for the Stokes+IB operator chain in serial (`np=1`) and verifies:
1. Lagrangian positions `X`.
2. Lagrangian force Jacobian `A` and force vector `F = A * X`.
3. Interpolation/spreading operators `J` and `S`.
4. Composite `SAJ` and velocity block embedding `A00`.

This phase is test/reference only and does not introduce production solver API changes.

## MATLAB-to-IBAMR Mapping
1. MATLAB structure setup (`alpha`, `beta`, `r_cyl`, `ds` rule, point ordering)
   - IBAMR: `IBRedundantInitializer` callbacks in `implicit_stokes_ib_m4_single_level_parity_01.cpp`.
2. MATLAB elastic operator and force
   - IBAMR: `IBMethod::constructLagrangianForceJacobian()` and `MatMult(A, X, F)`.
3. MATLAB interpolation/spreading operators
   - IBAMR: `IBMethod::constructInterpOp()` for `J`, and `MatTranspose(J, ...)` for `S`.
4. MATLAB `SAJ` and Eulerian embedding
   - IBAMR: `MatPtAP(A, J, ...)` for `SAJ`, then `constructPatchLevelMACStokesOp()` + `constructA00VelocitySubmatrix()` for `A00`.

## External-Reference Artifacts
Generated locally (developer workflow) under:
`<build_dir>/tests/external-reference/cav-stokes-ib/reference-generated/single_level_case_01/`

The following files are generated and compared in external-reference runs:
1. `X.ref`
2. `A.ref`
3. `F.ref`
4. `S.ref`
5. `J.ref`
6. `SAJ.ref`
7. `A00.ref`
8. `manifest.json`

Format:
1. Matrix: header `nrows ncols nnz`, then sorted `row col value` triplets.
2. Vector: header `n`, then `idx value`.

Comparison policy:
1. Index/shape/nnz equality must be exact.
2. Numeric values must satisfy `|err| <= 1e-12`.

## Current Status
1. New deterministic test added: `tests/IB/implicit_stokes_ib_m4_single_level_parity_01.cpp`.
2. CMake test registration added.
3. CMake attest self-contained runs pass for:
   - `IB/implicit_stokes_ib_m4_single_level_parity_01.multilevel_l2.input`
   - `IB/implicit_stokes_ib_m4_single_level_parity_01.multilevel_l3.input`
4. Multilevel parity support is implemented and passing for deterministic 2-level and 3-level cases:
   - `IB/implicit_stokes_ib_m4_single_level_parity_01.multilevel_l2.input`
   - `IB/implicit_stokes_ib_m4_single_level_parity_01.multilevel_l3.input`
5. Matched-case CAV subdomain parity support is implemented and passing for deterministic single-level checks:
   - `IB/implicit_stokes_ib_m4_cav_subdomain_parity_01.input`
6. Level-wise references are generated and compared for velocity-block operators in external-reference workflows:
   - `SAJ.level*.ref`
   - `A00.level*.ref`
7. Level-wise `SAJ` parity is transfer-consistent with the MATLAB/Octave baseline by using RT0 transfer semantics (`R * SAJ_u * P`) for coarse levels.
8. Reference generation now uses IB4 consistently in the baseline path.
9. Committed `.ref` datasets are intentionally excluded from CI-oriented test coverage.

## Test Classification Audit (2026-03-04)
1. CI-only self-contained checks:
   - `tests/IB/implicit_stokes_ib_m4_single_level_parity_01.multilevel_l2.input`
   - `tests/IB/implicit_stokes_ib_m4_single_level_parity_01.multilevel_l3.input`
   - `tests/IB/implicit_stokes_ib_m4_cav_subdomain_parity_01.input`
   - `tests/navier_stokes/stokes_ib_cav_relaxed_semantics_parity_{2d,3d}.input`
   - `tests/navier_stokes/stokes_ib_cav_strict_semantics_parity_{2d,3d}.input`
2. Developer-only external-reference checks:
   - `tests/external-reference/cav-stokes-ib/implicit_stokes_ib_reference_single_level_01.input`
   - run through CMake target `run-external-reference-cav-stokes-ib`
   - require local Octave/MATLAB baseline tooling and generated reference files.
3. Naming policy status:
   - `matlab` appears only in external-reference tooling/docs, not in self-contained CI test names.
4. Data policy status:
   - committed reference `.ref` files are not used; references are generated locally under build-tree external-reference paths.

## Milestone 4 Closeout
1. Milestone 4 structural parity goals are complete on this branch for the planned single-level and multilevel matched cases.
2. CI coverage remains self-contained while external-reference parity remains available as a developer-only workflow.
3. Iteration-trace/performance equivalence is deferred to Milestone 5 (smoother application semantics).
