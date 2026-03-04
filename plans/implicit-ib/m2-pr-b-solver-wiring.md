# Milestone 2 PR-B: Solver Wiring and Mode Selection

## Scope
1. Wire ASM subdomain construction mode selection into `StaggeredStokesPETScLevelSolver`.
2. Select exactly one path per solve setup:
   - `GEOMETRICAL`
   - `COUPLING_AWARE`
3. Keep CAV subdomain construction A00-driven via the solver-owned level matrix.
4. Add deterministic navier_stokes integration coverage for solver-mode ASM subdomain generation.

## Implemented Changes
1. Added `ASMSubdomainConstructionMode` enum and `string_to_enum` / `enum_to_string` conversions in:
   - `include/ibamr/ibamr_enums.h`
2. Extended `StaggeredStokesPETScLevelSolver` state and input parsing for:
   - `asm_subdomain_construction_mode`
   - `coupling_aware_asm_seed_axis`
   - `coupling_aware_asm_seed_stride`
   - `coupling_aware_asm_closure_policy`
   - `log_ASM_subdomains`
3. Replaced unconditional geometric subdomain construction in `generateASMSubdomains()` with a `switch` on mode.
4. In coupling-aware mode:
   - extract velocity submatrix from solver-owned `d_petsc_mat` with `constructA00VelocitySubmatrix()`,
   - build patch-level closure map data,
   - construct coupling-aware overlap/nonoverlap subdomains.
5. Added deterministic integration test:
   - `tests/navier_stokes/stokes_petsc_level_solver_asm_modes.cpp`

## Validation
1. `cmake --build <build_dir> -j10 --target tests-navier_stokes`
2. `<build_dir>/attest -j8 --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R "navier_stokes/stokes_(ib_cav|petsc_)|IB/implicit_stokes_ib_01.baseline.input"`

## Deferred to Later Milestones
1. Parallel validation and MPI-coverage invariants for CAV/ASM.
2. MSM/GS-like application-semantics tuning and performance optimization.
3. MATLAB parity studies that depend on smoother application semantics.
