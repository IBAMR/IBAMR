# Matrix-free experiment on PR #1997

2026-09-25. The experiment is rebased onto PR #1997,
`e80aa6cae02f054a733a65d64b9c016bf0b89b70`, from GitHub Stack #2026.
This base provides the current kernel evaluator, dispatch and matrix-assembly
APIs. The experiment requires IBTK and does not depend on the later implicit
Jacobian, Stokes-IB or FAC changes.

## Preserved history and scope

The eight experiment commits after qualified base
`c37b111fdabf267e83229c069698eca1a152bf1e` were replayed onto #1997.
The pending independent-review documentation was committed separately before
the replay and carried forward. Recovery refs preserve both states:

- `codex/ibkernels-matrix-free-pre-rebase-20260925`:
  original measured/report head `d4287f92a8aee00f6e467b5aaae9a8343d9760cc`.
- `codex/ibkernels-matrix-free-review-snapshot-20260925`:
  the same source with the previously uncommitted review documentation.

The original initial/qualified-base refs and existing evidence remain in place.
Snapshots, file hashes, replay comparison and validation logs are under
`evidence/experiments/ibkernels-matrix-free-20260915/rebase-1997-20260925/`.
All seven historical raw-sample files and eight preserved binaries/libraries
match their original recorded hashes.
No upstream stack branch or GitHub PR was changed or published.

## API adaptation

The tensor-product header conflicts were resolved by retaining #1997's
constructors, stencil concepts, weight restrictions and evaluation code, then
adding the experiment's owning `evaluateFactors()` interface. Its constraints
now name `NormalEvaluator` and `TransverseEvaluator`.

Experimental callers use `ib_kernel_evaluators.h`, the `IBKernelEvaluators`
namespace, and `ib_kernel_weights_value_t` / `ib_kernel_weights_extent_v`.
The current base requires array weight storage; the experiment already uses
arrays. The supported scalar kernel formulas are mathematically unchanged.
No coupling loop, traversal order, clipping rule, scaling or accumulation
policy was changed during the port.

The 3D Fortran IB5 x-index correction remains included. It is absent from #1997
and #1998; the same correction appears later in Stack #2026 in PR #1965.

The Silo-enabled Debug build exposed inherited unused `SILO_NAME_BUFSIZE`
constants in `LSiloDataWriter.cpp` and `IBInstrumentPanel.cpp`: the base had
already removed the buffers that used them. Both declarations and their
no-Silo `NULL_USE` markers were removed locally to retain the strict warning
settings. These four deleted lines do not change Silo output or numerical
behavior.

## Validation and measurement boundary

Correctness validation uses a separate persistent Debug build at
`build/ibkernels-matrix-free-pr1997/Debug`, with the Apple developer toolchain,
existing Debug dependencies, shared ccache and the applicable macOS CI warning
flags. The original `build/ibkernels-matrix-free/` artifacts are preserved.
`configure.sh` now defaults to the new build root and accepts
`EXPERIMENT_BUILD_ROOT` as an override.

The focused native selection covers 2D/3D matrix-free coupling, kernel
evaluators, operator builders and interpolation-matrix assembly, including the
periodic matrix case. Both benchmark executables are compiled for API coverage.
Build, discovery and test results are recorded in the evidence directory.

- Strict Debug build: passed for `tests-matrix_free`, `tests-IBTK` and both
  benchmark executables, using Apple Clang 17 and GNU Fortran 15.2.
- Native discovery: 11 tests. Initial focused run: **11/11 passed**. Repeat:
  **11/11 passed**, with retained work directories and output comparisons.
- Matrix-free coverage includes BS2–6, all ten CBS(k+1)k/CBSk(k+1) combinations
  for k=1–5, IB4/IB5, factor ownership, all three application modes, generic
  fallback, clipping, indexing, shifts and the adjoint check. The existing
  tolerances and expected outputs were retained.
- The selection includes two-rank evaluator and matrix-assembly tests, plus
  the expected-error evaluator cases.
- Required clang-format 16 formatting and `git diff --check`: passed.

The linker reported duplicate-library/rpath and dependency deployment-version
warnings: some PETSc dependencies target macOS 15.7 while this build targets
15.0. They did not prevent linking or the focused tests. No warning suppression
or deployment-target changes were added for this migration.

Independent source review found no port defect or documentation-provenance
error. Separate strict compiler syntax checks passed with Silo disabled for
both affected source files in 2D and 3D; these are not a full no-Silo build.

Release validation, benchmarks and profiling remain deferred. Historical
performance reports retain their original measured revisions and do not
establish performance of the rebased implementation. The deferred cache and
CBS65 anomaly investigations remain unchanged.
