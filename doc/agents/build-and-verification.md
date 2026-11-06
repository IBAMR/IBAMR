# Build efficiently without weakening checks

- Use Debug compilation and tests for routine implementation, review maintenance,
  and restacking verification. Complete the assigned dimensions, native test
  discovery, focused tests and reruns, and any assigned broad suite in Debug.
  Defer Release builds, tests, and runtime-budget measurements to the final
  working-stack integration and performance phase, or an explicit targeted request.
  Pending Release qualification does not block implementation-task completion,
  handoff, or otherwise-authorized publication.
- Reuse compatible out-of-tree builds when available; configure a new build when
  needed. Verify source paths, compilers, dependencies, and effective compile
  commands. "Fresh results" means results from the current source. When persistent
  storage is available, retain useful builds and caches across related work. Keep
  verification results for later qualification and handoff.
- When `ccache` is available and suitable for the toolchain, consider using it to
  speed up repeated builds. Configure the relevant compiler launchers
  (`CMAKE_C_COMPILER_LAUNCHER`, `CMAKE_CXX_COMPILER_LAUNCHER`) and check activity
  with `ccache -p` and `ccache -s -v` after a build that compiles files. Use valid
  default cache locations, or set `CCACHE_DIR` and `CCACHE_TEMPDIR` to writable
  locations suited to the environment. For cross-worktree reuse, investigate
  debug-path effects and validate `CCACHE_BASEDIR` and any compiler-specific
  debug-prefix mapping together. Preserve useful cache contents and correctness
  settings when tuning reuse.
- Match applicable CI warning/error settings in both library and test builds,
  in each configuration when run. Read the current workflow rather than freezing
  its flags here. Verify effective flags, not just requested CMake values.
  Compiler-specific flags and dependency exceptions need equivalent intent,
  not blind copying between Clang, GCC, and Fortran.
- The hosted jobs build with GCC on Linux and with Clang on macOS and treat
  warnings as errors, and each compiler reports problems that the other does
  not. Where a compiler of the family that the local build does not use is
  installed, consider compiling the changed C++ sources with it before
  publishing:
  `scripts/maintenance/check_pull_request.py --compiler-check <build directory>`
  does a syntax-only compilation with the warning flags of the hosted builds,
  using GCC for a Clang build and Clang for a GCC build. This is a
  recommendation, not a requirement: the tools are not available on every
  system, and the hosted jobs remain the check that decides.
- Fix warnings introduced by the change. Report inherited warnings separately;
  agree on a narrow, visible exception when necessary instead of silently dropping
  `-Werror` or turning a feature PR into whole-tree cleanup. Do not enable
  `-Weverything` as a substitute for selecting relevant diagnostics.
- Serialize edits, configuration, formatting, and builds using the same build
  directory. Stop an owned build promptly if a required configuration change
  invalidates it; preserve unrelated work. Rebuild dependent executables after
  an ABI/header change before diagnosing apparent runtime regressions.
- Run the affected tests on the final formatted code and repeat sensitive cases
  to check determinism. Include relevant setup/teardown and integration coverage
  for ownership or solver-composition changes. Run broader suites when warranted;
  focused passes do not replace an assigned broad run such as `attest -R 'IB/'`.
  Do not run the whole test suite locally as a routine check: run the test
  groups the change touches, and by name any test that CI excludes and the
  change plausibly affects, and rely on CI for the rest.
- Report commands, configuration, tested revision, and limitations concisely.
  A 3D library build is not a 3D test run; a simulated compatibility compile is not
  testing an old dependency installation. Local passes are not hosted CI passes.
  Report deferred Release qualification as pending, not as passed based on Debug.
  Record test duration for practicality, not as solver-performance evidence.

## Configure before expensive builds

Read [doc/cmake.md](../../doc/cmake.md), the current CMake options, and the
[CI workflow](../../.github/workflows/push-pull.yml).
Use a compiler toolchain that meets the project's language and dependency
requirements. Choose diagnostics and profiling tools suited to the environment
and the question being investigated. Use the CI matrix to identify relevant
compatibility checks, adapting commands and paths to the local environment.

Use compatible dependency installations. Installed headers and version-matched
documentation can establish interface behavior; when dependency source is
available, use it to investigate implementation details relevant to the change.
Choose evidence appropriate to the claim and report verification limits.
Locate the actual dependency roots and library variants, and verify paths and
options from any existing cache before reusing them. Choose compiler and warning
flags up front. Changing flags changes cache keys, so avoid unnecessary
full rebuilds while retaining strict checks. Where a change touches dependency
APIs, compare local versions with the CI matrix and account for API differences.
Do not turn an inherited warning from one machine into a standing suppression.

## Formatting in worktrees

- Before finalizing source commits, run the repository formatter through the
  configured build's `indent` target, for example `make indent` for Make builds
  or `cmake --build "$ibamr_build" --target indent` for CMake builds. The source
  entry point below is also available. Use the clang-format version pinned by
  IBAMR and required by CI. The authoritative version check is in
  [scripts/formatting/indent_common.sh](../../scripts/formatting/indent_common.sh).
  IBAMR provides `scripts/formatting/download-clang-format` and
  `scripts/formatting/compile-clang-format` to obtain a compatible formatter when
  needed. See [scripts/formatting/README.md](../../scripts/formatting/README.md); the
  source-tree formatting entry point is `scripts/formatting/indent`.
  Inspect all formatter changes and keep unrelated edits out. Markdown-only
  changes need document checks, not a library rebuild or a claim that source
  formatting was run.

Ignored formatter binaries under `scripts/formatting/programs/` are not copied
into new worktrees. Locate a compatible existing binary and put it on `PATH`, or
use the download/build helpers. Run the source entry point from the worktree root
so it uses that checkout's configuration. Check the required version in
`indent_common.sh` rather than assuming the newest installed formatter is suitable.
The changed-file formatter also includes untracked files.

`scripts/formatting/check-indentation` runs `indent-all` and can modify files;
use it in a controlled checkout and inspect the changes. A compatible
`clang-format --dry-run --Werror` on changed C++ files can supplement review but
does not replace the repository formatting pipeline.
