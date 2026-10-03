# Focused, trustworthy tests

- Include a focused native `attest` regression in the same PR as each feature or
  bug fix. First look for a suitable executable on the actual target branch and
  add another case when initialization, lifecycle, and linkage fit. An executable
  on an unmerged sibling branch is not available coverage. Explain the need for a
  new executable in the task report and, when material to review, the PR. Do not
  add a separate CTest test in place of native coverage.
- Build the executable and use the normal fixture-linking targets. A CMake target
  alone, or a hand-copied test collection, is not integrated coverage. Verify
  `attest -N` discovers the intended input/output pairs before running them.
- Test the public operation using real library objects. An independent small
  reference calculation can check the answer, but a surrogate implementation
  cannot stand in for the IBAMR operation under test. Do not expose production
  internals or add general export/replay machinery just to simplify a test. Drive
  real setup and lifecycle operations instead of fabricating internal state. A
  subclass exercising a supported extension contract is appropriate.
- Keep cases deterministic, repository-local, network-free, and small. Use legal
  SAMRAI patch geometry and no host-specific output. New cases should each take
  less than a minute in Release; verify this budget during the deferred Release
  phase in [build-and-verification.md](build-and-verification.md). Cover distinct
  contracts and meaningful edge cases, not every imaginable internal state.
- Each case should cover a contract that no other case in the change covers. Do
  not repeat a case across dimension, process count, periodicity, grid size, or
  an option that the changed code does not read. An expected output identical to
  a sibling's is evidence that a case is redundant. Each expected-error case
  should reach a different check.
- A change that should not alter results, such as removing redundant work,
  normally needs no new case when existing cases reach the changed code and
  their expected outputs do not change. Do not add a case that counts internal
  operations such as ghost fills, messages, or timer calls. For changed branches
  that no committed case reaches, compare outputs with the base branch locally
  using temporary inputs, and report that comparison instead of committing the
  cases.
- Add an expected-error case when the error path contains logic that could
  plausibly be wrong, such as detecting a nonfinite or unconverged result. A
  precondition or argument check that is evident by inspection does not need one.
- For numerical results, normally write the values to `output` through `plog` and
  supply the matching expected output. `attest` uses `numdiff` to compare numbers
  with tolerances that allow small differences between processors, compilers,
  and dependency versions. Prefer this to replacing useful numerical values with
  Boolean `PASS` markers or duplicating the comparison in a custom test helper.
- Keep numerical output short, labeled, and readable by a reviewer. Select a few
  meaningful values or norms instead of dumping large vectors and matrices. Print
  enough digits to expose relevant differences. Do not write the input database
  to `output`, so that editing an input's comments does not change its expected
  output. Use an additional numerical check
  when a required property is not adequately checked by the output comparison.
- Use existing assertions for suitable non-numerical conditions instead of a
  bespoke Boolean-checking wrapper. Ensure required checks remain active in the
  tested configuration: a Debug-only assertion is not Release coverage. An
  assertion-only case can have empty expected output only when its actual logging
  is empty; the run must still create `output`. With `AppInitializer`, normally
  set `log_file_name = "output"` in the input file's `Main` database. In a test
  without `AppInitializer`, follow the existing `std::ofstream` output pattern.
  Do not remove real integrator logging to obtain an empty baseline or add a
  transcript of `PASS` messages.
- An expected-error case must compare the actual, stable diagnostic against a
  nonempty expected output. Reuse logging helpers such as `TestAppender` in
  [tests/tests.h](../../tests/tests.h) where suitable. If the operation unexpectedly
  succeeds, return normally so `expect_error=true` rejects the run.
- Keep test-only instrumentation in tests, not the library. For example, an
  allocation counter may use PETSc logging when enabled, but must not report zero
  allocations as measured when logging is unavailable; numerical checks still run.

For example, `plog << "interpolated_velocity = " << velocity << '\n';` leaves a
useful numerical result for `numdiff` and a human reviewer. Printing every entry
of an interpolation matrix usually does not. Choose output precision appropriate
to the quantities being compared, and investigate changed results before updating
expected-output files or tolerances.

For an expected-error case, call the real operation and let unexpected
continuation return `EXIT_SUCCESS`. Do **not** print the anticipated diagnostic
yourself and then abort with "the operation should have failed": that can pass
even when the intended validation is absent. Check a new failure mechanism once
with a small, recoverable negative control, then restore the fixture. For a bug
fix, normally demonstrate that the regression fails without the fix and passes
with it; report both results, or explicitly explain why that comparison is
impractical. A before/after comparison that exercises the failure mechanism is
sufficient; do not duplicate it with another mutation check. Use an isolated,
recoverable comparison that preserves other work. Do not require mutation testing
of every new test or introduce a permanent mutation-testing framework.

Example, with `ibamr_source` and `ibamr_build` set to the intended source and
configured build paths:

```sh
cmake --build "$ibamr_build" --target tests-IBTK
cd "$ibamr_build"
"$ibamr_source/attest" -N -R '^IBTK/petsc_options_file\.'
"$ibamr_source/attest" -R '^IBTK/petsc_options_file\.'
```

Run from the build root containing `attest.conf` and the normal `tests/` links;
the runner's source location does not determine discovery. Select the collection
appropriate to the change. Use joined job syntax such as `-j2` if needed.

## Fixture discovery and failure checks

- Use descriptive `name.case[.mpirun=N][.expect_error=true].input` and matching
  `.output` names. The prefix before the first dot selects the executable.
  `SETUP_2D` and `SETUP_3D` add `_2d` and `_3d` to executable names; ensure the
  target and linkage actually exist for each fixture's dimension.
- Build the collection target, such as `tests-IB` or `tests-IBTK`, to run its
  fixture-linking step. Building only `tests-IB_name` does not run that step.
  The manual link command takes the two test directories, not source/build roots:
  `bash "$ibamr_source/tests/link-test-files.sh" "$ibamr_source/tests/IB" "$ibamr_build/tests/IB"`.
  Prefer the collection target. Before using the manual command, inspect the
  target branch's `tests/CMakeLists.txt` for additional setup or exclusions that
  the command would bypass. Confirm the intended cases appear in `attest -N`;
  investigate missing cases rather than assuming they were deliberately excluded.
  A successful discovery command alone does not prove that cases were found.
- `attest` compares the actual `output` file even for an expected-error case.
  Generate new expected output from a real run and inspect it. Preserve required
  diagnostic whitespace; explain any necessary, narrowly scoped exception to
  `git diff --check` instead of exempting all expected-error files.
- In IBSAMRAI2, `TBOX_ASSERT` is unconditional; `TBOX_CHECK_ASSERT` depends on
  `DEBUG_CHECK_ASSERTIONS`. Check the installed definitions and surrounding
  preprocessor guards when relying on Release coverage. Use a condition with
  `TBOX_ERROR` when an explicit runtime failure is needed.
- Ensure failed conditions actually fail the test. Reject nonfinite values when
  a norm or maximum could hide a NaN. Keep checks proportional to the contract.

For hierarchy setup and norm comparisons, follow the ghost-width registration
and composite-norm masking guidance in [Errors and numerics](errors-and-numerics.md).

## Input file comments

- Begin a new input file with at most one or two comment lines that say what the
  case sets up, and comment the setting that distinguishes it from its siblings.
  Many existing input files have few comments; add only what a reader of that
  file needs.
- Do not name library functions, describe library internals, narrate the defect
  or the change, or explain at length what the executable checks. Explain the
  checks in the test source.
- When copying an input file, update or remove its comments, and remove keys the
  new case does not read and settings that are commented out. Use the same
  wording for the same thing in sibling input files.

## Expected-output changes

Never replace an expected result or loosen a tolerance merely to obtain a pass.
Trace each numerical baseline change to an intended algorithm or behavior change
and confirm it with an independent calculation or control run. Review new
baselines with the same care. List changed baselines and reasons in the task
report, grouping equivalent changes; distinguish numerical changes from diagnostic
or other text-only changes. Reviewed generation and deliberate text edits are
allowed when they reflect verified behavior.
