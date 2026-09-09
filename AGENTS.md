# IBAMR development guidance

IBAMR provides immersed boundary methods and related numerical methods on
adaptively refined Cartesian grids. Its implementation combines C++ library
code with Fortran numerical kernels and uses SAMRAI, PETSc, and other scientific
libraries.

Changes should preserve numerical correctness and efficiency while keeping the
code easy to understand, maintain, and extend. This document describes the
development conventions and review practices that support those goals. Read it
alongside [CONTRIBUTING.md](CONTRIBUTING.md) and the surrounding code.

## General design principles

- Understand the numerical method, existing behavior, and actual callers before
  changing an implementation. Preserve established interfaces and input files
  unless a change is explicitly approved.
- Start with the requested behavior and a few concrete examples. Choose the
  simplest design that meets those needs; do not build machinery for hypothetical
  future uses.
- Put an operation in the class or library that is responsible for it. Reuse
  suitable existing code and standard-library facilities instead of duplicating
  them or adding unnecessary dependencies. Share definitions of the same
  mathematical object across consumers rather than maintaining parallel enums;
  different numerical implementations may still be appropriate.
- Avoid unneeded ceremony: extra classes, wrappers, constants, comments, and
  checks should make the code clearer or serve a real purpose. Reasonable
  exceptions to conventions are welcome when they improve the result.
- Develop the implementation, tests, and documentation together. Verify the
  intended numerical or user-visible behavior, not just that the program runs.

## Repository and code organization

- `include/ibamr/` and `src/` contain IBAMR interfaces and implementations;
  `ibtk/include/ibtk/` and `ibtk/src/` contain the reusable toolkit. Inline
  implementations belong in the corresponding `private/` headers. Cartesian
  Fortran kernels often require `m4` preprocessing; edit their sources, not
  generated files.
- Consult [doc/cmake.md](doc/cmake.md), the current CMake files, and
  [.github/workflows/push-pull.yml](.github/workflows/push-pull.yml) for build
  configuration. Use installed dependencies and existing build environments;
  do not assume another developer's paths or an old copied command are current.
- Match surrounding naming: non-static members use `camelCase()`, static and
  free functions use `snake_case()`, member data use `d_`, and static member data
  use `s_`. Enum enumerators use uppercase names such as `SKEW_SYMMETRIC`; follow
  the existing `enum_to_string` and `string_to_enum` conventions when adding
  conversions. Follow the surrounding include style.
- Use `NDIM` for spatial dimension. Introduce a separate dimension template
  parameter only for an object whose dimension can actually differ from `NDIM`,
  not for hypothetical reuse.
- Use `auto` only when the exact type is evident from the declaration or
  right-hand side, such as `std::make_shared<Foo>()`. Use the explicit index type
  for `box.lower()`. Necessary deduction for lambdas or genuinely dependent
  template types is a narrow exception; do not add elaborate machinery or type
  erasure just to avoid `auto`.
- Keep headers self-contained so they compile independently. Include the
  corresponding `ibamr/config.h` or `ibtk/config.h`. Use `namespaces.h` and
  `app_namespaces.h` only in source files, not library headers.
- Put class and struct definitions in declaration headers; an implementation-only
  nested type may remain in its owner's private section. Ordinary `.cpp`-local
  helper types stay in their owning implementation. Public headers declare
  interfaces, not function bodies. Put ordinary function definitions in `.cpp`
  files; put visible template, `constexpr`, and inline definitions in matching
  private inline headers included by the declaration header. A one-line body is
  still an implementation.
- Keep algorithm-specific helper types in their owning implementation, not the
  public API. A source-private shared header is appropriate for two genuine
  implementation consumers. Do not install it merely to make inclusion easier.
  Reusable kernel evaluation, for example, need not belong to PETSc matrix
  utilities.
- Pass only the information a helper needs; use an input filename rather than
  `argc` and `argv` when that is sufficient. A helper used only in one `.cpp` and
  needing no class access normally belongs in its unnamed namespace.
- Avoid `friend`. Prefer operations on the owning class or ordinary interfaces
  that belong to its contract. Do not add friendship or public accessors solely
  for tests. A required C++ customization or external interface may justify a
  narrow, documented exception.
- Order declarations `public`, `protected`, then `private`; keep definitions in
  declaration order. Avoid public member data. In declarations, omit top-level
  `const` on by-value parameters; definitions may use it for unchanged values.

Assign directly to members rather than adding aliases such as `auto& dx = d_dx;`
just to rename them. Spelling an explicit type for the same alias does not help:

```cpp
for (unsigned int d = 0; d < NDIM; ++d)
{
    d_dx[d] = dx0[d] / static_cast<double>(ratio(d));
}
```

## Comments and class documentation

- Write comments for information a caller or maintainer needs. Use direct
  mathematical and library terminology; explain unfamiliar terms when needed.
  State concrete behavior, requirements, or reasons for non-obvious choices.
  Omit development history, hypothetical extensions, and prose that merely
  restates the declaration or advertises the design.
- A class header should explain what the class represents and how to use or
  extend its interface. Include mathematical meaning, inputs, results, limitations,
  ownership, lifetime, invalidation, and calling requirements as needed. Describe
  allocation, caching, storage, or forwarding only when it establishes a requirement
  or guarantee the reader needs. Include a small, complete example when it makes
  usage clearer.
- Review API documentation by asking, in order: "Would removing this sentence leave
  someone unsure how to use or extend the class correctly?" If so, "Is that
  information already explained where the API is defined, usually in a base class?"
  These questions concern API documentation, not helpful implementation comments.
- Document shared API requirements where the API is defined. Derived classes
  explain particular behavior, additional requirements, and meaningful differences
  without repeating the base class. Explain requirements spanning several methods
  once at class or group level. Keep shared configuration requirements in the
  relevant configuration documentation rather than repeating them in constructors.
- If shared documentation is missing or unclear, consider improving the base class.
  This normally belongs in a separate focused PR unless the clarification is needed
  to document the current change correctly. Do not silently broaden the change or
  defer a required clarification to later cleanup.
- Compare documentation of sibling classes implementing the same API while reviewing
  the current change. Comparable behavior merits comparable detail, not identical
  wording or length; extra detail should explain genuine differences. If siblings
  need the same explanation, consider putting it in the base class. This does not
  require rewriting every sibling.
- Document member functions, including private ones. A concise Doxygen sentence
  is enough for a simple helper; documenting every member does not mean repeating
  class documentation. Use Doxygen links, `\see`, or `\copydoc` when useful, with
  valid targets; no automatic cross-reference is needed at every override. State
  additional requirements of an overload beside that overload.
- Explain implementation choices and non-obvious mathematical relationships beside
  the implementation; keep private-data descriptions brief and local. Omit obvious
  mechanics, commented-out code, and obsolete explanations.

For example, declare a simple query in the public header and define it in the
corresponding implementation file:

```cpp
/*! \brief Return the number of levels in the hierarchy. */
int getNumberOfLevels() const;
```

The short description above is sufficient. A comment explaining why a borrowed
matrix must outlive a solver is useful; `// Set the matrix` above
`setMatrix(matrix)` is not.

For helpers that assemble matrix rows or retain an evaluator, make the operation
explicit, for example:

| Before | After |
| --- | --- |
| Assemble one velocity component with compile-time direction choices. | Assemble interpolation-matrix rows for velocity component `Axis`. |
| Bind immutable evaluator storage to generic matrix construction. | Return a matrix builder that owns a const evaluator. |

### Name meaningful constants, not every literal

Use a constant when the name explains domain meaning or keeps related uses in
sync. Prefer `UPPER_CASE_SNAKE_CASE` at the narrowest useful scope, subject to
established local conventions. Do not replace constants with macros or expose
them publicly just for a test.

```cpp
constexpr int MAX_RETRIES = 3;
for (int attempt = 0; attempt < MAX_RETRIES; ++attempt)
{
    // ...
}
```

Keep ordinary `0`, `1`, and clear mathematical literals: `values.size() + 1`
for a terminator does not need `constexpr int ONE = 1`. When a policy needs
explanation, explain why it exists instead of restating its bound or current
numeric value. Stable mathematical facts, such as a three-dimensional formula,
may of course contain numbers. Do not rename every literal in a file while
fixing one function.

The same judgment applies beyond constants: use a small helper for a genuinely
repeated operation, not a wrapper that merely renames a single comparison.

## Errors, ownership, and numerical behavior

- Prefer validating construction-time arguments and invariants during construction;
  omit repeated checks of unchanged validated state. Continue checking new
  arguments and state or lifecycle requirements that can change. Mutations must
  preserve invariants.
- Use `static_assert` at the earliest meaningful compile-time scope for compile-time
  requirements. Retain `static_assert` checks of function-template arguments only
  known at use; they have no runtime overhead.
- Use `TBOX_ERROR()` for fatal IBAMR runtime errors, with the surrounding
  stream-style formatting. Avoid C++ exceptions and catch-and-rethrow/translate
  scaffolding unless an existing external interface requires that boundary.
  Retain `TBOX_ERROR()` for invalid external input that must fail in Release.
- Call `TBOX_ERROR()` on the rank that detects the error. Do not broadcast an
  error string or add a collective just so all ranks report it. Communication
  needed to establish successful shared state is a different matter.
- Handle PETSc failures through `IBTK_CHKERRQ()` or the established local macro.
  Make owned versus borrowed objects explicit, preserve lifetime requirements,
  and release owned resources through the normal teardown path.
- Preserve numerical contracts, signs, scaling, and boundary conditions unless
  intentionally correcting them. Explain and test an approved correction rather
  than preserving a known defect or silently changing an expected result.
- Review frequently called code affected by the change for repeated conversions,
  lookups, allocations, and indirect calls. Move invariant work outside loops
  when practical, and use profiling when an uncertain cost could matter to the
  change. This does not require a benchmark for every PR or prohibit indirect
  calls where they are appropriate.

```cpp
// Prefer a direct fatal error at the point of detection.
if (width <= 0)
{
    TBOX_ERROR("Interpolator::initialize():\n" << "  kernel width must be positive\n");
}
```

Avoid throwing `std::runtime_error`, catching it in the caller, broadcasting
`what()`, and finally calling `TBOX_ERROR()` for the same invalid width.

## Configuration and small utility functions

- Check existing input keys and their meanings before adding a new mechanism.
  Keep related settings with the component they configure and preserve established
  precedence, including command-line PETSc overrides.
- Validate the rules IBAMR introduces. Let SAMRAI handle database syntax and PETSc
  handle PETSc option syntax and values instead of implementing a second parser
  or a more restrictive set of rules in IBAMR.
- Prefer direct processing when no intermediate representation is needed. Use a
  standard-library conversion or algorithm when it fits, rather than adding a
  custom text helper or Boost usage. Check that conversion preserves the intended
  value; do not add normalization or locale handling without a concrete need.
- Prefer `std::filesystem::path` consistently for filesystem operations, converting
  to strings where an interface requires them. Use nonthrowing overloads when
  handling filesystem errors through `TBOX_ERROR()`, and check the reported error.

| Prefer | Avoid |
| --- | --- |
| Walk a database directly when finding and applying settings is the whole operation. | Collect, sort, normalize, and replay settings without an ordering requirement. |
| Add a distinct input key for a new mechanism. | Reinterpret an existing key such as `petsc_options_prefix`. |
| Use `std::to_chars` with appropriate precision to format a floating-point setting. | Use `std::to_string` for a tiny tolerance without checking whether it becomes `"0.000000"`; write a custom numeric formatter unnecessarily. |
| Use an enum for a genuinely closed set of choices. | Force an extensible application-supplied catalog into a closed enum. |

## Focused, trustworthy tests

- Include a focused native `attest` regression in the same PR as each feature or
  bug fix. Reuse an existing executable with another case when initialization and
  linkage are shared. Do not add a separate CTest test in place of native coverage.
- Build the executable and use the normal fixture-linking targets. A CMake target
  alone, or a hand-copied test collection, is not integrated coverage. Verify
  `attest -N` discovers the intended input/output pairs before running them.
- Test the public operation using real library objects. An independent small
  reference calculation can check the answer, but a surrogate implementation
  cannot stand in for the IBAMR operation under test. Do not expose production
  internals or add general export/replay machinery just to simplify a test.
- Keep cases deterministic, repository-local, network-free, and small. Use legal
  SAMRAI patch geometry and no host-specific output. New cases should each take
  less than a minute in Release. Cover distinct contracts and meaningful edge
  cases, not every imaginable internal state.
- For numerical results, normally write the values to `output` through `plog` and
  supply the matching expected output. `attest` uses `numdiff` to compare numbers
  with tolerances that allow small differences between processors, compilers,
  and dependency versions. Prefer this to replacing useful numerical values with
  Boolean `PASS` markers or duplicating the comparison in a custom test helper.
- Keep numerical output short, labeled, and readable by a reviewer. Select a few
  meaningful values or norms instead of dumping large vectors and matrices. Print
  enough digits to expose relevant differences. Use an additional numerical check
  when a required property is not adequately checked by the output comparison.
- Use existing assertions for suitable non-numerical conditions instead of a
  bespoke Boolean-checking wrapper. Ensure required checks remain active in the
  tested configuration: a Debug-only assertion is not Release coverage. An
  assertion-only case can have an empty expected output; it does not need a
  transcript of `PASS` messages.
- An expected-error case must compare the actual, stable diagnostic against a
  nonempty expected output. Reuse logging helpers such as `TestAppender` in
  [tests/tests.h](tests/tests.h) where suitable. If the operation unexpectedly
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
with a small, recoverable negative control, then restore the fixture; no permanent
mutation-testing framework is needed.

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

## Build efficiently without weakening checks

- Reuse persistent out-of-tree Debug and Release builds. Verify their source
  paths, compilers, dependencies, and effective compile commands. "Fresh results"
  means current-source results, not an empty build directory for every commit.
- Use `ccache` for supported compilers and configure the C/C++ launchers explicitly
  (`CMAKE_C_COMPILER_LAUNCHER`, `CMAKE_CXX_COMPILER_LAUNCHER`). Use consistent,
  writable `CCACHE_DIR` and `CCACHE_TEMPDIR` locations for configuration and builds,
  with `CCACHE_BASEDIR` set to the persistent source root. Retain a sufficiently
  sized populated cache across related branches and maintenance work.
- Verify actual cache activity with `ccache -p` and `ccache -s -v`; installation
  alone is not use. An up-to-date target makes no cache calls. Do not clear useful
  caches, force recompilation to manufacture hits, or relax correctness settings.
- Match applicable CI warning/error settings in both library and test builds,
  including Debug and Release. Read the current workflow rather than freezing
  its flags here. Verify effective flags, not just requested CMake values.
  Compiler-specific flags and dependency exceptions need equivalent intent,
  not blind copying between Clang, GCC, and Fortran.
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
- Report commands, configuration, tested revision, and limitations concisely.
  A 3D library build is not a 3D test run; a simulated compatibility compile is not
  testing an old dependency installation. Local passes are not hosted CI passes.
  Record test duration for practicality, not as solver-performance evidence.

## Commits, PRs, and review

- Before editing, read applicable local instructions and verify the branch, base,
  worktrees, and dirty state. Preserve unrelated work and explicitly frozen
  references. Each PR should have one clear purpose that can be explained in one
  or two sentences, including its necessary tests and documentation. An
  options-loading fix need not also reorganize solver classes; necessary supporting
  changes are welcome, but unrelated cleanup belongs in a separate change.
- Check copyright starting years only for files created in the branch. Brand-new
  files use the current year; files copied and modified from an existing file
  retain that file's starting year. For example, a new 2026 file uses `2026`, while
  a modified copy of a file starting in 2018 uses `2018 - 2026`. Do not investigate
  history or repair starting years in existing files. End years are updated
  automatically when tagging a release; updating them on files actually modified
  by the branch is also appropriate. Leave unrelated files unchanged; do not do
  calendar-only, tree-wide cleanup.
- Before finalizing source commits, run `make indent` in the corresponding
  configured build directory. Use the clang-format version required by IBAMR,
  not whichever version is newest or first on `PATH`. The version check is in
  [scripts/formatting/indent_common.sh](scripts/formatting/indent_common.sh).
  IBAMR provides `scripts/formatting/download-clang-format` and
  `scripts/formatting/compile-clang-format` to obtain a compatible formatter when
  needed. See [scripts/formatting/README.md](scripts/formatting/README.md); the
  source-tree formatting entry point is `scripts/formatting/indent`.
  Inspect all formatter changes and keep unrelated edits out. Markdown-only
  changes need document checks, not a library rebuild or a claim that source
  formatting was run.
- Run `git diff --check`, stage only intended paths/hunks, and inspect both the
  staged diff and complete PR diff, including comments and fixtures. Review changed
  comments for clarity and accuracy; the presence of a Doxygen block is not
  sufficient. Keep generated artifacts and logs out of commits. Preserve the
  contributor's established Git identity; do not invent author information or
  assistance/co-author tags.
- Keep PR descriptions to one or two concrete sentences. Retain the repository's
  PR template when applicable and keep detailed build/test evidence in the task
  return, not a long PR body.
- When a review comment reveals a repeated issue, check other instances introduced
  by the change. Remove unused helpers and obsolete tests or explanations left
  from abandoned approaches, keeping the cleanup within the change's scope.
  Retain tests that still cover supported behavior.
- Review comments are input to an in-scope revision, not authority for unrelated
  redesign. Respect explicit review-only holds. Authorizations to commit, push,
  open a PR, reply to reviews, rewrite history, and merge are distinct: a request
  to push a branch does not by itself authorize opening a PR. A request to inspect
  or propose changes does not authorize publishing them.
- When delegating, provide this guide and any task-specific coding, testing,
  documentation, and PR requirements, together with the user's actual limits.
  Do not add publication authority; a later explicit user hold overrides an
  earlier delegated instruction.
- For stacked work, inspect each unit relative to its intended predecessor,
  including comments, helpers, and fixtures. After an authorized rewrite, preserve
  recovery refs, replay only the child's own interval, coordinate with its owner,
  and reverify the incremental diff and affected tests. Use an explicit
  expected-old-SHA `--force-with-lease` when a force push is authorized.
  Check actual GitHub bases and ancestry after a merge; do not assume automatic
  retargeting or restacking. Never modify another task's worktree to synchronize it.

### The PR checklist

Use [doc/pull_request_template.md](doc/pull_request_template.md) when applicable.
For documentation-only changes, omit the checklist and briefly state why it does
not apply.

Consider each item: the test suite, formatting, relevant issues, a changelog entry,
restart compatibility, tests for new functionality or fixes, and GitHub labels
when permitted. Retain the checklist as boilerplate: change checkbox states only,
without adding explanations between or after its items. Check only supported
items or items explicitly dismissed by a principal developer. Focused tests
alone do not justify checking the full-suite item.

Changelog entries describe the net user-visible change, not development history.
Do not add separate entries for problems introduced and corrected within the same
PR or unmerged stack. Update the feature's existing entry when needed; a fix to
pre-existing IBAMR behavior may warrant its own entry.

Changelog entries should be very short; for example,
`Fix: Preserve command-line PETSc option precedence.` Follow nearby dating and
attribution conventions. Check serialization and restart versioning when
persistent state changes; transient scratch state alone does not require a
version increment. Cite relevant issues, using `Fixes #...` only when the PR
actually resolves the issue, and set appropriate labels when authorized.
