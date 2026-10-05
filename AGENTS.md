# IBAMR development guidance

IBAMR provides immersed boundary methods and related numerical methods on
adaptively refined Cartesian grids. Its implementation combines C++ library
code with Fortran numerical kernels and uses SAMRAI, PETSc, and other scientific
libraries.

Changes should preserve numerical correctness and efficiency while keeping the
code easy to understand, maintain, and extend. This document describes the
development conventions and review practices that support those goals. Read it
alongside [CONTRIBUTING.md](CONTRIBUTING.md) and the surrounding code.

## Workflow checklist

- Read the topic files below that apply to the task; ordinary links are not a
  guarantee that an agent has loaded their contents. Re-read them when scope
  expands. Read [CONTRIBUTING.md](CONTRIBUTING.md) alongside this guide.
- Before editing, confirm the intended checkout, branch or revision, base, scope,
  and dirty state. Preserve unrelated work, other tasks' worktrees, and explicitly
  frozen refs.
- Reuse a suitable native test or explain why a new executable is needed. Verify
  fixture discovery and assigned Debug coverage. For bug fixes, normally show
  failure without the fix and success with it; explain impractical comparisons.
- Investigate changed numerical baselines; never replace expected results merely
  to make a test pass. Report changed baselines and their justification.
- For source changes, read [scripts/formatting/README.md](scripts/formatting/README.md),
  use the formatter required by `scripts/formatting/indent_common.sh`, and inspect
  the changes. Markdown-only work needs document checks, not library builds.
- When a changelog is needed, first read
  [doc/news/changes/README.md](doc/news/changes/README.md). Credit human contributors.
- Run `git diff --check`; inspect the intended staged diff and complete change,
  including comments, fixtures, copyright, and attribution. Stage only intended
  paths or hunks. Keep logs and generated artifacts out of commits.
- Run `scripts/maintenance/check_pull_request.py` on the branch and the intended
  title and description, and resolve or be ready to explain what it reports. It
  is a screen for common problems, not a substitute for reading the change.
- Before any publication, verify that an explicit user request covers the action
  and destination. Use concise, accurate commit and PR text and the applicable
  complete [.github/pull_request_template.md](.github/pull_request_template.md).
- Check conditional checklist items after verifying they do not apply; check the
  issue-citation item when no relevant issue exists. Leave the full-suite item
  unchecked until CI passes on the anticipated final PR revision.
- After an authorized push, check CI on that exact commit and report its state
  accurately. After a restack, verify every affected head, not only the tip.

The task report is the account of the work that an agent gives the person who
requested it; it is not part of the commit or the PR.

## Authorization and attribution

- Opening or updating a PR requires an explicit user request. Do not infer
  permission from a coding task, a review finding, account access, or having
  opened the PR earlier. A request to prepare or inspect a change is not a request
  to publish it.
- A request to open, update, or maintain a PR covers relevant existing labels,
  body corrections, and normal follow-up pushes within its stated scope. Before
  each action, confirm that the request covers it; do not treat a completed request
  as indefinite maintenance authority. Honor existing explicit grants without
  asking again for already-authorized actions. Explicit review-only holds prevail.
- Commit, push, open/update a PR, reply to reviews, rewrite history, and merge are
  distinct authorizations. A push request alone does not authorize opening a PR.
  History rewriting, force pushing, and merging require separate authorization,
  including when the agent created the branch or PR. Preserve recovery refs and
  use `--force-with-lease` pinned to the expected old SHA for any authorized force
  push. Stack relinking also needs authorization for the affected PR changes.
- Use descriptive `agent/...` names for new agent branches; preserve existing
  names. The publication authorization must specify the remote: do not assume
  upstream or a fork.
- Preserve established Git identity and genuine human attribution. Do not add
  automated `Co-Authored-By` trailers or `Generated with ...` credits to commits
  or PR bodies; this rule overrides optional tool defaults. Add `Signed-off-by`
  only when explicitly required or requested. Do not invent author information.

## Required topic guidance

Read every applicable file before doing the listed work. These files preserve
this guide's detailed conventions; they are not optional background reading.
`CLAUDE.md` is a relative symlink to this file so both entry points share one text.

For reviews, read the topic files applicable to the changes being reviewed and
[Review and stacked work](doc/agents/review-and-stacked-work.md).

| Before you... | Read... |
| --- | --- |
| Edit C++, headers, Fortran, or m4 kernels | [Code style and organization](doc/agents/code-style.md) |
| Change numerical behavior, ownership, errors, configuration, or utilities | [Errors and numerics](doc/agents/errors-and-numerics.md) |
| Write or review API comments and documentation | [Documentation](doc/agents/documentation.md) |
| Add or modify tests, fixtures, or expected output | [Testing](doc/agents/testing.md) |
| Configure, build, format, or report verification | [Build and verification](doc/agents/build-and-verification.md) |
| Prepare commits, changelogs, or PRs | [Commits and pull requests](doc/agents/commits-and-pull-requests.md) |
| Delegate a review or work on dependent branches/history rewrites | [Review and stacked work](doc/agents/review-and-stacked-work.md) |

## Verification defaults

Use a supported compiler toolchain and compatible dependencies for the available
environment. Existing builds, compiler caches, and dependency sources can help
when available; see [Build and verification](doc/agents/build-and-verification.md)
for suggested workflows. Use the clang-format version pinned by the repository
and required by CI.

Use Debug for routine implementation, review maintenance, and restacking. Complete
assigned dimensions, discovery, focused tests and reruns, and assigned broad
suites. Defer Release builds, tests, and runtime-budget measurements to final
working-stack integration/performance work or an explicit targeted request.
Retain existing Release builds, caches, and results. Pending Release qualification
does not block implementation handoff or otherwise-authorized publication; report
it as pending, not passed. Local success does not establish hosted CI success.

## General design principles

- Understand the numerical method, existing behavior, and actual callers before
  changing an implementation. Preserve established interfaces and input files
  unless a change is explicitly approved.
- Start with the requested behavior and a few concrete examples. Choose the
  simplest design that meets those needs; do not build machinery for hypothetical
  future uses.
- Where existing conventions are inconsistent, prefer modern C++ within the
  supported language standard and choices that improve clarity and efficiency.
  Preserve required interfaces and keep modernization within the change's scope.
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
