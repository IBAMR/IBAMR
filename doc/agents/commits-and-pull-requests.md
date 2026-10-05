# Commits, changelogs, and pull requests

Read the [root authorization rules](../../AGENTS.md#authorization-and-attribution)
before publishing. An explicit request must cover opening or updating a PR,
including follow-up pushes and metadata edits; preparing a change alone does not.
Use descriptive `agent/...` names for new agent branches, preserve existing branch
names, and use the remote specified in the publication authorization.

## Scope and review

- Before editing, read applicable local instructions and verify the branch, base,
  worktrees, and dirty state. Preserve unrelated work and explicitly frozen
  references. Each PR should have one clear purpose that can be explained in one
  or two sentences, including its necessary tests and documentation. An
  options-loading fix need not also reorganize solver classes; necessary supporting
  changes are welcome, but unrelated cleanup belongs in a separate change.
- Run `git diff --check`, stage only intended paths/hunks, and inspect both the
  staged diff and complete PR diff, including comments and fixtures. Check naming
  against the established convention in the same class, struct, or namespace,
  using the defaults when none applies. Check declarations and corresponding
  definitions and call sites; out-of-class definitions omit `static`, so inspecting
  function bodies alone can miss static member functions. Review changed
  comments for clarity and accuracy. Apply the
  [two documentation questions](documentation.md) to the actual changed
  API documentation, comparing it with documentation where the API is defined and
  in relevant sibling classes. Remove unnecessary inherited or internal repetition
  while preserving needed differences and caller/subclass guarantees. A Doxygen
  block alone is not sufficient. Keep generated artifacts and logs out of commits.
  Follow the attribution rules below.
- When a review comment reveals a repeated issue, check other instances introduced
  by the change. Remove unused helpers and obsolete tests or explanations left
  from abandoned approaches, keeping the cleanup within the change's scope.
  Retain tests that still cover supported behavior.

Record unrelated defects with their location, trigger or reproducer, and impact
in the task report; propose a separate follow-up instead of expanding this PR.
Publishing an issue also requires authorization. If the requested repair proves
insufficient, explain the additional repair needed for the same behavior. Include
limitations in the PR when they affect evaluating or using the change.

## Commit messages and attribution

Use a short imperative subject. A simple fix normally needs only the subject;
add a brief body when the rationale or constraints need explanation. Preserve
human authorship and genuine human co-author information. Do not add automated
`Co-Authored-By` trailers or `Generated with ...` footers. Configure optional tool
defaults to follow this rule. Follow [CONTRIBUTING.md](../../CONTRIBUTING.md),
including its Developer Certificate of Origin. Do not add `Signed-off-by` unless
explicitly required or requested; the DCO reference alone is not an instruction
to add a trailer.

## Copyright

- Check copyright starting years only for files created in the branch. Brand-new
  files use the current year; files copied and modified from an existing file
  retain that file's starting year. For example, a new 2026 file uses `2026`, while
  a modified copy of a file starting in 2018 uses `2018 - 2026`. Do not investigate
  history or repair starting years in existing files. End years are updated
  automatically when tagging a release; updating them on files actually modified
  by the branch is also appropriate. Leave unrelated files unchanged; do not do
  calendar-only, tree-wide cleanup.

Do not copy a duplicated year range such as `2026 - 2026` from a neighboring
file. Follow the repository's source-header form for new source files; do not
add source-code banners to `.input` or `.output` fixtures that do not use them.

## Titles and descriptions

Use short imperative PR titles, at most 72 characters, that name the fix or
addition. Start the title with `Fix`, after any tag, if and only if the change
corrects behavior that is wrong, and give such a PR the `Bug` label. A title
may begin with a bracketed tag that tells
reviewers which area the PR belongs to, such as `[FAC]`, `[HYPRE]`, or `[MPI]`.
Use a tag that open or recently merged PRs already use, and ask before
introducing a new one. The imperative title follows the tag. Begin the body
with one or two imperative sentences that state the change, in about 60 words
or fewer. Add only details a reviewer needs, such as material limitations,
API incompatibilities, stack dependencies/base, or an unusual validation limit.
Keep logs and detailed test reports in the task report. Check every claim against
the final diff and evidence, and update stale text within authorized PR work.

State the change itself, not the history of the defect or of the branch. Do not
list or describe the tests in the body, because the checklist covers them; the
exception is a sentence giving the reason for a new test executable. For a
stacked PR, end the opening text with `Stacked on #N.`, naming the PR it is based
on, and omit notes about stack position that a restack would invalidate. When a
PR adds something that is first used by a later PR in the stack, say in one clause
why it exists and name that PR, so that it can be reviewed on its own.

Use [.github/pull_request_template.md](../../.github/pull_request_template.md), which
links to the canonical template in `doc/`. Preserve the entire template, including
its HTML comment; change checkbox states only, without inserting explanations
between or after its items. Put necessary details before the template. For
example, a compact bug-fix body begins with:

> Fix component bounds in 3D specialized-linear side refinement.

Then append the complete template. For documentation-only changes, omit the
checklist and briefly state why it does not apply.

## Checklist decisions

| Item | When to check it |
| --- | --- |
| Test suite | CI has passed on the commit anticipated to be the final PR version. Focused or local full-suite runs alone do not satisfy this policy. A new revision needs CI evidence for that revision. |
| Formatting | The required repository formatting has run and its changes have been reviewed. |
| Issues | Relevant issues have been cited, or there are no relevant issues. Use `Fixes #...` only when the PR actually resolves the issue. |
| Changelog | An appropriate entry is present, or review confirms the conditional requirement does not apply. |
| Restart | Required restart handling/versioning is correct, or review confirms the conditional requirement does not apply. Transient scratch state alone does not require a version increment. |
| Tests for changes | Appropriate regression/feature coverage is present, or review confirms the conditional requirement does not apply. Follow the deferred Release policy for runtime qualification. |
| Labels | Relevant existing labels have been applied within the authorized PR work, or the account lacks permission. |

A principal developer may explicitly dismiss an item. Do not infer evidence or a
dismissal. After an authorized push, inspect CI for the exact pushed commit and
report pending or failed checks accurately. Leave the test-suite box unchecked
until the final-version CI condition above is met.

When authorized to open or update the PR, apply relevant existing labels: for
example, `Bug` for a bug fix, `Tests` for test changes, and `C++`, `Fortran`,
`Build System`, `CI`, or `Documentation` for the affected area. Query the current
repository label list; do not create labels or infer publication permission from
account access.

## Changelog entries

Read [doc/news/changes/README.md](../../doc/news/changes/README.md) before writing an
entry. It is the source of truth for categories, filenames, dates, and attribution.
Choose a filename that the target branch does not already use, and check again
after a rebase, since another merged PR may have taken it.
Use the actual human contributor names, not an agent name or invented identity.
Use an established prefix such as `Fixed:`, `Improved:`, or `New:`; for example,
`Fixed: Preserve command-line PETSc option precedence.` Keep the entry very short,
name the relevant qualified symbol when useful, and omit implementation detail.

Describe the net user-visible change, not development history. Do not add separate
entries for problems introduced and corrected within the same PR or unmerged
stack. Update the feature's existing entry when needed; a fix to pre-existing
IBAMR behavior may warrant its own entry. An approved incompatible public-API
change needs an `incompatibilities/` entry and a concise explanation in the PR.
Check serialization and restart versioning whenever persistent state changes.
