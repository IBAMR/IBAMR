# Review, delegation, and stacked work

## Review and delegation

- When delegating, provide this guide and any task-specific coding, testing,
  documentation, and PR requirements, together with the user's actual limits.
  Do not add publication authority; a later explicit user hold overrides an
  earlier delegated instruction.

For an independent review, provide an exact revision or immutable snapshot,
scope, settled decisions, and explicit limits. Reviewers should remain read-only
unless editing is explicitly part of their assignment. Verify actionable findings
against the code before implementing them; a reviewer suggestion is not authority
for unrelated redesign. Save a detailed report when useful, with a concise task
summary.

## Stacked changes

- For stacked work, inspect each unit relative to its intended predecessor,
  including comments, helpers, and fixtures. After an authorized rewrite, preserve
  recovery refs, replay only the child's own interval, coordinate with its owner,
  and reverify the incremental diff and affected tests. Use an explicit
  expected-old-SHA `--force-with-lease` when a force push is authorized.
  Check actual GitHub bases and ancestry after a merge; do not assume automatic
  retargeting or restacking. Never modify another task's worktree to synchronize it.

Verify each affected head's incremental diff, formatting, fixture discovery, and
assigned tests after a restack; testing only the tip does not establish intermediate
heads' correctness. Ensure each expected output matches the code at that head.
State dependencies and the intended base in the PR and confirm that CI actually
runs on the affected PRs.

Keep a bug discovered during stack work as a focused fix with its own regression
and changelog when needed. Coordinate dependency changes with owners. For a fix
landing on the base branch first, inspect conflicts read-only using isolated
snapshots; patch applicability alone does not establish semantic compatibility.
Do not modify another task's worktree. Choose restack commands from the actual
graph and ownership constraints rather than applying a universal script.
