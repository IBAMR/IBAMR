# Audit Lane E: Concrete Extensibility

**Prompt version:** `cav-audit-v5/E` (2026-08-20)

Use with `cav-independent-audit-prompt.md`. This lane is independent and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

Use the common protocol and implementation-plan Section 7 for the canonical focused-versus-broad `attest` evidence contract; do not redefine or substitute those scopes in this lane.

## Question

Does the candidate support the known planned variations through clear, bounded seams without building a hypothetical framework?

## Required evidence

Assess only these axes:

- RELAXED versus STRICT patch construction;
- global multiplicative versus macro-subdomain CAV-RAS composition;
- supported local-solver backends and their failure policies;
- forward/reverse traversal and palindromic use;
- eventual serial versus distributed execution.

For each in-scope axis, identify its selection point, interface/data contract, state ownership, tests, and whether adding the next planned variant requires modifying unrelated layers. Verify the four-layer decomposition: construction, local solve, correction composition, FAC/V-cycle orchestration.

Perform a structured, no-code change exercise for each concrete axis: list the minimal files/interfaces/tests needed for the named next variant. Also assess a seeded cross-layer-entanglement example supplied in the disposable audit fixture or protocol exercise and require the lane to reject it. Record the exercise version, activation, unmutated control, and result.

For serial/distributed evolution, inspect whether ownership and residual-gather contracts are explicit without requiring current implementation of a later distributed milestone.

Verify that public API and nearby comments accurately document only the concrete supported seams, ownership/failure behavior, and extension boundaries; speculative framework commentary or a documented seam absent from code is a defect.

## Exclusions

- No hypothetical plugin framework, registry, policy hierarchy, or public API without a named planned consumer.
- Do not require current implementation of later milestones.
- Do not count arbitrary configurability as extensibility.
- Do not excuse duplicated semantic paths merely because they are selectable.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires narrow tested seams for every in-scope concrete axis, a credible bounded next-change exercise, successful rejection of seeded entanglement, and no unused general framework. A planned axis entangled across unrelated layers, semantic duplication, or substantial hypothetical machinery is **FAIL**. Missing interfaces/tests/exercises or an inactive structured challenge is **INCOMPLETE**.

## Lane-specific report schema

After the common schema, include: axis-by-seam matrix; exact next-change exercises; selection/state/test locations; cross-layer edits required; seeded-entanglement result; serial/distributed boundary; overgeneralization findings; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
