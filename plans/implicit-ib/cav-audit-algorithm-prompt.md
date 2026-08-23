# Audit Lane A: Algorithm Semantics

**Prompt version:** `cav-audit-v8/A` (2026-08-21)

Use with `cav-independent-audit-prompt.md`. This lane is independent and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

Use the common protocol and implementation-plan Section 7 for the canonical focused-versus-broad `attest` evidence contract; do not redefine or substitute those scopes in this lane.

## Question

Does the candidate implement the mathematically correct action of global multiplicative CAV and, when in scope, macro-subdomain CAV-RAS?

Independently derive both algorithms before tracing code. The sandbox is executable evidence, not the definition of mathematical correctness.

## Required evidence

- Derive the patch-restriction/injection form of multiplicative CAV, including original residual initialization, current-residual evolution after every patch, relaxation, forward/reverse/palindromic traversal, and pressure gauge.
- Derive macro-subdomain CAV-RAS with one common outer residual, an independent zero-state multiplicative CAV sweep on each overlapping solve set, restriction of the completed macro correction to its unique owned set, and outer additive assembly.
- Trace every derived step through production code and focused native tests.
- Use the formal authorities pinned in the common protocol: `docs/CAV_ALGORITHMS.tex`, `docs/CAV_RAS_PARALLEL_NOTES.md`, and the named manuscript sections at their exact commits/blobs.
- Verify the implemented resolution of the seed mismatch: the formal RELAXED path must be pressure-cell seeded, use only row-or-column `E_h` adjacency, retain exact standard patches when no coupling leaves `U_i`, and use the full original saddle-point matrix only for local solves/residual updates. Audit the declared pressure-seeded STRICT extension separately and ensure any foundation full-`A00` velocity-seeded compatibility path is explicitly named legacy behavior. Do not presume family or sweep equivalence.
- Show fresh `b-Ax` checks at declared boundaries and prove that pressure normalization changes only the pressure gauge.
- For RAS, verify macro geometry/order, solve and owned sets, exact ownership coverage, overlap, patch order within each macro, local residual/state traces, restricted completed correction, and common-state outer assembly.
- Cover advertised MPI scope and boundary/no-IB cases where partitioning could alter semantics, using only hierarchy patches legal for the neighboring SAMRAI/IBAMR path. In 2D every hierarchy patch extent is at least `4x4`; preserve the intended semantic condition within that supported geometry. Do not require a CAV-specific sub-minimum path. This platform constraint does not restrict algebraic CAV-patch membership.
- Exercise and detect activated stale-residual, patch-order, pressure-gauge, RAS restriction/ownership, overlap, and inner-additive-for-inner-multiplicative mutations, with unmutated controls.
- Verify that comments on mathematical state transitions, ordering, restriction, residual evolution, and gauge/nullspace behavior accurately describe the implemented invariants; materially misleading statements are defects, not style preferences.

## Exclusions

- Do not infer correctness from final convergence alone.
- Do not define correctness by copying sandbox or historical code.
- Do not evaluate container taste or performance except where it changes mathematical action.
- Do not accept flat additive-then-restrict patch processing as CAV-RAS.
- Do not treat normal SAMRAI/IBAMR rejection of an invalid hierarchy-patch size as an algorithm defect or accept padding, splitting, synthetic levels, alternate initialization, or a special CAV path as an algorithmic extension.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires a complete independent derivation, seed-formulation reconciliation, production trace, focused tests, supported-scope checks, and all applicable sensitivity controls with no unresolved semantic gap. Any confirmed stale residual, wrong traversal, gauge corruption, incomplete/overlapping ownership, inner-additive substitution, or mislabeled algorithm is **FAIL**. Missing formal reconciliation, inactive sensitivity controls, or unavailable required evidence is **INCOMPLETE**.

## Lane-specific report schema

After the common schema, include: derived equations/pseudocode; authority comparison; formal-pressure/legacy-velocity construction table; `E_h` provenance and row-or-column trace; code-to-step trace; invariants; supported dimension/MPI scope; traversal/composition cases; sensitivity outcomes; earliest semantic mismatch; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
