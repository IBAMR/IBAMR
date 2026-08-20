# Audit Lane S: Simplicity and Minimality

**Prompt version:** `cav-audit-v3/S` (2026-08-18)

Use with `cav-independent-audit-prompt.md`. This lane is independent, deletion-oriented, and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

## Question

Is the candidate the smallest clear implementation that preserves required semantics, evidence, performance, and concrete extension axes?

## Required evidence

- Audit every production-code unit where complexity is introduced, including `cav-2`, `cav-4`, `cav-6a`, `cav-6b`, `cav-7`, `cav-8a`, and `cav-8b`, plus `cav-1a/1b`, `cav-3a/3b`, or `cav-9a/9b` wherever those units add APIs, state, options, or duplicate mechanisms.
- Inventory every new abstraction, class/helper, state variable, cache/view, generic mechanism, public/private option, compatibility path, and duplicate traversal/solve/export path.
- Tie each retained element to a requirement, focused test, measured hot-path benefit, or concrete planned extension.
- Look for dead state, redundant representations, pass-through wrappers, premature genericity, duplicate setup/application code, and audit-only mechanisms leaking into production APIs.
- Recheck the canonical historical ledger and reject reintroduced duplication or mutable cache internals.
- Explicitly ask whether roughly 30% of the new code could be deleted while preserving requirements, tests, fixed performance criteria, and planned axes. The percentage is a challenge, not a quota.
- Sketch the smallest credible alternative for every material complexity finding.
- Assess a seeded redundant abstraction/state/option in a disposable review fixture or structured challenge and demonstrate that it is identified while the unmutated control is not falsely failed.

## Exclusions

- Do not require deletion that breaks measured cache behavior, precise semantics, focused testing, or a named planned axis.
- Do not treat short code as inherently clearer.
- Do not fold optional style cleanup into a blocking finding.
- Do not redesign hypothetical future systems.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires specific necessity for material complexity, no unjustified duplicate path or exposed option, completed per-unit and 30% challenges, canonical-ledger preservation, and a valid sensitivity control. Avoidable complexity that materially obscures review or creates parallel semantic paths is **FAIL**. Missing production-unit coverage, complexity inventory, deletion exercise, or active control is **INCOMPLETE**.

## Lane-specific report schema

After the common schema, include: per-review-unit complexity inventory; deletion candidates with estimated lines/state/options removed; smallest alternatives; 30% challenge result; ledger-regression check; seeded-control result; blocking versus optional cleanup split; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
