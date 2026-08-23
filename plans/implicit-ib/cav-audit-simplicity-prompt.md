# Audit Lane S: Simplicity and Minimality

**Prompt version:** `cav-audit-v8/S` (2026-08-21)

Use with `cav-independent-audit-prompt.md`. This lane is independent, deletion-oriented, and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

Use the common protocol and implementation-plan Section 7 for the canonical focused-versus-broad `attest` evidence contract; do not redefine or substitute those scopes in this lane.

## Question

Is the candidate the smallest clear implementation that preserves required semantics, evidence, structurally efficient data paths, and concrete extension axes?

## Required evidence

- Audit every production-code unit where complexity is introduced, including `cav-2`, `cav-4`, `cav-6a` through corrective `cav-6f`, `cav-7`, `cav-8a`, and `cav-8b`, plus `cav-1a/1b`, `cav-3a/3b`, or `cav-9a/9b` wherever those units add APIs, state, options, or duplicate mechanisms.
- Inventory every new abstraction, class/helper, state variable, cache/view, generic mechanism, public/private option, compatibility path, and duplicate traversal/solve/export path.
- Tie each retained element to a requirement, focused test, measured hot-path benefit, or concrete planned extension.
- Look for dead state, redundant representations, pass-through wrappers, premature genericity, duplicate setup/application code, and audit-only mechanisms leaking into production APIs.
- Recheck the canonical historical ledger and reject reintroduced duplication or mutable cache internals.
- Explicitly ask whether roughly 30% of the new code could be deleted while preserving requirements, tests, structural-efficiency invariants, and planned axes. The percentage is a challenge, not a quota.
- Sketch the smallest credible alternative for every material complexity finding.
- Assess a seeded redundant abstraction/state/option in a disposable review fixture or structured challenge and demonstrate that it is identified while the unmutated control is not falsely failed.
- Challenge every new test executable, helper, fixture, and audit hook: prefer an existing common executable with separate cases unless a documented technical seam requires another target, and reject duplicated harness logic or production APIs that become a testing framework.
- Challenge generated-looking, excessive, redundant, obvious, speculative, or drift-prone comments while distinguishing materially misleading statements from nonblocking comment-density preferences.
- Explicitly challenge and delete any padding, synthetic splitting or hierarchy level, alternate initialization, shim, compatibility option, or special solver path whose only purpose is to admit a hierarchy patch smaller than SAMRAI supports. Normal validation owns invalid requests; naturally small algebraic CAV patches are not such a workaround.

## Exclusions

- Do not require deletion that breaks measured cache behavior, precise semantics, focused testing, or a named planned axis.
- Do not treat short code as inherently clearer.
- Do not fold optional style cleanup into a blocking finding.
- Do not redesign hypothetical future systems.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires specific necessity for material complexity, no unjustified duplicate path, exposed option, test executable/helper, audit hook, or comment burden, completed per-unit and 30% challenges, canonical-ledger preservation, and a valid sensitivity control. Avoidable complexity that materially obscures review or creates parallel semantic/test paths is **FAIL**. Missing production-unit coverage, test/hook/comment inventory, deletion exercise, or active control is **INCOMPLETE**.

## Lane-specific report schema

After the common schema, include: per-review-unit production/test/hook/comment complexity inventory; executable/fixture reuse assessment; deletion candidates with estimated lines/state/options removed; smallest alternatives; 30% challenge result; ledger-regression check; seeded-control result; blocking versus optional cleanup split; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
