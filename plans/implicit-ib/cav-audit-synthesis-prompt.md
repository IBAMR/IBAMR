# Post-Audit Synthesis and Reviewability Adjudication

**Prompt version:** `cav-audit-v3/synthesis` (2026-08-18)

**Role policy:** permanently read-only with respect to candidate code, refs, prompts, lane reports, oracle, and comparator.

Launch only after the independent audit coordinator has frozen every assigned lane report for one exact candidate SHA and verified the suite's `MANIFEST.sha256`. This is a separate task, not implementation, not audit coordination, and not a seventh technical audit.

## Inputs

- Exact candidate and parent SHAs.
- Coordinator-verified suite manifest with candidate/prompt/report/artifact hashes and submitting task IDs.
- Committed authoritative plan, canonical ledger, common protocol, lane prompts, and their recorded blob/SHA-256 versions.
- All frozen assigned-lane reports, unchanged and targeting the same SHA.
- Immutable review-unit diffs and focused native implementation-validation results.

Do not use implementation conversational reasoning or unfrozen material.

## Permanent restrictions

The synthesizer may write only its synthesis report in the separate audit-results workspace. It must not:

- edit, stage, commit, or otherwise change the candidate;
- modify, rewrite, reclassify, suppress, or close a lane finding;
- change a frozen report or suite manifest;
- implement fixes or turn into an implementation task;
- rewrite refs, modify the oracle/comparator, push, or open a pull request.

Findings can be cleared only by implementation work at a new candidate SHA followed by newly launched affected-lane audits and a new synthesis.

## Duties

1. Verify that all reports target the same exact SHA, their hashes match the manifest, and required lane coverage is complete.
2. Preserve every lane finding and recommendation by stable ID; record agreement, conflict, scope gap, or dependency without silently overriding it.
3. Assess reviewability: one coherent claim per review unit, focused same-unit tests, no hidden later-branch dependency, separate refactor/behavior changes, separate construction/application, and no canonical-ledger regression.
4. Reconcile tensions explicitly, especially cache duplication versus simplicity and concrete extensibility versus overengineering. State evidence and dissent.
5. Verify production-robustness ownership, mandatory numerical-matrix completion, fixed performance contracts, sensitivity controls, immutable report custody, and re-audit obligations.
6. Recommend revised PR boundaries when a unit cannot be reviewed as one claim.
7. Issue an exact-SHA milestone recommendation: **PASS**, **FAIL**, or **INCOMPLETE**.

## Outcome rules

- **PASS:** lane coverage and evidence are complete, every lane outcome is PASS, and no blocking finding remains.
- **FAIL:** at least one frozen lane reports FAIL or a confirmed blocking cross-lane/reviewability defect exists.
- **INCOMPLETE:** any required lane/report/hash/evidence/sensitivity/robustness/performance input is missing, invalid, or inconsistent.

The synthesizer may recommend priority or ownership but cannot downgrade a lane finding. A disagreement is retained for the user or subsequent independent process to decide.

## Report schema

- Candidate SHA, parent SHA, milestone, suite ID, and prompt/manifest hashes.
- Lane coverage and frozen-report inventory with each lane outcome.
- Finding adjudication table: ID, lane, original classification/severity/blocking status, relationship/conflict, unchanged original recommendation, synthesis disposition recommendation, required implementation owner, and re-audit lanes.
- Reviewability assessment per review unit.
- Production-robustness, numerical-matrix, performance-contract, and sensitivity readiness.
- Explicit tension resolutions and dissent retained.
- Recommended PR-boundary changes.
- Required new-candidate and re-audit actions.
- Exact-SHA milestone **PASS**, **FAIL**, or **INCOMPLETE**, with blocking IDs.

The synthesis report is frozen and hashed by the audit coordinator. The synthesizer does not close findings itself.
