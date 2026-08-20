# Audit Lane C: Cache and Data-Layout Efficiency

**Prompt version:** `cav-audit-v3/C` (2026-08-18)

Use with `cav-independent-audit-prompt.md`. This lane is independent and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

## Question

Is the measured hot application path organized for reasonable locality and low avoidable overhead without speculative optimization?

## Required evidence

- Identify the actual hot path for patch solves, residual updates, macro-local sweeps, restriction, outer assembly, and V-cycle use.
- Record traversal order, memory layout, indirection depth, temporary allocation/copying, gather/scatter patterns, factor/matrix reuse, and cache invalidation.
- Run the common protocol's exact Release workloads P1/P2 and, at M9, P3 with the stated hardware, warmups, repetitions, normalization, baselines, memory metrics, and material-regression rule.
- Provide reproducible profiles, counters, timings, or allocation measurements tied to live application access patterns and exact builds.
- Compare setup and apply separately; normalize by sweeps/cycles, active patches/macros, and DOFs.
- Quantify whether duplicate cached views improve the hot path enough to justify synchronization/state burden.
- Check forward/reverse/palindromic and RAS traversal for material layout effects.
- Support D/N lifecycle checks with peak/retained memory, NaN/resource tooling where available, and cache invalidation observations.
- Activate a disposable cache-disable, per-patch allocation/copy, or locality-degrading traversal challenge; prove activation and show that the selected measurements distinguish it from the unmutated control when the induced effect meets the material threshold.

## Exclusions

- Do not require optimization based only on intuition.
- Do not fail debug-only audit/export paths as production hot-path regressions unless enabled in production.
- Do not remove caches solely for aesthetic simplicity.
- Do not use unrelated microbenchmarks.
- Do not compare different hardware/toolchains as if results were paired.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires valid fixed-workload evidence, sound cache reuse/invalidation, proportionate temporaries/indirection, no material regression by the common criterion, and a valid sensitivity control. A measured material avoidable hot-path regression, unsafe/stale cache, or unbounded per-apply overhead is **FAIL**. Noncomparable environments, missing baselines/measurements, or inactive sensitivity controls is **INCOMPLETE**; intuition alone cannot produce FAIL.

## Lane-specific report schema

After the common schema, include: hot-path map; exact build/hardware/workload manifest; baseline and candidate medians/MADs; normalized timing; allocation/copy and peak-memory results; locality/indirection analysis; cache benefit/cost/invalidation assessment; sensitivity control; limitations; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
