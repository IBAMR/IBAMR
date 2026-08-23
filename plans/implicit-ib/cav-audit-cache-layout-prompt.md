# Audit Lane C: Cache and Data-Layout Efficiency

**Prompt version:** `cav-audit-v9/C` (2026-08-21)

Use with `cav-independent-audit-prompt.md`. This lane is independent and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

Use the common protocol and implementation-plan Section 7 for the canonical focused-versus-broad `attest` evidence contract; do not redefine or substitute those scopes in this lane.

## Question

Is the production application hot path organized for reasonable locality and low avoidable overhead without speculative optimization or repeated representation translation?

## Required evidence

- Statically map the actual production hot path for patch solves, residual updates, macro-local sweeps, restriction, outer assembly, and V-cycle use.
- Record traversal order, memory layout, global/local index gathers, indirection depth, temporary allocation/copying, gather/scatter patterns, factor/matrix reuse, residual recomputation, ownership, and cache invalidation. Distinguish setup-time construction from every apply-time operation.
- Use the common protocol's exact W1/W2 and, at M9, W3 structural scenarios as reproducible trace carriers through an existing live production executable/case. Do not manufacture a benchmark harness.
- Apple Developer Tools are available. Record `xcrun xctrace list templates`, then independently verify any selected template against the exact candidate. Command-line Instruments tooling may provide allocation, retained-memory, copy, call-structure, repeated-work, cache/memory-access, or other non-timing counter/trace evidence. Time Profiler or `sample` may locate call structure only. Label every Instruments observation as non-timing structural evidence and ignore all elapsed-time, duration, percentage-time, throughput, and relative-speed fields contained in a trace. Compiler optimization records/remarks may supplement source inspection.
- Retain reproducible counter, allocation/copy, compiler-diagnostic, call-path, or equivalent direct observations tied to the live application path and exact build. Record exact executable, input/selector, tool/template/options, environment, exclusions, and raw trace hashes.
- Use a fully optimized production-equivalent candidate and linked libraries for tool-backed inspection when supported, while ignoring elapsed-time fields and making no performance comparison.
- Compare setup and apply structurally. Relate counters, allocations, copies, translations, indirections, and recomputations to sweeps/cycles, active patches/macros, and DOFs without converting them into speedup or throughput claims.
- Quantify whether duplicate cached views improve the hot path enough to justify synchronization/state burden.
- Check forward/reverse/palindromic and RAS traversal for material layout effects.
- Support D/N lifecycle checks with peak/retained memory, NaN/resource tooling where available, and cache invalidation observations.
- Analyze supported legal hierarchy geometries only. Flag any padding, splitting, synthetic-level, alternate-initialization, or translation path introduced solely for sub-minimum hierarchy patches as avoidable application/setup complexity; do not manufacture a sub-minimum profiling case. This does not prohibit naturally small algebraic CAV patches.
- Activate a disposable cache-disable, per-patch allocation/copy, or locality-degrading traversal challenge; prove activation and show that static review, allocation tracing, or cache/profile counters distinguish it from the unmutated control. No wall-clock materiality threshold is required for this sensitivity exercise.
- Verify that comments identifying hot paths, cached representations, copying/allocation constraints, traversal, and invalidation accurately match measured production behavior.

## Exclusions

- Do not require optimization based only on intuition.
- Do not fail debug-only audit/export paths as production hot-path regressions unless enabled in production.
- Do not remove caches solely for aesthetic simplicity.
- Do not use unrelated microbenchmarks.
- Do not run or require production timing; report speedup, throughput, medians/MADs, or percentage timing regressions; infer cache behavior from elapsed time; or require a measured speedup before reporting a structurally proven avoidable copy/allocation/translation.
- Do not compare different hardware/toolchains as if results were paired.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires a complete static hot-path/layout map, appropriate diagnostic counter/allocation/compiler/call-path evidence when useful, sound cache reuse/invalidation, proportionate temporaries/indirection, no structurally demonstrated avoidable hot-path defect, and a valid sensitivity control. Hardware counters being unavailable is a recorded limitation, not automatically `INCOMPLETE`, when complete static and other direct evidence answer the lane question; any claim specifically about cache misses still requires counters or cache simulation. A structurally proven repeated translation/copy/allocation, unsafe/stale cache, or unbounded per-apply overhead is **FAIL**. A missing hot-path map, raw timing substituted for structural evidence, unsupported cache-miss claim, unusable evidence provenance, or inactive sensitivity control is **INCOMPLETE**; missing timing is never a gap, and intuition alone cannot produce `FAIL`.

## Lane-specific report schema

After the common schema, include: static hot-path and representation map; exact build/hardware/structural-scenario manifest; tool/template/options and raw-trace hashes; counter/cache-simulation evidence when available; allocation/copy evidence; locality/indirection/translation analysis; cache benefit/cost/invalidation assessment; concrete improvement candidates with structural justification; sensitivity control; limitations; an explicit statement that no raw timing informed the outcome; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
