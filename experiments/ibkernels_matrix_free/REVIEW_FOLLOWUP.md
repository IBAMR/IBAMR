# Independent-review follow-up

2026-09-20. Reviewed source: `d4287f92a8aee00f6e467b5aaae9a8343d9760cc`;
measured implementation: `ca0dc4b521827ae2dc8e84dd6e601eeac66d31f2`.
This follow-up corrects the study's interpretation using source and saved
evidence. It introduces no implementation changes or new performance results.
No builds, tests, benchmarks, Instruments captures or hardware queries were run.
Live measurement remains on hold because the system is in other use.

Evidence paths below are relative to
`evidence/experiments/ibkernels-matrix-free-20260915/` in this checkout.
The independent review and its scripts were preserved unchanged under
`independent-review/received-20260920/`, with hashes in `custody.json`.
A later revision appeared in the external review directory during this pass;
it is preserved separately in `independent-review/received-20260920-update/`.
That revision acknowledges user confirmation of the chip and adds tool comments.
The original review snapshot and its scripts remain unchanged.

## Conclusions retained and corrected

- The independent reviewer found no correctness defect in the coupling contract
  from source and saved outputs. This is independent review, not an independent
  test run. Historical local test passes are not hosted CI evidence.
- The reviewer recomputed the material timing ratios from all 9,828 samples.
  The numerical tables in [FACTORIZATION.md](FACTORIZATION.md) remain unchanged.
- The 2D CBS65 spread result is an unresolved contracted-path anomaly. For the
  small ordered case, contracted CBS65 takes 80.8 microseconds and contracted
  CBS56 takes 57.9; expanded times match across the two kernels, as do factorized
  times. Mirrored stencils and nearly identical assembly make simple arithmetic
  count an inadequate explanation. Pitch/alignment and SIMD update behavior
  remain hypotheses. FACTORIZED stays as a diagnostic comparison pending
  resolution, without a permanent kernel-specific dispatch rule.
- C++/Fortran ratios compare Apple Clang with GNU Fortran and different loop
  structures. They do not isolate algorithm or language effects.
- The large 3D workload is shuffled over four x-cell bands. An ordered version
  of the same large workload is missing; additional spatial distributions will
  be needed before generalizing to other marker sets.

## Cache footprint and the LRU model

The original 8.548 MiB CBS65 figure counts unique double entries across all
components. The review correctly identifies cache-line rounding as important.
A separate source-derived calculation reproduces its aggregate line counts,
assuming each component array starts on a cache-line boundary:

| Kernel | Assumed line size | Axis 0 MiB | Axis 1 MiB | Axis 2 MiB | Sum MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| CBS65 | 64 B | 6.961 | 4.286 | 4.286 | 15.533 |
| CBS65 | 128 B | 11.215 | 8.572 | 8.572 | 28.359 |
| CBS56 | 64 B | 7.155 | 4.348 | 4.348 | 15.852 |
| CBS56 | 128 B | 11.536 | 8.697 | 8.697 | 28.930 |

These are unique field-line unions, not measured cache occupancy or misses.
The calculation excludes positions, shifts, indices, marker values, stack
storage, prefetch, associativity, cache history and interference. Scripts and
exact counts are in `independent-review/followup-20260920/line_footprint.py`
and `line-footprint.json`. This short static calculation was run offline;
the reviewer's full LRU simulation was not rerun.

**The reviewer's approximately 15x L2-miss reduction cannot yet be applied to
this implementation.** Its `footprint_independent.py` simulates all components
for marker 0, then all components for marker 1, and so on. Both
`tests/matrix_free/benchmark.cpp::component_operation` and
`tests/matrix_free/coupling-inl.h::couple` instead finish all markers for one
component before starting the next. That changes reuse distances substantially.
Its shuffled order also uses Python's random generator, rather than the exact
benchmark permutation.

Although the aggregate exceeds the assumed 12 MiB capacity, each component's
field-line union is smaller than 12 MiB in this calculation. This establishes
neither that the real working set fits nor that capacity misses dominate.
A revised model must reproduce component order, marker permutation, other
arrays and repeated-call history. Cache capacities remain assumptions until
verified. Deterministic binning is still a useful experiment; its benefit is
unmeasured and the 15x estimate is not a quantitative forecast for this code.

## What the saved CPU Counters trace establishes

The successful smoke capture is for the small 3D CBS65 case: 16 cubed cells,
4,096 ordered markers, approximately 0.330 MiB allocated field, 200 calls per
sample and three repeats. It records cycles and instruction bottleneck metrics,
not cache misses. It does not characterize the large-field cache behavior.

The reviewer's offline reconstruction suggests approximately 4.13 cycles per
stencil entry for expanded gather and 1.05 for contracted gather, with
factorized/contracted spread near 1.3. These are **approximate reconstructed
phase estimates**, not directly delimited per-operation counter measurements.
`align_phases.py` fits a start offset and constant inter-phase gap using the
variation of the Useful metric, then uses cumulative CSV durations. Its close
end alignment is useful corroboration but does not independently verify every
boundary; reset/checksum work and partial intervals add uncertainty.

The gather estimates are consistent with the saved serial FMA dependency chain
and shorter contracted reductions. They do not by themselves establish a sole
bottleneck. Similar spread costs also do not prove that multiplication is
irrelevant or that load/store updates are the limiting resource.

A thread-aware read of the export finds 25,361 main-thread cycle intervals.
Adjacent recorded main-thread intervals change CPU in 6,764 of 25,360 pairs
(26.67%). The main thread visits all eight recorded P cores, as well as E cores
in the whole-process trace. This supports caution about migration, but is not
a complete migration count or a measured cross-cache-cluster migration rate.
The review's 1.79 GHz fifth percentile is derived from cycles divided by
interval duration; it is not a separately recorded clock-frequency measurement
and may be affected by off-CPU time. Thread-aware counts are saved in
`independent-review/followup-20260920/counter-thread-summary.json`.

Hardware identity is better documented than the review states:
`environment-preparation.json` records a successful query on September 15 for
Apple M1 Max, 64 GiB and 10 physical/logical cores. The failed later query in
`environment-native.json` does not remove that evidence. Neither file verifies
cache sizes, line sizes or the CPU-to-shared-cache-cluster mapping.

## Regression coverage

The review understates the explicit assertions slightly: `compare_mode` checks
both differences against CONTRACTED and the EXPANDED/FACTORIZED adjoint errors
at 1e-11. Default CONTRACTED reference/adjoint and Fortran errors are printed
and compared through `numdiff -r 1e-6 -a 1e-10`. Float rounding errors are also
part of expected output. Explicit scale-aware error bounds would improve these
checks without requiring exact floating-point reproducibility.

Empty clipped stencils, partial/anisotropic ghost widths and scalar component
stride one remain useful native regression additions. Stride one already has
full-output checks in the benchmark; it is absent from the native regression.
The custom Cartesian evaluator's lack of a factor interface selects the
fallback at compile time. No new correctness defect was established here.

## Deferred experiments, in order

1. Add operation/axis selection and signposts around repeated calls, excluding
   initialization, validation, resets and checksums from analyzed intervals.
   Add an ordered large-field case using the exact same marker set as shuffled.
2. Resolve the 2D CBS65/56 contracted-spread anomaly with separate axes and
   controlled pitch/alignment changes, including neighboring grid sizes.
   Consider retiring FACTORIZED only after this comparison is understood.
3. Compare deterministic spatial binning with ordered and shuffled large-field
   traversal. Keep the existing banded set and add distributed marker sets.
   Account for binning cost, reuse across calls and changed spread summation
   order; repair the offline access model before using its predictions.
4. Try fixed-order SIMD or multiple-accumulator gather reductions. Deterministic
   reassociation is allowed; validate numerical error bounds after changes.
5. Use CPU Profiler and CPU Counters bottleneck/cache-miss modes supported by
   this hardware, including `l1d_miss_sampling` if available. Normalize cycles
   and events per processed stencil entry. Retain thread/core identity, verify
   cache topology, and flag migration and scheduling/frequency effects per
   interval. Confirm improvements with unprofiled timings once conditions permit.

This is a deferred plan, not authorization to resume measurements during the
hold. The outstanding performance causes remain hypotheses.
