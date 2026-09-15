# Matrix-free kernel improvements — 2026-09-15

## Assessment

The generic C++ implementation improves broadly after exposing complete stencil
bounds and simplifying row indexing. Paired native median-time ratios reach
**1.78×** relative to the original C++ loops. Gains remain in the larger 3D
working set. Wider-kernel gathering and several 2D spreading cases are close to
ties. The clear regression is **small-case 2D IB5 spreading, about 7% slower**;
its repeat confirms the result. The large 2D IB5 spread difference is about 3%.

The corrected Fortran implementation is still faster for 2D IB5 and for IB5
spreading in the large 3D case. The optimization is therefore retained as an
opt-in generic experiment with this measured tradeoff, not a universal
performance claim. All requested BS/CBS forms pass numerical checks, as does
the separately corrected 3D Fortran IB5 spread. The tables below report both
C++-before/C++-after and Fortran/C++ comparisons for every requested form.

## Changes and source revisions

The implementation remains source-private and generic over the owning-output
Cartesian evaluator API. When a stencil fits inside allocated side data,
including allocated ghosts, its loop bounds are compile-time constants. A
separate instantiation retains the original clipped bounds. Row addresses use
an invariant base offset plus the x-index. Stencil placement, accumulation
order, cell-volume scaling and owning coefficient storage are unchanged.
There are no kernel-name branches, compiler-specific attributes, fast-math,
coordinate copies or per-weight callbacks in the numerical loops.

- Expanded baseline: `92c0a79009098c618afd47a19a835d2ae581c249`.
  This includes all requested kernels and the corrected Fortran IB5 routine,
  with the original C++ numerical loops.
- First fixed-bound candidate: `b139e1341f0b18f98a5499db9f9c10c95a14c003`.
- Final measured implementation: `1ed9b79dd33098045a54190af1ef3158940d5241`.
- Qualified evaluator base: `c37b111fdabf267e83229c069698eca1a152bf1e`.
- Initial three-kernel report and results remain preserved at
  [RESULTS.md](RESULTS.md), from source `7367b53e17033823ff820de4d21ba7ff2c82ccae`.

The first candidate improved most pilot cases but slowed small shuffled 3D BS5
gathering by about 12%. Changing the row-index expression removed that regression
in the focused probe. Both preliminary measurements are retained separately
from the final six-case comparison. Their Fortran libraries already contain the
IB5 correction; defective and corrected IB5 timings are never compared as if
they performed the same operation.

## IB5 correction and correctness coverage

`lagrangian_ib_5_spread3d` omitted `ic0 = ic_lower(0) + i0` from its innermost
loop. The one-line correction restores spreading to each x-index. Its previous
volume-scaled whole-field discrepancy was 0.2632433047652; the corrected 3D
fixture reports **2.602085213965e-17** against the C++ spread result. Independent
placement and scalar-delta reference checks also pass. Corrected IB5 is timed
in every final benchmark configuration.

Native attest discovers two fixtures, one per dimension. Each now contains 37
labeled result records: 18 forms with full and clipped stencils, plus an indexed
float-coefficient case. The forms are IB4, IB5, BS2-6, CBS(k+1)k and CBSk(k+1) for
k=1,...,5, and an application-defined cosine kernel. **CBSnm means normal width
n and tangential width m**; both orientations are explicit in fixture output.
The BS1 factor selects the upper nearest-grid index at an exact tie.

The final source passes **2/2 plus a 2/2 rerun in both Debug and native Release**.
The expanded baseline passed discovery, focused and rerun checks in both
configurations. Final fixture values match the reviewed expanded-baseline
fixtures. Strict warnings and the actual `make indent` target with clang-format
16.0.6 pass. These are local focused checks, not hosted CI or a full IBAMR suite.

Across the full and clipped double-precision cases, the maximum independent
errors are 8.882e-16 for gathering and 2.221e-16 for volume-scaled spreading;
the largest adjoint error is 3.908e-14. The indexed float-coefficient gathering
error is at most 2.891e-8. Cases use real SAMRAI side arrays, nonzero patch
indices, anisotropic spacing, side-centering ties, ghosts, overlapping markers,
additive output, selected/repeated indices, shifts and preserved inputs.
B-spline reference values use an independent truncated-power formula.

CBS12 has **no matching existing Fortran backend**. It has independent
correctness checks and C++ direct/patch timings; Fortran comparison is unavailable.
CBS21 uses the existing discontinuous-linear Fortran backend. Clipped fixtures
validate C++ against the independent reference and adjoint identity; existing
Fortran wrappers are not called outside their ghost-width preconditions.

## Measurement method and environment

Both implementations use the user's rebuilt optimized dependencies from
`/Users/boyceg/code/autoibamr/opt/configuration/enable.sh`. The installed
IBSAMRAI2-2025.10.29 works for these Apple Clang C++20 builds and executions;
no additional dependency rebuild or installation was needed.

The hardware and flags match the initial comparison: Apple M1 Max, 10 CPU cores,
64 GiB RAM, macOS 15.7.7; Xcode 26.3 (17C529), Apple Clang 17.0.0 selected by
`xcrun` with its SDK, and GNU Fortran 15.2.0. C, C++ and Fortran use
`-O3 -mcpu=native -fno-fast-math`. C++ uses C++20 and NDEBUG. No LTO or extra
inlining flag is used. Native resolves to Apple M1. Inherited duplicate-library
and macOS deployment-target linker warnings remain as recorded in the initial
report; no compiler warnings or build errors were introduced.

| Dimension | Cells per side | Markers | Iterations per sample | Samples | Order |
| --- | ---: | ---: | ---: | ---: | --- |
| 2 | 32 | 1,024 | 200 | 9 | ordered |
| 2 | 32 | 1,024 | 200 | 9 | shuffled |
| 2 | 512 | 65,536 | 5 | 9 | shuffled |
| 3 | 16 | 4,096 | 30 | 9 | ordered |
| 3 | 16 | 4,096 | 30 | 9 | shuffled |
| 3 | 128 | 65,536 | 3 | 9 | shuffled |

Each executable checks complete outputs before timing. Three warmups precede
nine samples, operation order rotates across samples, and baseline/final program
order alternates across cases. Reset and consumed checksums are outside timing.
Spreading accumulates equally many calls on each path. There is one MPI rank
and one CPU thread; OMP, OpenBLAS and Accelerate thread limits are set to one.
No builds or assembly-generation jobs run during timing. Before/after load
snapshots accompany every invocation on this interactive desktop host.

Direct C++ and Fortran calls use the same scalar component buffers, selected
indices, positions and explicit zero-shift buffers. C++ caches patch geometry,
strides and inverse cell volume; the existing Fortran prologue remains timed.
Both evaluate coefficients during each call. Patch timings describe different
entry points: C++ constructs its geometry object and accesses interleaved values
with an empty shift span; LEInteractor also selects indices and copies
components. Subtracting direct from patch times does not isolate geometry cost,
and these comparisons do not establish whole-solver or production-integration
speedups.

## Original C++ / optimized C++

Each cell gives the range of median-time ratios across the three configurations
and the retained repeats for that dimension. A ratio above 1 favors optimized
C++; values near 1 should be treated as ties. These are direct numerical calls.

| Kernel | 2D gather | 2D spread | 3D gather | 3D spread |
| --- | ---: | ---: | ---: | ---: |
| IB4 | 1.42–1.47 | 1.64–1.75 | 1.30–1.32 | 1.29–1.30 |
| IB5 | 1.14–1.17 | 0.93–0.97 | 1.07–1.09 | 1.27–1.45 |
| BS2 | 1.39–1.67 | 1.42–1.72 | 1.28–1.62 | 1.28–1.68 |
| BS3 | 1.42–1.48 | 1.47–1.54 | 1.38–1.49 | 1.33–1.42 |
| BS4 | 1.35–1.41 | 1.46–1.57 | 1.09–1.21 | 1.29–1.35 |
| BS5 | 1.06–1.17 | 0.99–1.10 | 1.04–1.06 | 1.28–1.46 |
| BS6 | 1.20–1.26 | 1.00–1.05 | 1.00–1.03 | 1.42–1.68 |
| CBS21 | 1.42–1.59 | 1.22–1.57 | 1.08–1.45 | 1.06–1.41 |
| CBS32 | 1.49–1.67 | 1.52–1.74 | 1.43–1.67 | 1.44–1.78 |
| CBS43 | 1.32–1.37 | 1.50–1.61 | 1.28–1.31 | 1.06–1.26 |
| CBS54 | 1.32–1.41 | 1.26–1.48 | 1.10–1.15 | 1.38–1.44 |
| CBS65 | 1.11–1.15 | 1.03–1.05 | 1.01–1.03 | 1.33–1.57 |
| CBS12 | 1.44–1.61 | 1.24–1.59 | 1.24–1.62 | 1.29–1.56 |
| CBS23 | 1.47–1.62 | 1.50–1.66 | 1.41–1.64 | 1.44–1.59 |
| CBS34 | 1.37–1.37 | 1.52–1.62 | 1.30–1.30 | 1.12–1.32 |
| CBS45 | 1.31–1.41 | 1.27–1.46 | 1.11–1.11 | 1.32–1.45 |
| CBS56 | 1.11–1.17 | 1.00–1.01 | 1.00–1.02 | 1.36–1.62 |

## Fortran / optimized C++

Each cell gives the range of median-time ratios across the three configurations
and the retained repeats for that dimension. A ratio above 1 favors optimized
C++; values near 1 should be treated as ties. These are direct numerical calls.

| Kernel | 2D gather | 2D spread | 3D gather | 3D spread |
| --- | ---: | ---: | ---: | ---: |
| IB4 | 1.77–1.86 | 2.03–2.12 | 1.70–1.86 | 1.38–1.47 |
| IB5 | 0.96–0.96 | 0.89–0.92 | 1.03–1.19 | 0.90–1.03 |
| BS2 | 1.72–5.29 | 1.62–4.96 | 1.61–3.96 | 1.72–4.32 |
| BS3 | 2.41–2.90 | 2.59–3.02 | 1.58–2.08 | 1.71–2.18 |
| BS4 | 1.99–3.29 | 2.07–3.27 | 1.21–1.40 | 1.59–2.14 |
| BS5 | 1.76–1.79 | 1.94–2.11 | 1.21–1.40 | 1.46–1.97 |
| BS6 | 1.50–1.83 | 1.03–1.47 | 1.05–1.08 | 1.73–2.33 |
| CBS21 | 1.68–4.25 | 1.38–3.80 | 1.18–3.19 | 1.04–2.86 |
| CBS32 | 3.37–6.40 | 3.83–6.32 | 2.88–6.28 | 3.26–6.83 |
| CBS43 | 2.41–3.81 | 2.42–3.68 | 1.88–2.44 | 1.75–2.31 |
| CBS54 | 2.54–3.38 | 2.96–3.89 | 1.81–2.10 | 2.24–2.82 |
| CBS65 | 1.71–2.03 | 1.13–1.60 | 1.47–1.61 | 2.07–2.70 |
| CBS12 | unavailable | unavailable | unavailable | unavailable |
| CBS23 | 3.20–5.57 | 3.75–5.88 | 2.02–3.46 | 2.24–3.59 |
| CBS34 | 2.38–4.32 | 2.43–4.16 | 1.65–2.09 | 1.64–2.27 |
| CBS45 | 2.67–3.06 | 3.05–3.52 | 1.57–1.82 | 1.79–2.29 |
| CBS56 | 1.70–2.19 | 1.08–1.68 | 1.23–1.29 | 1.91–2.55 |

## Variation and repeat checks

The six paired configurations retain **14,256 raw samples**. Repeating both
small 2D configurations with reversed baseline/final order adds **4,752**,
for **19,008 final-comparison samples**. The tables include both original and
repeated median ratios; no samples are discarded or pooled across invocations.
Each executable emits 1,188 samples: eight paths for each of 16 matching forms,
and four C++ paths for CBS12, with nine samples per path.

CV is sample standard deviation divided by mean. The median CV across all
2,112 operation series is **0.37%**. The ordered 2D repeat reduces its maximum
CV from 12.59% to 3.74%. The shuffled repeat still has outliers, including a
65.19% CV in the baseline CBS32 LEInteractor spread series; its baseline C++
spread series has a 29.26% CV. The largest CV in the original matrix is 25.77%.
The repeat does not make every series stable. Across repeated direct-loop
ratios, the median change is 0.38% and the largest change is 4.32%. Larger gaps
are supported by these checks; near-ties and precise wrapper ratios remain
uncertain on this interactive host.

| Configuration | Original largest CV | Repeat largest CV |
| --- | ---: | ---: |
| 2D small ordered | 12.59% | 3.74% |
| 2D small shuffled | 25.77% | 65.19% |
| 2D large shuffled | 5.49% | — |
| 3D small ordered | 8.23% | — |
| 3D small shuffled | 10.50% | — |
| 3D large shuffled | 12.86% | — |

Patch-entry C++ before/after ratios span 0.96–1.79; the small 2D IB5 spread
regression also appears there. Raw wrapper medians and variation are in the
summaries. Fortran control before/after ratios span 0.92–1.13, which is another
reason to avoid interpreting small differences as precise speedups.
Before/after snapshots for the main matrix show 81–96% aggregate CPU idle.
All 54 recorded dependency files, both corrected IBTK libraries, the baseline
executables and the final executables retain their expected hashes through
measurement.

## Compiler observations and limits

The final assembly and optimization remarks confirm that Clang sees constant
bounds and unrolls the complete stencil loops. For example, the inspected 3D
BS5 axis-zero gather changes from a 25-term body with a five-plane loop in the
first candidate to a fully unrolled 125-term sequence with paired field loads.
The arithmetic still uses one ordered accumulation. This is static evidence
of changed generated code; it does not isolate which instructions account for
the measured difference.

Large evaluator calls remain out of line, including 3D IB4 and IB5 and BS5/BS6
in both dimensions. The current change does not eliminate all evaluator overhead.
The benchmark exercises complete stencils; clipping is correctness-tested but
its performance is not characterized. No conclusion here covers other CPU
architectures, whole-hierarchy marker selection, ghost synchronization or
concurrent spreading. The experiment remains opt-in; production coupling
selection is unchanged.

## Reproduction and evidence

Use [README.md](README.md) for the configuration, strict build and native attest
commands. The native build is `build/ibkernels-matrix-free/Release-native`.
The final timing plan uses the table above, invoking
`tests/benchmark-matrix-free-<dimension>d` with cells, markers, iterations,
samples and shuffle flag. The optional sixth argument selects a kernel.

All evidence is retained under
`evidence/experiments/ibkernels-matrix-free-20260915/improvement/`:

- `baseline.json`, `baseline-binaries/`: expanded baseline revision, executable
  and corrected IBTK library copies and hashes.
- `fixed-bounds.json`, `fixed-bounds-binaries/`, `fixed-bounds.patch`, `pilot/`:
  first candidate and its paired small-shuffled comparison.
- `row-offset.json`, `row-offset.patch`, `row-probe/`: focused indexing probe.
- `final-source.json`, `final-optimization.patch`: final measured source and
  binary hashes. Baseline and final executables use identical corrected IBTK
  libraries; their hashes are checked before and after measurement.
- `coverage-attest.json`, `coverage-*.log`, `coverage-*-*d/output`,
  `fixed-bounds-attest.json`, `final-attest.json`, `final-*.log`,
  `build-*.log`, `indent-*.log`: baseline and candidate build, formatting and
  correctness records.
- `final-plan.json`, `final-matrix/runs.json`, `final-matrix/[23]d-*.csv`,
  matching stderr and per-invocation load snapshots: complete paired matrix.
- `repeat-plan.json`, `repeat/runs.json`, `repeat/[23]d-*.csv`: retained full-case
  repeats of the two noisier 2D small cases, with reversed variant order.
- Each measurement directory's `summary.csv` gives sample counts, median,
  minimum, maximum, CV, baseline/final ratios and matching Fortran/C++ ratios.
  `comparison-analysis.json` records the report's aggregates and exceptions.
- `final-{2,3}d.s`, `final-assembly-command-{2,3}d.json`,
  `final-{2,3}d-remarks.txt`, `final-assembly-excerpts-{2,3}d.json`: compiler
  inspection using the exact native benchmark compilation flags.
- `dependency-stability-final.json`, `final-custody.json`: dependency hashes and
  final local source/evidence custody.

`run-comparison.py` records exact sequential commands and load snapshots.
`summarize-comparison.py <measurement-directory>` summarizes only raw
`[23]d-*.csv` inputs, retaining every sample. All elapsed CSV values are raw
seconds for the recorded iteration count; summary times are seconds per call.
No dependency installation, hosted CI, publication, merge, CAV integration or
independent audit is claimed.
