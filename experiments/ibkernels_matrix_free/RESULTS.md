# Native matrix-free comparison — 2026-09-15

## Assessment

The generic C++ side-coupling implementation is competitive with the existing
Fortran routines on this Apple M1 Max. Most measured direct-loop medians favor
C++, without fast-math. Two comparisons are effectively tied: ordered small
2D BS6 spreading (factor 1.003) and large shuffled 3D BS6 gathering (0.999).
There is no reproducible C++ slowdown in these cases. These results cover three
matching mathematical forms and six patch/marker configurations; they do not
establish whole-solver performance or performance on other architectures.

The rebuilt Autoibamr SAMRAI passed configuration, strict C++20 compilation and
real 2D/3D execution with Apple Clang. A separate IBSAMRAI2 rebuild was not needed
for this experiment. This does not qualify that SAMRAI release for every C++20
compiler or every IBAMR component.

## Source and environment

- Experiment branch: `codex/ibkernels-matrix-free-20260915`.
- Qualified R01d base: `c37b111fdabf267e83229c069698eca1a152bf1e`.
- Measured source: `7367b53e17033823ff820de4d21ba7ff2c82ccae`.
- Hardware: Apple M1 Max, 10 CPU cores, 64 GiB RAM; macOS 15.7.7.
- Apple developer tools: Xcode 26.3 (17C529), Apple Clang 17.0.0
  (`clang-1700.6.4.2`), compiler and macOS SDK selected explicitly with `xcrun`.
- Fortran: GNU Fortran 15.2.0 (Homebrew GCC 15.2.0_1).
- C, C++ and Fortran: `-O3 -mcpu=native -fno-fast-math`; C++20 and `NDEBUG`
  for the C++ Release targets. Native resolves to Apple M1. Applicable strict
  warnings remain enabled. No LTO or extra inlining variant was used.
- Dependencies: the user's completed Autoibamr optimized rebuild, loaded from
  `/Users/boyceg/code/autoibamr/opt/configuration/enable.sh`. SAMRAI 2025.10.29,
  PETSc 3.23.3, HDF5 1.12.2, libMesh 1.7.8 (opt), Boost 1.91.0 and Silo 4.11.
  SAMRAI's recorded C/C++/Fortran flags include `-O3 -march=native`, with debug
  macros disabled. PETSc records Apple Clang 17, GNU Fortran 15.2 and native flags.
- The native build is persistent at `build/ibkernels-matrix-free/Release-native`.
  Shared ccache activity was verified: 304 cacheable calls, 13 hits, 291 misses.
  Hashes and modification times of all 54 recorded dependency files remained
  unchanged through the measurements.

The linker reported inherited duplicate-library warnings and a macOS 15.0
link target against METIS/ParMETIS built for 15.7. Execution occurred on 15.7.7.
There were no introduced compiler warnings or build errors.

## Correctness

Native attest discovered and passed both dimensional fixtures, then passed both
again: **2/2 + 2/2**. The focused run took 1.32 s and its rerun 0.45 s in total.
Earlier Debug validation also passed 2/2 plus a repeat, and both dimensional
Debug benchmark smoke checks passed. Debug timings are not performance evidence.

The compact fixtures cover IB4, IB5, BS3, BS6, both normal/tangential 3/2 and 2/3
compositions, and a compiled application-defined cosine kernel. They exercise
real SAMRAI side arrays, nonzero patch indices, anisotropic spacing, ghost
access, exact placement ties, overlapping markers, additive spreading, selected
and repeated indices, shifts, clipping and preserved inputs.

Across the native double-precision cases, independent gather errors were at
most 4.45e-16, volume-scaled spread errors at most 1.43e-16, and adjoint errors
at most 2.49e-14. The float-coefficient indexed gather error was at most 2.90e-8.
All seven benchmark invocations additionally passed complete output comparisons
for all eight timed paths before measurement.

The inherited 3D IB5 spread routine omits its inner x-index update. Its scaled
error is reported as 0.2632433047652; Fortran remains unchanged, and that defective
spread is not timed. IB5's scalar Fortran delta supplies its reference values,
while placement is independently constructed. Other reference formulas are
independent of the evaluator implementation.

## Direct numerical loops

Each entry below is **Fortran median time / C++ median time**. A factor above
1 favors C++; 1 is a tie. Gather and spread are reported separately. CV is the
sample standard deviation divided by the mean; the last column is the largest
CV among all 24 operation/kernel series in that case, including wrappers.

| Dimension / case | IB4 gather / spread | BS3 gather / spread | BS6 gather / spread | Largest CV |
| --- | ---: | ---: | ---: | ---: |
| 2D small, ordered | 1.33 / 1.24 | 1.64 / 1.70 | 1.21 / 1.00 | 1.20% |
| 2D small, shuffled | 1.33 / 1.24 | 1.96 / 1.96 | 1.48 / 1.39 | 1.62% |
| 2D large, shuffled | 1.24 / 1.17 | 1.82 / 1.80 | 1.45 / 1.40 | 2.93% |
| 3D small, ordered (initial) | 1.45 / 1.05 | 1.41 / 1.64 | 1.07 / 1.27 | 90.49% |
| 3D small, ordered (repeat) | 1.45 / 1.07 | 1.41 / 1.64 | 1.07 / 1.27 | 1.94% |
| 3D small, shuffled | 1.45 / 1.07 | 1.40 / 1.62 | 1.07 / 1.40 | 3.01% |
| 3D large, shuffled | 1.28 / 1.14 | 1.09 / 1.14 | 1.00 / 1.25 | 4.72% |

The first six invocations produced 1,296 retained timing samples. Their median
series CV was 0.85%; five series in the initial small ordered 3D BS3 case had
CVs from 6.20% to 90.49%. Repeating that entire case added 216 samples, reduced
its largest CV to 1.94%, and reproduced the median ratios closely. Both runs
are shown and retained; no individual samples were removed.

No experiment builds or compiler-inspection jobs ran during timing. This was
an interactive host with background indexing: before/after snapshots showed
roughly 70–80% aggregate CPU idle. It was not an otherwise idle dedicated
benchmark host. The repeat and generally small variation support the larger
observed gaps; near-ties and small differences should be treated conservatively.

The larger 3D working set narrows the BS3 advantage and removes the small BS6
gather advantage. The data establish that trend, but do not isolate its cause.

## Patch entry-point costs

These factors compare `LEInteractor` against `cpp_patch`, across the recorded
cases and repeat:

| Dimension / kernel | Gather factor range | Spread factor range |
| --- | ---: | ---: |
| 2D IB_4 | 1.37–1.60 | 1.36–1.57 |
| 2D BSPLINE_3 | 1.80–2.43 | 1.92–2.31 |
| 2D BSPLINE_6 | 1.25–1.58 | 1.03–1.52 |
| 3D IB_4 | 1.36–1.48 | 1.10–1.23 |
| 3D BSPLINE_3 | 1.18–1.50 | 1.25–1.71 |
| 3D BSPLINE_6 | 1.04–1.09 | 1.26–1.40 |

These are costs of the respective entry points. `LEInteractor` selects indices
and copies components; the C++ patch path receives selected indices, constructs
the geometry object, and accesses interleaved values directly. Thus these factors
are not a claim for a complete drop-in replacement including marker selection.

The direct C++ and Fortran paths use identical scalar component buffers,
positions, zero-shift buffers and indices. C++ caches strides and inverse cell
volume in its prepared geometry; the existing Fortran routine retains its fixed
prologue. Neither path precomputes interpolation coefficients. The C++ patch
path uses interleaved values and an empty shift span, so subtracting its timing
from `cpp_loop` would not isolate geometry-construction cost.

## Compiler inspection and next steps

The final native assembly has out-of-line BS6 evaluator calls in both dimensions,
and out-of-line IB4 evaluator calls in 3D. Owning return values and visible
concrete types have not eliminated every evaluator call.

Clang reports a vectorized stencil-row loop with width 2 and interleave count 4.
However, the inspected 3D IB4 gather branch requires eight entries to enter that
vector path, while an IB4 row has at most four. That remark alone is not evidence
of SIMD execution for the row. GNU Fortran also reports vectorized loops,
including BS6 spreading. Assembly and remarks are retained for both languages.

A useful next experiment would separate fully interior stencils from clipped
stencils so that the compiler sees fixed row bounds, then test a targeted
inlining change. Each should be measured independently against this baseline.
Whole-hierarchy selection, ghost synchronization and concurrent spreading are
additional integration work before any production replacement decision.

## Reproduction and evidence

Run the commands in [README.md](README.md) to configure and build. From the
native build root, use the source `attest -N -R '^matrix_free/'`, followed by
`attest --keep-work-directories -R '^matrix_free/'`. Benchmark arguments are:

| Dimension | Cells per side | Markers | Iterations | Repeats | Shuffled |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 | 32 | 1024 | 200 | 9 | 0 |
| 2 | 32 | 1024 | 200 | 9 | 1 |
| 2 | 512 | 65536 | 5 | 9 | 1 |
| 3 | 16 | 4096 | 30 | 9 | 0 |
| 3 | 16 | 4096 | 30 | 9 | 1 |
| 3 | 128 | 65536 | 3 | 9 | 1 |

The repeated case uses the same arguments as small ordered 3D. Set
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1` and invoke
`build/ibkernels-matrix-free/Release-native/tests/benchmark-matrix-free-<dimension>d`.
Three warmups precede nine samples; path order rotates by repetition. Reset,
validation and consumed checksums are outside the timer. Spread accumulates the
same number of calls on each path. Marker positions are deterministic cell-order
positions or a fixed shuffle of the same data, not a random spatial cloud.

Evidence is retained in the experiment checkout under
`evidence/experiments/ibkernels-matrix-free-20260915/`:

- `experiment-source.patch`, `source-native.json`: exact source diff and hash.
- `environment-native.json`, `native-binaries.json`,
  `dependency-stability-native.json`: toolchain, flags, linkage and file hashes.
- `configure-release-native.log`, `build-release-native.log`, `indent-native.log`,
  `ccache-native-statistics.txt`: build and required clang-format 16.0.6 evidence.
- `native-discovery.log`, `native-focused.log`, `native-rerun.log`,
  `attest-native-results.json`, `attest-work/`: executable correctness evidence.
- `native-run-plan.json`, `native-benchmarks/runs.json`,
  `native-benchmarks/rerun.json`: exact commands, timestamps and case definitions.
- `native-benchmarks/[23]d-*.csv` and matching `.stderr`: all raw timing samples;
  `summary.csv` is the original six-case summary, `summary-rerun.csv` the repeat.
  Each raw run has before/after machine-load snapshots in its matching directory.
- `native-{cpp,fortran}-{2,3}d.s`, matching command/diagnostic files and
  `native-optimization-summary.json`: final native compiler inspection.

To summarize any explicit set of raw CSV files, use
`python3 experiments/ibkernels_matrix_free/summarize.py <raw CSV paths>`.
Do not include derived summary CSVs as inputs. Logs, binaries and raw evidence
remain outside the source commits. Changes are local to this experiment;
there is no hosted-CI, CAV integration, publication or merge qualification.
