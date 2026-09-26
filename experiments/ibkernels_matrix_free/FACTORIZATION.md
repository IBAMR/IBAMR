# Factorized CBS coupling: native performance results

Keeping one-dimensional factors and contracting the tensor product substantially
improves high-order CBS gathering. Against the preceding expanded C++ algorithm,
CBS56/65 gather is **3.9–4.4× faster on the small 3D field** and **2.5–2.7× faster
on the large shuffled 3D field**. Their large-field spreading improves about
**1.5×**. These are patch coupling measurements on one Apple M1 Max CPU thread.

Use **CONTRACTED** as the general experimental default. Retain **FACTORIZED**
for comparison while investigating an unresolved **2D CBS65 contracted-spread
slowdown**. Factorized CBS65 spreading improves the expanded path by
**1.60–1.63×**, versus **1.15–1.17×** for contraction, but the mirrored CBS56
case does not show the same slowdown. These timings do not establish an intrinsic
advantage for factorized CBS65 spreading. The lowest-order CBS12/21 cases are
essentially ties. No kernel-specific dispatch or machine-specific selection
threshold is introduced.

The [independent-review follow-up](REVIEW_FOLLOWUP.md) qualifies the cache model,
counter attribution and validation coverage. Live measurements remain on hold.

Measured implementation: `ca0dc4b521827ae2dc8e84dd6e601eeac66d31f2`.
The clean starting revision was `24f9ae65366dce550db25db1f30f85ed1c913d78`;
its expanded optimization was measured previously at
`1ed9b79dd33098045a54190af1ef3158940d5241` in [IMPROVEMENTS.md](IMPROVEMENTS.md).
This report compares all three modes within the same new executable, preserving
the old reports and binaries. CBSnm means normal width n and tangential width m.

## Gains over expanded C++

Each entry is a range of **expanded median time / contracted median time**
over the configurations below. Values greater than one mean faster contraction.
The repeated small shuffled 3D invocation is kept separate, not pooled into the
original invocation. Gather and spread are the direct component-buffer paths.

| Kernel | 2D gather | 2D spread | 3D gather | 3D spread |
| --- | ---: | ---: | ---: | ---: |
| CBS21 | 0.98–1.01 | 0.95–0.99 | 1.00–1.02 | 0.98–1.00 |
| CBS12 | 0.98–1.06 | 0.95–1.05 | 0.99–1.05 | 0.98–1.01 |
| CBS32 | 1.10–1.14 | 1.03–1.07 | 1.03–1.29 | 1.04–1.13 |
| CBS23 | 1.10–1.19 | 1.07–1.09 | 1.08–1.43 | 1.02–1.22 |
| CBS43 | 1.20–1.25 | 1.15–1.16 | 1.45–2.15 | 1.35–1.68 |
| CBS34 | 1.24–1.27 | 1.13–1.17 | 1.49–2.12 | 1.32–1.60 |
| CBS54 | 1.30–1.31 | 1.15–1.21 | 1.93–2.86 | 1.31–1.40 |
| CBS45 | 1.31–1.34 | 1.18–1.20 | 1.86–2.54 | 1.36–1.48 |
| CBS65 | 1.44–1.47 | 1.15–1.17 | 2.48–3.98 | 1.33–1.48 |
| CBS56 | 1.45–1.47 | 1.53–1.61 | 2.72–4.36 | 1.35–1.47 |

Retaining factors alone already helps, but the successive sums account for much
of the 3D gather improvement. The following table compares both approaches on
the repeated small shuffled 3D case (16 cells per side, 4,096 markers):

| Kernel | Factorized gather | Contracted gather | Factorized spread | Contracted spread |
| --- | ---: | ---: | ---: | ---: |
| CBS21 | 1.00× | 1.00× | 0.99× | 0.99× |
| CBS12 | 1.02× | 1.01× | 1.03× | 1.01× |
| CBS32 | 1.04× | 1.28× | 1.08× | 1.13× |
| CBS23 | 1.02× | 1.39× | 1.15× | 1.22× |
| CBS43 | 1.24× | 2.14× | 1.60× | 1.68× |
| CBS34 | 1.21× | 2.10× | 1.50× | 1.58× |
| CBS54 | 1.11× | 2.85× | 1.34× | 1.38× |
| CBS45 | 1.08× | 2.54× | 1.34× | 1.36× |
| CBS65 | 1.42× | 3.98× | 1.32× | 1.34× |
| CBS56 | 1.34× | 4.35× | 1.36× | 1.36× |

### Absolute times and the Fortran comparison

Median milliseconds per call, including all components. Small 3D uses the
separate repeat above; large 3D uses 128 cells per side and 65,536 shuffled markers.
Fortran uses the matching existing direct routine and identical component buffers.

| Case | Kernel / operation | Expanded | Factorized | Contracted | Fortran |
| --- | --- | ---: | ---: | ---: | ---: |
| Small 3D | CBS65 gather | 2.452 | 1.721 | 0.616 | 3.947 |
| Small 3D | CBS65 spread | 1.003 | 0.762 | 0.749 | 2.713 |
| Small 3D | CBS56 gather | 3.104 | 2.319 | 0.714 | 3.989 |
| Small 3D | CBS56 spread | 1.103 | 0.811 | 0.810 | 2.822 |
| Large 3D | CBS65 gather | 48.255 | 37.034 | 19.466 | 68.800 |
| Large 3D | CBS65 spread | 28.398 | 20.187 | 19.187 | 57.986 |
| Large 3D | CBS56 gather | 58.758 | 47.266 | 21.572 | 71.390 |
| Large 3D | CBS56 spread | 31.537 | 22.463 | 21.404 | 62.773 |

For these high-order 3D kernels, contracted C++ is 5.59–6.40× faster than direct
Fortran gathering on the repeated small case and 3.31–3.53× on the large case.
Large-case spreading is 2.93–3.02× faster than direct Fortran.
All nine supported CBS Fortran forms are retained in the raw comparison;
CBS12 has no existing matching Fortran backend.

These are implementation comparisons: C++ uses Apple Clang and Fortran uses
GNU Fortran, with differently structured loops and timed prologues. The ratios
do not isolate language or tensor-contraction effects. The expanded C++ path
provides the same-compiler baseline for evaluating the algorithm change.

C++ patch-entry timings show the same high-order 3D benefit: CBS56/65 gather
improves 3.95–4.35× on the small cases and 2.41–2.72× on the large case; large-case
spread improves 1.50–1.51×. The summary also retains LEInteractor timings.
Patch C++ and LEInteractor perform different setup and layout work, so their
ratio does not isolate numerical-loop or geometry-construction cost.

## Implementation and numerical behavior

`IBKernelEvaluatorTensorProduct::evaluateFactors<Axis, Coefficient>(r)` returns
an owning tuple of one-dimensional arrays in coordinate order. The normal
kernel supplies direction `Axis`; the tangential kernel supplies the others.
The expanded-output `evaluate` API remains available. The experimental
`SideCoupling` detects the optional factor interface at compile time; evaluators
without it use the Cartesian path.

Three compile-time modes are available:

- **EXPANDED:** construct the complete coefficient tensor, then apply it.
- **FACTORIZED:** retain one-dimensional arrays and form products at the point of use.
- **CONTRACTED:** gather using successive row and plane sums; spread by scaling
  the marker value once per plane and row. This is the experimental default.

For example, a 3D gather computes each x-row dot product, weights those row sums
in y, then weights the plane sums in z. CBS65 retains 16 coefficients instead of
150 per component; CBS56 retains 17 instead of 180. Every stencil field entry
still requires a load or update.

Contractions intentionally change floating-point association. The selected
coefficient precision applies to the one-dimensional factors; field arithmetic
uses double, including with float factors. Deterministic throughput is the
priority, and bitwise cross-implementation agreement is not a requirement.
Placement, clipping, selected/repeated indices, shifts, additive spreading and
cell-volume scaling retain their contracts. There is no coefficient cache across
calls, marker-coordinate copy, heap allocation in the component loop, parallel
reduction, or concurrent scatter.

### Generated code

Assembly from the actual benchmark translation unit, compiled with the measured
Apple flags, confirms that CBS65 3D axis-0 gathering uses unrolled x/y row
reductions and a five-plane loop. Its total stack frame drops from 1,504 bytes
for expanded weights to 496 bytes for factors and 448 bytes for contraction.
The evaluator remains an out-of-line call in these instantiations: owning return
values do not imply that evaluation is fully inlined.

Both factorized and contracted CBS65 2D axis-0 spread paths emit two-double SIMD
updates. Contraction uses fewer multiplications but is slower in the measured
all-component operation. Arithmetic count alone does not explain that result;
the mirrored CBS56 contracted paths have nearly identical assembly yet avoid
the large slowdown. In the small ordered 2D case, contracted CBS65 takes
80.8 microseconds versus 57.9 for contracted CBS56, while the two kernels have
similar expanded times and similar factorized times. Per-axis comparisons and
array-pitch/alignment variations are deferred discriminating checks. Exact
assembly, extracted functions and compiler vectorization remarks are retained
with the evidence; no microarchitectural cause has been established.

## Correctness

Debug and native Release each discover and pass both dimensional native attest
fixtures, followed by a passing 2/2 rerun. Each fixture contains 40 labeled
records covering BS2–6, both CBS(k+1)k and CBSk(k+1) for k=1,...,5, IB4/IB5,
and a custom Cartesian evaluator without a factor interface. Full and clipped
stencils compare all three modes. Additional cases cover float factors, indexed
and shifted markers, repeated spread indices, untouched outputs and independent
ownership of returned factors.

The maximum double-precision error against the independent reference is
2.665e-15 for gathering and 2.221e-16 for volume-scaled spreading. The maximum
adjoint discrepancy across modes is 3.908e-14. Float-factor indexed checks have
errors below 3e-8. The previous 3D Fortran IB5 correction remains intact.
Every benchmark invocation checks complete outputs before timing. Within every
kernel/operation series, all nine repetitions produce the same recorded checksum.

These maxima describe the saved validation runs. The regression explicitly
checks mode differences and the expanded/factorized adjoint discrepancies at
1e-11. Reference, Fortran and default contracted adjoint errors are printed and
checked through the normal output comparison (`numdiff -r 1e-6 -a 1e-10`).
Float-error values are also stored in expected output; explicit scale-aware
bounds would better express their intended tolerance. Empty clipped stencils,
partial ghost widths and scalar component stride one need native regression
cases. Stride one is already exercised by the benchmark's pre-timing checks.

## Measurement method and variation

Apple M1 Max, macOS 15.7.7, Xcode 26.3/Apple Clang 17 via `xcrun`, GNU Fortran
15.2, and the user's optimized Autoibamr dependencies. C, C++ and Fortran use
`-O3 -mcpu=native -fno-fast-math`; C++ uses C++20. All 54 recorded dependency
files and both IBTK libraries still match their pre-experiment hashes.
No dependency rebuild was needed. Strict-warning builds and `make indent`
with clang-format 16.0.6 pass; inherited linker warnings are unchanged.

The earlier `environment-preparation.json` records a successful hardware query
identifying Apple M1 Max, 64 GiB and 10 physical/logical cores. The later failed
query in `environment-native.json` does not invalidate that saved record.
Cache capacities, line sizes and the mapping of CPUs to shared-cache clusters
were not established by those queries.

Three warmups precede nine samples, rotating operation order. Reset and consumed
checksums are outside timing but can affect initial cache state. Each invocation
uses one MPI rank and single-threaded numerical operations; runtime helper
threads may also exist. OMP, OpenBLAS and Accelerate thread limits are one.
Builds and assembly inspection run outside the timing interval. The seven invocations retain
**9,828 raw timing samples**; all complete successfully.

| Dimension / case | Cells per side | Markers | Iterations per sample |
| --- | ---: | ---: | ---: |
| 2D small ordered | 32 | 1,024 | 200 |
| 2D small shuffled | 32 | 1,024 | 200 |
| 2D large shuffled | 512 | 65,536 | 5 |
| 3D small ordered | 16 | 4,096 | 30 |
| 3D small shuffled, plus separate repeat | 16 | 4,096 | 30 |
| 3D large shuffled | 128 | 65,536 | 3 |

CV is sample standard deviation divided by mean elapsed time. Each row below
summarizes all 156 operation series in that invocation, including Fortran:

| Invocation | Median series CV | Maximum series CV |
| --- | ---: | ---: |
| 2d-small-ordered | 1.85% | 3.84% |
| 2d-small-shuffled | 1.89% | 7.50% |
| 2d-large-shuffled | 2.85% | 37.19% |
| 3d-small-ordered | 0.85% | 3.28% |
| 3d-small-shuffled | 2.50% | 240.89% |
| 3d-small-shuffled-repeat | 0.95% | 5.27% |
| 3d-large-shuffled | 2.27% | 19.65% |

The original small shuffled 3D run contains a severe CBS12 outlier and noisier
CBS45/56 samples. The full repeat reduces variation and reproduces the
high-order gains. Every sample is preserved. This is an interactive desktop;
near-ties and percent-level differences are inconclusive. Reduced gains on the
large 3D field are consistent with a greater memory-access cost, but no hardware
counter measurement isolates that cost here.

All large-field invocations use shuffled markers; no ordered large-field
comparison was measured. The 3D marker placement covers four x-cell bands.
The 8.548 MiB CBS65 unique-entry count is not a cache footprint: under aligned
array assumptions, the union is 15.533 MiB with 64-byte lines or 28.359 MiB with
128-byte lines. Components execute sequentially, so neither aggregate alone
establishes capacity misses. See the follow-up for component footprints and
the limits of the reviewer's cache model.

Direct paths use the same component buffers, marker data, selected indices and
explicit zero shifts. C++ caches strides and inverse cell volume in its geometry
object; Fortran retains its existing timed prologue. Patch C++ constructs
geometry and uses interleaved values and empty shifts; LEInteractor also selects
indices and copies components. Clipping is correctness-tested; the timed cases
use complete stencils. These results do not establish whole-solver speedups.

## Reproduce and inspect

Use the native configuration and build commands in [README.md](README.md).
The optional final `CBS` argument selects all ten composites; omitting it runs
the complete IB/BS/CBS benchmark collection. For example:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  build/ibkernels-matrix-free/Release-native/tests/benchmark-matrix-free-3d \
  128 65536 3 9 1 CBS
```

Select a mode explicitly at the component entry point when needed:

```cpp
coupling.spreadAxis<0, double, IBTK::Experimental::TensorProductMode::FACTORIZED>(
    kernel, field, positions, indices, shifts, values);
```

Evidence root in this checkout:
`evidence/experiments/ibkernels-matrix-free-20260915/factorization/`.

- `source.json`, `source.patch`, `setup.json`: exact revision, binary hashes and baseline copies.
- `pilot/` and `matrix/`: raw CSVs, commands, run statuses and machine-load snapshots.
- `summary-by-run.csv`: all medians, ranges, CVs and matching Fortran/expanded ratios.
  Separate summaries preserve the repeat as its own invocation.
- `attest.json`, dimensional outputs, build logs and `repeat-checksums.json`: validation.
- `native-cpp-{2,3}d.s`, compiler commands/remarks and `assembly-observations.json`: generated code.
- `artifact-verification.json`, `final-custody.json`: dependency and final artifact checks.

Run `summarize.py` on each raw invocation separately to reproduce its summary;
combining repeated configurations in one invocation pools their samples.
All source and report changes are local to the experimental branch.
