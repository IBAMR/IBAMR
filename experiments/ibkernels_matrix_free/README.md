# Matrix-free IBKernels experiment

This opt-in experiment uses the owning-output evaluator API from qualified
R01d `c37b111fdabf267e83229c069698eca1a152bf1e`. It is independent of the CAV
implementation stack. Fortran kernels and production defaults are unchanged.

See [RESULTS.md](RESULTS.md) for the native Apple-toolchain comparison.

## Source

- `ibtk/src/lagrangian/experimental/SideCoupling.h` and its inline header provide
  a patch-geometry object with generic component gather/scatter operations.
  The header is source-private and is not installed.
- `tests/matrix_free/side_coupling.cpp` is the shared native 2D/3D regression.
- `tests/matrix_free/benchmark.cpp` is a standalone benchmark with no timing
  fixtures or CI performance threshold.
- `tests/matrix_free/coupling.h` demonstrates all-component use and a compiled
  application-defined cosine kernel. `fixture.cpp` contains independent scalar
  reference formulas, real SAMRAI patch construction and direct Fortran symbols.

The field has depth one and all side directions. Marker coordinates and vector
values are interleaved. Index lists are explicit; shifts are indexed by list
position. Gather overwrites selected marker components. Scatter adds values
divided by cell volume; the marker input must already contain any desired
Lagrangian quadrature factors. Stencils clip at allocated ghost bounds, without
renormalization. The caller owns selection, ghost synchronization, physical
boundary treatment and concurrent-spread synchronization.

The numerical operation uses exact stencil widths and requests
`IBKernels::Weights<Coefficient, stencil_size>` by value. Default coefficients
and field/marker arithmetic are double. The regression also requests float
coefficients while retaining double coordinates and field accumulation.
There are no marker-coordinate copies, matrices, per-weight callbacks or heap
allocations in the component loops. The nonoverlap contract is documented at
the entry point; there are no `restrict` promises inferred from concepts.

## Correctness

The compact native cases cover IB4, IB5, B-spline 3 and 6, both 3/2 and 2/3
normal/tangential composition, and an application-defined cosine functional
form. They exercise nonzero patch indices, anisotropic spacing, side centering,
face/center ties, ghosts, overlapping markers, additive output, selected and
repeated indices, shifts, clipping, and unchanged marker inputs. Dense reference
calculations use physical grid coordinates; B-splines use a truncated-power
formula independently of the evaluator recurrence. They also check the
cell-volume-scaled gather/scatter adjoint identity.

IB5's scalar Fortran delta supplies its reference values. The unmodified 3D
IB5 *spreading* routine omits the x-index update in its innermost loop and
accumulates contributions on the center x plane. Its nonzero discrepancy is
reported explicitly by the regression. Its interpolation and scalar delta are
valid comparisons. IB5 spreading is excluded from performance comparisons.

## Reproduce on this machine

The configuration script records installed dependency paths and uses the
existing shared ccache. It neither installs dependencies nor fetches anything.
Release loads `/Users/boyceg/code/autoibamr/opt/configuration/enable.sh` and uses
its rebuilt dependencies, including SAMRAI and Boost. An explicit
`EXPERIMENT_SAMRAI_ROOT` can select a newer SAMRAI installation.
C and C++ use the active Apple developer toolchain and SDK selected by `xcrun`;
Fortran uses GNU Fortran. Both configurations use applicable strict macOS
warnings. Debug uses `-O1` C++/`-O2` Fortran plus debug flags; Release uses
`-O3 -mcpu=native` for C, C++ and Fortran, with `NDEBUG` for C/C++.
Fast-math is explicitly disabled. Only Release runs support performance claims.
The native Release build has a separate persistent directory so earlier results
remain distinguishable from results with the rebuilt dependencies.

```sh
bash experiments/ibkernels_matrix_free/configure.sh Debug
export CCACHE_DIR=/Users/boyceg/Library/Caches/ccache
export CCACHE_BASEDIR="$PWD"
export CCACHE_TEMPDIR="$PWD/.cache/ibkernels-ccache-tmp"
cmake --build build/ibkernels-matrix-free/Debug --target indent
cmake --build build/ibkernels-matrix-free/Debug --target tests-matrix_free -j4
cd build/ibkernels-matrix-free/Debug
../../../attest -N -R '^matrix_free/'
../../../attest -R '^matrix_free/'
```

From the source root, build optimized targets:

```sh
bash experiments/ibkernels_matrix_free/configure.sh Release
cmake --build build/ibkernels-matrix-free/Release-native \
  --target tests-matrix_free benchmark-matrix-free-2d benchmark-matrix-free-3d -j4
```

Run the benchmark from a durable working directory after builds finish and
machine load is suitable. Its arguments are cells per side, marker count,
iterations per sample, repeated samples, and shuffled marker order (0 or 1):

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  build/ibkernels-matrix-free/Release-native/tests/benchmark-matrix-free-3d \
  16 4096 20 9 0
```

Each CSV row records raw elapsed seconds and a consumed checksum. Three warmups
precede measurements. Implementation order rotates across repeats; output
reset and checksum evaluation occur outside timing. Spreading accumulates for
the same number of iterations in both implementations. Full output comparisons
precede timing, including direct Fortran ABI calls.

`cpp_loop` and `fortran_loop` use preselected indices, precomputed patch geometry,
identical scalar component buffers, double precision and the same positions and
zero shifts. The C++ geometry object also caches strides and inverse cell volume;
its construction cost is included in `cpp_patch`. The Fortran routine retains
its existing fixed prologue inside the timed call. Kernels and stencil locations
are evaluated during every call; there is no coefficient precomputation.
`cpp_patch` includes construction of
the geometry object and operates on interleaved values with caller-selected
indices. `LEInteractor` additionally includes its own index selection and
component copies. Those wrapper timings therefore describe their respective
entry-point costs, not identical wrapper internals. Each operation is one CPU
thread; there is no concurrent scatter.
The C++ patch path also uses an empty shift span instead of the direct paths'
explicit zero-shift buffer. Its layout and shift handling differ from
`cpp_loop`, so subtracting those timings does not isolate geometry setup cost.

Machine, compiler versions, effective compile commands, exact revisions,
native outputs, raw CSVs and result assessment are retained under
`evidence/experiments/ibkernels-matrix-free-20260915/` in the experiment checkout.
