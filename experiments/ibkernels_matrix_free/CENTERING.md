# Cartesian centering in the matrix-free experiment

2026-09-26. The experiment uses one `CartesianCoupling<C>` implementation for
cell, node, side, face, and edge fields. The centering types and geometry come
from `CartesianCentering<C>` in PR #1985. `SideCoupling` is an alias for the
side specialization, preserving the existing experimental callers.

## Dependency and history

The four centering commits after `4de063aa3` through
`88d1773b2053608faf86cb13cb8da7290c220ead` were incorporated on top of the
validated #1997-based experiment at `ba75c65f0`. Their local counterparts are
`3accac972`, `287171234`, `0beb67a83`, and `6a1bc0a7f`; `git range-diff`
reports identical patches. Original authorship is retained.

Recovery ref `codex/ibkernels-before-centering-20260926` preserves `ba75c65f0`.
The existing historical Release builds, measurements, and recovery refs remain
in place. This change is local; no upstream branch or PR was updated.

## Component interface

The constructor copies patch geometry and allocation bounds. The caller selects
one depth plane with the native SAMRAI pointer accessor. `Axis` selects the data
direction: normal for side/face, tangent for edge, and zero for cell/node.
The optional fourth template argument, `KernelAxis`, selects the evaluator's
distinguished Cartesian direction independently of depth and centering.

For a face field with at least two coordinate directions:

```cpp
using namespace IBTK;
Experimental::CartesianCoupling<DataCentering::FACE> coupling(patch, field);
coupling.interpolateAxis<1>(kernel,
                           field.getPointer(1, depth),
                           positions,
                           indices,
                           shifts,
                           values);
```

Here `values` is a scalar marker buffer. For interleaved marker values, pass
the selected component pointer and the number of components as `marker_stride`.
For cell data with the kernel oriented along coordinate one, select
`interpolateAxis<0, double, Experimental::TensorProductMode::CONTRACTED, 1>`
and pass `field.getPointer(depth)`.

Cell and node data have one geometric array; staggered data have one per
allocated direction. Construction skips missing side directions, and an
operation requesting one fails before entering the marker loop. Data depth
does not determine kernel orientation or the CBS family.

SAMRAI face storage places the normal coordinate first, followed cyclically by
the remaining coordinates. Compile-time traversal and factor selection follow
that storage order to retain contiguous field access. Evaluator coordinates and
expanded coefficient indices remain Cartesian. Side, cell, node, and edge
storage use Cartesian traversal. No per-entry dispatch or allocation was added.

Spreading remains additive and scaled by inverse cell volume; interpolation
overwrites selected marker values. Clipping does not renormalize coefficients.
The existing marker selection, list-position shifts, and nonoverlap contracts
remain in force. Contraction order is deterministic; different modes and
storage orderings may differ by floating-point rounding.

## Verification

The existing persistent `build/ibkernels-matrix-free-pr1997/Debug` build used
Apple Clang 17, GNU Fortran 15.2, the installed Debug dependencies, shared
ccache, and strict macOS warning flags. The repository's required formatter ran
before compilation. Both benchmark executables compiled but were not executed.

- Native discovery found 19 cases. The focused run and repeat each passed 19/19.
- The original side/IB5 regressions retained their expected outputs.
- The two new fixtures reuse the existing 2D/3D executable. They cover every
  requested BS2–6 and CBS(k+1)k/CBSk(k+1) pair, all five centerings, all axes,
  unequal box extents and ghost widths, shifted indices, depth isolation,
  scalar marker stride, indexed shifts, clipping, and empty stencils.
- Representative asymmetric CBS cases exercise all three application modes
  and kernel orientation independent of the geometric axis. Face fallback and
  partial side allocations have additional cases.
- Dense scalar references use SAMRAI iterators and indexed access rather than
  the implementation's raw strides, stencil placement, or centering helpers.
  Maximum new-case gather, spread, and adjoint errors were approximately
  `1.11e-15`, `1.11e-16`, and `7.77e-16`, respectively.
- New expected outputs were generated from actual runs and inspected against
  those references. Each contains 22 compact numerical records. No tolerance
  was changed. The inherited pointwise-function fixture changes belong to the
  unchanged #1985 patches; all six corresponding cases passed here.
- Independent source and regression reviews found no correctness defects.
  One marker-component documentation example was clarified.

The focused selection also includes kernel evaluators and interpolation-matrix
construction, including two-rank cases. Logs, effective compiler commands,
source hashes, generated outputs, and dependency comparisons are under
`evidence/experiments/ibkernels-matrix-free-20260915/centering-20260926/`.
Stale generated links for fixtures removed or renamed by #1985 were removed
after inspecting their targets. The known duplicate-library/rpath and macOS
deployment-version linker warnings remain unchanged.

Release qualification, performance measurements, and profiling remain deferred.
The previous timing reports describe their original revisions and do not
establish performance of this extension. The generalized coupling remains a
source-private experiment; production spreading/interpolation defaults are unchanged.
