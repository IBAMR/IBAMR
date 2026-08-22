# CAV multiplicative paired replay

This directory contains external-reference tooling for the mandatory N0, N1,
and N2 multiplicative-CAV comparisons. It is not part of the normal `attest`
suite and does not replace the native CAV regressions.

`run_multiplicative_paired_replay.sh` requires:

- a clean IBAMR worktree at the exact candidate commit;
- a build configured from that same worktree;
- a separate clean sandbox worktree detached at
  `5b77344db6746269f8c77695c99e9043907ba74b`; and
- MATLAB with the pinned sandbox dependencies available.

The runner regenerates every raw bundle. The candidate uses live IBAMR objects
from the supported periodic `4x4 -> 8x8` hierarchy and observes the finest
`8x8` pressure-CAV level. N1 treats that live level as a transfer-free level;
it does not add a special FAC path for a one-level hierarchy. The runner covers
`K=1`, `1e2`, and `1e4`, RELAXED and STRICT construction, and forward and
reverse common-arithmetic sweeps. Reverse order is a sensitivity replay, not a
production option.

N2 uses a supported periodic `8x8 -> 16x16` hierarchy and the live RELAXED
pressure-CAV construction on its finest level. Its production V-cycle has a
forward pre-sweep and a forward post-sweep. A forward/reverse palindrome is a
sensitivity replay only. N2 records the explicit `RT0_REFINE`/`RT0_COARSEN`
velocity transfers, `LINEAR_REFINE`/`CONSERVATIVE_COARSEN` pressure transfers,
the SVD pressure-nullspace coarse solve, Richardson level smoothing, both FAC
levels, transfer-stage actions, V-cycle state, FGMRES state, and original
physical residual components.

The oracle exporter calls the pinned sandbox RELAXED construction and
multiplicative smoother without editing the oracle. Since STRICT is an IBAMR
extension rather than a sandbox algorithm, the consumer independently derives
its patches from the canonical `E_h` row-or-column graph and the complete
incident-velocity rule. It likewise derives RELAXED patches as a construction
control. Candidate patches always come from the live production builder; the
consumer never substitutes its reference construction for a candidate export.

The replay maps all bundles into the oracle's global velocity-pressure
ordering, applies the declared IBAMR pressure-equation-row sign map, and
executes candidate and reference sequences with one MATLAB arithmetic path on
both the live candidate and independently constructed oracle operator/RHS
inputs. It compares local matrices, local right-hand sides, local corrections,
global states, and residuals after every patch. Separately labeled candidate
and oracle native-smoother comparisons remain cross-runtime diagnostics, not
the common-arithmetic gate.

For N2, the actual IBAMR solve is accepted only when it produces complete
finite physical/action diagnostics and a freshly recomputed original relative
residual no larger than `1e-10`. Candidate/oracle physical scalar differences,
iteration-count/history differences, and the matrix-free-versus-assembled
`rho_IB` value remain mandatory first-mismatch diagnostics. The fixed per-`K`
table applies to the common-arithmetic replay, not to those separately
accumulated cross-representation diagnostics.

The generated matrices, vectors, reports, logs, and hashes remain in the
runner's temporary output directory and must not be committed. The runner
preserves and hashes a generated report even when a comparison fails, so a red
bundle remains reproducible diagnostic evidence.

Example:

```sh
tests/external-reference/cav-stokes-ib/reference-tools/run_multiplicative_paired_replay.sh \
  /tmp/ibamr-cav-durable-6f2c \
  /tmp/ibamr-cav-durable-6f2c-build-dbg \
  CANDIDATE_COMMIT_SHA \
  /tmp/implicit-ib-sandbox-cav-oracle-5b77344 \
  /Applications/MATLAB_R2025b.app/bin/matlab \
  /private/tmp \
  N1
```

The optional output-parent argument defaults to `/private/tmp`. The final scope
is `N0`, `N1`, or `N2` and defaults to `N1`; N0 runs `K=1`, both construction
policies, and forward traversal, N1 runs the complete single-level matrix
above, and N2 adds the live two-level transfer/V-cycle/FGMRES path. The runner
creates a new directory and never deletes or overwrites an existing evidence
bundle.
