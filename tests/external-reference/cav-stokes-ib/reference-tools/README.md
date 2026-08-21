# CAV N0 paired replay

This directory contains external-reference tooling for the mandatory N0
multiplicative-CAV comparison. It is not part of the normal `attest` suite and
does not replace the native CAV regressions.

`run_n0_paired_replay.sh` requires:

- a clean IBAMR worktree at the exact candidate commit;
- a build configured from that same worktree;
- a separate clean sandbox worktree detached at
  `5b77344db6746269f8c77695c99e9043907ba74b`; and
- MATLAB with the pinned sandbox dependencies available.

The runner regenerates both raw bundles. The candidate uses the live IBAMR
objects from the supported periodic `4x4 -> 8x8` hierarchy and observes the
finest `8x8` pressure-CAV level. The oracle exporter calls the pinned sandbox
construction and multiplicative smoother without editing the oracle. The
replay maps both bundles into the oracle's global velocity-pressure ordering
and applies the declared IBAMR pressure-equation-row sign map. It then executes
both exported patch sequences with one MATLAB arithmetic path on both the live
candidate and independently constructed oracle operator/RHS inputs.

The generated matrices, vectors, reports, logs, and hashes remain in the
runner's temporary output directory and must not be committed.

Example:

```sh
tests/external-reference/cav-stokes-ib/reference-tools/run_n0_paired_replay.sh \
  /tmp/ibamr-cav-durable-6f2c \
  /tmp/ibamr-cav-durable-6f2c-build-dbg \
  CANDIDATE_COMMIT_SHA \
  /tmp/implicit-ib-sandbox-cav-oracle-5b77344 \
  /Applications/MATLAB_R2025b.app/bin/matlab \
  /private/tmp
```

The final argument is an output parent. The runner creates a new directory
inside it and never deletes or overwrites an existing evidence bundle.
