# Common Independent Audit Protocol and Launch Index

**Protocol version:** `cav-audit-v5` (2026-08-20)

**Candidate policy:** permanently read-only

**Result ownership:** independent audit coordinator; never the implementation task

This is the common committed protocol for every CAV audit lane. Launch each assigned lane as a separate Codex task with this file and exactly one focused lane prompt:

| Lane | Prompt |
|---|---|
| Algorithm semantics | `plans/implicit-ib/cav-audit-algorithm-prompt.md` |
| Numerical reproduction/parity | `plans/implicit-ib/cav-audit-numerical-prompt.md` |
| Data-structure efficiency | `plans/implicit-ib/cav-audit-data-structures-prompt.md` |
| Cache/data-layout efficiency | `plans/implicit-ib/cav-audit-cache-layout-prompt.md` |
| Concrete extensibility | `plans/implicit-ib/cav-audit-extensibility-prompt.md` |
| Simplicity/minimality | `plans/implicit-ib/cav-audit-simplicity-prompt.md` |

After all assigned reports are frozen, launch `plans/implicit-ib/cav-audit-synthesis-prompt.md` as a separate adjudication task. The implementation task must never launch itself in any auditor, coordinator, or synthesizer role.

Uncommitted drafts are not launch inputs. The exact candidate commit must contain the authoritative plan, canonical ledger, this protocol, the assigned versioned lane prompt, and the synthesis prompt before the coordinator may launch a suite.

## Permanent read-only role boundary

An auditor is permanently read-only with respect to the candidate checkout, candidate refs, pinned source refs, sandbox oracle source, comparator source, committed prompts, and another lane's evidence. This is not a default that can be relaxed inside the audit task.

An auditor may create clean build directories, run commands, transform copies of raw exports in a disposable audit harness, and write its own report and generated evidence only in the coordinator-designated audit-results workspace. It must not:

- edit, stage, commit, or otherwise change the candidate;
- implement or propose a patch as audit output;
- modify the oracle, comparator, formal authorities, or pinned refs;
- rewrite refs, push, open a pull request, or resolve a finding;
- continue as the implementation task after finding a defect.

Every fix belongs to the separate implementation task and produces a new candidate SHA. A newly launched independent auditor then reruns every affected lane. The same auditor task must not switch roles to fix or re-audit its own findings.

At conclusion, resolve and record the exact historical branch-7/8 refs and the oracle detached HEAD again, verify that each matches its launch value, and record final dirty state. Any unexpected ref/source change is a blocking `S1` reproducibility defect. Generated build/run artifacts must remain outside the candidate checkout and be either included immutably in the coordinator manifest or explicitly identified as disposable and excluded; no auditor decides to retain them in candidate history.

## Audit coordinator and immutable custody

A separate audit coordinator, never the implementation task and never a lane auditor for the same suite, owns launch completeness and immutable report custody. The coordinator does not alter findings, decide technical merit, or edit implementation code. It must:

1. assign a suite ID and exact candidate SHA;
2. record the Git commit and Git blob ID plus SHA-256 for this protocol, the authoritative plan, canonical ledger, synthesis prompt, and every assigned lane prompt;
3. launch lanes in separate clean worktrees at the same exact candidate SHA;
4. withhold every lane report from the other lane auditors until all assigned reports are submitted;
5. freeze reports and raw evidence under `audit-results/implicit-ib/cav/<candidate-sha>/<suite-id>/` in a separate append-only audit-results checkout, never in the candidate checkout;
6. create `MANIFEST.sha256` over every prompt, report, raw export, comparator result, log, and environment manifest retained for the suite;
7. reject overwritten files: corrections create a new suite ID or a clearly versioned addendum whose old hashes remain present;
8. verify required lane coverage, then launch the separate synthesis task after reports are frozen.

The suite manifest records report filenames, submitting task IDs, UTC submission times, prompt versions, candidate and parent SHAs, oracle and manuscript pins, build identities, and artifact retention policy. A report is not frozen until its SHA-256 appears in the manifest. Candidate code and generated audit artifacts are never mixed in one commit.

## Common launch block

Audit lane **[LANE]** for milestone **[MILESTONE]** at exact candidate SHA **[CANDIDATE_SHA]**, with declared parent **[PARENT_SHA]**, suite ID **[SUITE_ID]**, and required endpoint **[ENDPOINT]**. This is an independent audit task and the candidate is permanently read-only.

Create a fresh, clean IBAMR worktree pinned directly to `[CANDIDATE_SHA]`. Record repository and worktree paths, requested/resolved SHA, parent and foundation SHAs, all relevant stacked SHAs, branch/upstream identities, submodules, and full porcelain dirty state before and after work. Do not audit a moving branch tip.

Create a second, separate clean sandbox-oracle worktree detached exactly at `5b77344db6746269f8c77695c99e9043907ba74b`. It must not be the existing checkout at `/Users/boyceg/code/implicit-ib-sandbox-cav-smoothing`, whose current HEAD and dirty state are not the oracle state. A local no-network clone followed by detached checkout is acceptable. Record the source repository path, resolved detached HEAD, full dirty state, and submodules. If an exact clean oracle worktree cannot be established, return **INCOMPLETE**.

Read, in order:

1. every applicable `AGENTS.md`;
2. committed `plans/implicit-ib/cav-implementation-plan.md`, especially Sections 2–12;
3. committed `plans/implicit-ib/cav-historical-review-findings-ledger.md`;
4. the recovery manifest in plan Section 11;
5. this common protocol;
6. the assigned lane prompt only.

Do not require or use the original attachment path. Do not receive implementation conversational reasoning or another lane's report.

## Pinned authorities and their distinct roles

- **IBAMR architectural specification:** cleaned foundation `0c8447e3593113a8b62f94df9482c60d896f0d3b`.
- **Formal prototype algorithms at the oracle pin:** sandbox `docs/CAV_ALGORITHMS.tex` (Git blob `0900b817bf5abd40116d8164b641e4c732619375`) and `docs/CAV_RAS_PARALLEL_NOTES.md` (Git blob `fb2cc3417bff136b09d4c8cf6bf9e522a5e0b8bc`), both read from `5b77344db6746269f8c77695c99e9043907ba74b`.
- **Pinned manuscript context:** `paper-coupling-aware-vanka` commit `9c523e0135b880758e773b373b9e6b787b7bb486`, `manuscript-body.tex` (Git blob `4955185c61cf02a0f1d3776b1f50b92d7b3d2b07`), especially “Coupling-aware Vanka patches” and “Multiplicative patch relaxation.” It is mathematical context for pressure-centered construction and multiplicative relaxation; it does not define CAV-RAS production details absent from that manuscript.
- **Executable algorithmic evidence:** the detached sandbox worktree at `5b77344d...`, including its validation tests. Executable agreement is evidence, not the definition of mathematical correctness.
- **Historical evidence only:** branch-7 commit `4f83a028013e6ce8d716ea74e25b4c194c9d587f` and branch-8 commit `d400a07ebe9809c581adbf0d3a3778239ae5d434`.
- **Recovery disposition authority:** the committed canonical ledger.

Verify every ref and record local/remote ambiguity. No moving branch name can substitute for an exact commit.

### Mandatory seed-formulation reconciliation

The formal algorithm and manuscript start from each pressure cell's standard Vanka velocity set and use the row-or-column graph of the Eulerian elasticity block. Cleaned-foundation IBAMR instead enumerates one selected velocity component in a declared logical order and constructs RELAXED or STRICT closures from rows of the full level `A00`, which contains fluid/Stokes terms as well as IB coupling. The implementation plan therefore requires a distinct pressure-seeded `E_h` formal path and retains the foundation path only as explicitly named legacy compatibility behavior. Their equivalence must not be assumed or manufactured by relabeling.

Before patch or smoother parity can pass, the assigned A/N audit evidence must state and verify:

1. that the formal path receives the live Eulerian elasticity block `E_h` separately and never uses fluid/Stokes nonzeros to enlarge patches;
2. that row-or-column adjacency is implemented without an unproved symmetry shortcut;
3. the exact pressure-cell seed, standard Vanka velocity set, logical order, periodic/boundary handling, seed stride, and duplicate policy;
4. that RELAXED exactly implements the pinned targeted-IB construction, including an exact standard patch when no `E_h` edge leaves `U_i`;
5. that STRICT is reported as the declared IBAMR extension and that any legacy full-`A00` velocity-seeded option remains separately named and tested;
6. that full original saddle-point matrices—not `E_h`—remain the source for local solves and original-residual updates.

A mismatch must be reported; it cannot be hidden by a final smoother norm.

## Evidence provenance and comparator independence

Rerun evidence from builds tied to the active audit checkout. Old logs and implementation claims are context, not current validation. Use live IBAMR objects for production operators, mappings, patches, transfers, smoothers, V-cycles, and residuals. Generated evidence records commit, build, dimension, grid, coefficients, MPI size, mapping version, gauge, stiffness, traversal, backend, precision, and exact regeneration command.

Numerical acceptance has two mandatory, separately reported tracks. The common-arithmetic track replays candidate and oracle semantics in one documented runtime on both independently exported operator/RHS sources, while consuming each side's own exported ordered construction/application artifacts. The live-production track executes the actual IBAMR path and reports fresh original momentum, divergence, and IB consistency residuals plus convergence and iteration history. Common replay cannot replace live construction or execution; raw cross-runtime entrywise differences remain first-mismatch diagnostics but are not, by themselves, the algorithm-parity definition. The exact formulas, fixed low-`K` targets, and pass/fail rules are in the versioned N-lane prompt.

The numerical lane must independently validate the exporter, mapping, and comparator before accepting parity:

- retain raw, pre-mapping exports from both live IBAMR and the detached oracle;
- hand-check a tiny live `4x4` 2D case, including selected velocity/pressure DOFs, one interior or periodic patch, local matrix rows, and one patch correction;
- run comparator unit tests on exact-match and independently produced raw files;
- mutate copies—not candidate/oracle sources—to exercise a DOF permutation, velocity or pressure sign flip, legal pressure-only constant gauge shift, illegal velocity shift, omitted row/DOF/patch, and patch reordering;
- require declared mappings and legal gauge normalization to pass, while undeclared permutations, sign errors, omissions, illegal shifts, and algorithmically meaningful reorderings are detected;
- record raw hashes, mapped hashes, comparator version/hash, mutation activation proof, expected result, and actual result.

A comparator or mapping that has not passed these controls makes the numerical lane **INCOMPLETE**.

## Native CI and source-conformance evidence

Section 7 of `cav-implementation-plan.md` is the canonical same-unit native-test, runtime-budget, executable/fixture reuse, audit-hook, and IBAMR comment/style contract. Lane N owns actually configuring/building the active checkout and running both the focused feature cases and the appropriate broad candidate-checkout regression through the repository-root `attest` workflow. The canonical broad form is `/path/to/attest --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R IB/`; recorded commands use the checkout's actual runner path and preserve explicit tool paths required by the environment. For each scope, N records the selector, exact command, discovered cases, input/output fixtures, expected status/exit semantics, budget, measured wall time, deterministic rerun, unmutated result, meaningful failure sensitivity, and exclusions. A focused pass or built CMake target alone is not broad CI-integration evidence. External parity is separate evidence and cannot replace native `attest` runs. Assigned A/D/C/E lanes verify that comments accurately state their mathematical, ownership/invalidation, hot-path, and concrete-extension invariants. Lane S challenges unnecessary executables, duplicated harness/helper code, elaborate test APIs, and excessive or redundant comments. Synthesis verifies that feature and test were not split across PR boundaries and that the focused and broad `attest` sets are discoverable and CI-suitable.

A missing, unregistered, nondeterministic, nonsensitive, or unjustifiably slow required native test is blocking or makes the responsible lane **INCOMPLETE** under the shared taxonomy. A materially misleading comment is a defect; optional comment-density preference is not.

## Standard finding taxonomy

Classification and severity are separate fields.

| Classification | Meaning |
|---|---|
| `confirmed defect` | Reproduced incorrect behavior, violated contract, or demonstrated material regression. |
| `untested risk` | Required evidence is absent or inconclusive; correctness is not established. |
| `optional cleanup` | Non-required improvement with no demonstrated acceptance impact. |

| Severity | Gate effect |
|---|---|
| `S0 critical` | Wrong algorithm, corrupted result/state, or unsafe resource behavior; always blocking. |
| `S1 high` | Required correctness, reproducibility, supported-scope robustness, or material performance criterion fails; always blocking. |
| `S2 medium` | Material scoped weakness or test/evidence gap; blocking when it touches a milestone requirement, otherwise synthesis records its disposition. |
| `S3 low` | Local cleanup or clarity issue with no acceptance impact; nonblocking unless it exposes a broader defect. |

Lane outcomes are:

- **PASS:** all required evidence is complete and valid, and no blocking finding remains;
- **FAIL:** at least one confirmed blocking defect or mandatory criterion failure exists;
- **INCOMPLETE:** required evidence, environment, scope, or sensitivity control is missing or invalid, without enough evidence to assert failure.

`INCOMPLETE` is never a pass and blocks the applicable gate.

## Common finding and conclusion schema

Report findings in priority order. Each finding includes:

- stable finding ID and lane;
- classification and `S0`–`S3` severity;
- blocking status and affected milestone/review unit;
- exact candidate SHA and file/line or immutable artifact/hash;
- reproduction command and required environment;
- expected versus actual behavior;
- impact within this lane;
- recommended disposition, without implementing or closing it.

Conclude with the lane question answered; scope and exclusions; evidence/tests/measurements rerun; sensitivity exercises and controls; limitations and unresolved risks; refs/worktrees/artifacts state; highest parity layer reached; and exact-SHA **PASS**, **FAIL**, or **INCOMPLETE**.

## Production-robustness ownership within the six lanes

Production robustness is mandatory cross-lane work, not an unowned seventh concern. The audit coordinator verifies coverage; the named primary lane owns the conclusion and supporting lanes supply their scoped evidence.

| Concern | Primary owner | Supporting lane(s) |
|---|---|---|
| Same-unit fast native test registration, execution, runtime, determinism, and sensitivity | N | S |
| Debug and Release 2D/3D build and focused smoke tests | N | D |
| Advertised MPI ranks and rank-consistent results | N | A, D |
| initialization, deallocation, repeated reinitialization | D | N, C |
| hierarchy regridding, cache invalidation, stale-source rejection | D | N, C |
| local-solver singular/ill-conditioned/failure behavior | N | A, D |
| NaN/Inf propagation, sanitizer/resource/leak checks | N | D, C |
| empty/tiny/boundary/periodic/no-IB and ownership boundary cases | N | A, D |

An unsupported configuration must be explicitly rejected or excluded, not silently treated as passing. The multiplicative M5 gate supports full 2D numerical parity on MPI rank 1; MPI rank 2 and 3D are native robustness scopes unless the candidate explicitly advertises broader parity. The serial CAV-RAS M9 gate is rank 1. Distributed CAV-RAS requires its later independent gate.

## Executable performance contract

The C and D lanes use the following reproducible contract at production-code units and all-six gates:

- **Build:** `Release`, same compiler, flags, linked PETSc/BLAS/LAPACK, MPI implementation, and assertions for candidate and baseline; record all versions and the complete configuration command.
- **Hardware:** same host, CPU model, physical/logical cores, RAM, OS, power/frequency state, process affinity, and MPI binding; no mixed-host comparison.
- **Workload P1:** 2D periodic membrane, single `64x64` level, `K=1e4`, RELAXED, forward global multiplicative sweep, canonical parity backend, one rank; 5 warmups and 20 measured sweeps.
- **Workload P2:** 2D periodic membrane, three levels `16x16`, `32x32`, `64x64`, `K=1e4`, one pre- and one post-sweep plus one V-cycle, one rank; 3 warmups and 15 measured cycles.
- **Workload P3 at M9:** P1 geometry with CAV-RAS, pressure-cell macro blocks `8x8`, standard seed ghost width 0, IB seed ghost width 2, forward inner traversal, one rank; 5 warmups and 20 measured sweeps.
- **Setup workload:** 10 clean initialize/deallocate cycles on P1; report setup time, peak RSS, retained bytes, and allocation count where tooling supports it.
- **Baselines:** exact parent SHA for the review unit; cleaned `0c8447e...` for unchanged foundation operations; and the last frozen all-six tip for later cross-cutting changes. Missing algorithms get absolute measurements, not invented ratios.
- **Metrics:** median and MAD of wall time, time per sweep/cycle, time per active DOF, allocations and bytes per apply, setup time, peak RSS and bytes per DOF, plus scaling observations at `N=32,64,128` when a scaling claim is made.
- **Material regression:** blocking when median normalized time or peak memory is more than 10% worse and the delta exceeds three times the larger baseline/candidate MAD (or timer resolution), when new per-apply allocation/recomputation scales with patches/DOFs without necessity, or when an asserted linear-scaling path shows a reproducible superlinear trend. Report smaller changes without declaring a defect.

If build/hardware/workload parity cannot be established, the performance conclusion is **INCOMPLETE**, not speculative.

## Audit-suite sensitivity controls

The candidate and oracle remain untouched. Use precommitted fault fixtures or mutations of disposable harness state/raw copies. Every exercise records mutation ID/hash, activation proof, unmutated control, expected detector, and actual detector.

| Lane | Required sensitivity exercise |
|---|---|
| A | stale residual, patch reorder, gauge error, RAS ownership/overlap error, and inner additive substitution |
| N | the A mutations plus mapping permutation/sign/gauge/omission/reordering controls and first-mismatch localization |
| D | disposable forced per-apply map rebuild, duplicated full-domain ownership metadata, or equivalent lifecycle/recomputation challenge |
| C | disposable cache-disable/per-patch-allocation or locality-degrading traversal challenge measured against its unmutated control |
| E | structured change exercises for each concrete axis, including a seeded cross-layer-entanglement example the lane must reject |
| S | structured deletion challenge plus a seeded redundant abstraction/state/option example the lane must identify |

Failure to activate a mutation is **INCOMPLETE**. Failure to detect an activated required mutation is a blocking audit-suite-readiness failure.

## Milestone entry points and all-lane gates

Replace every placeholder with an exact SHA before launch.

| Entry point | Candidate/baseline scope | Required lanes |
|---|---|---|
| Foundation restack | `[FOUNDATION_TIP]` versus declared stack parents and canonical ledger | Plan Section 12.4 rows for affected foundation units |
| Recovered multiplicative review units | each exact lettered tip from `cav-1a` through `cav-6b` versus its exact parent | Assigned plan-matrix row only |
| **Working multiplicative production gate** | exact durable `cav-6b-fac-vcycle-fgmres-parity` tip after **M5**, also demonstrated equivalent to the disposable integration tip | **all six independently at the same SHA** |
| Review-stack equivalence | durable M5 tip versus `[INTEGRATION_TIP]` | N mandatory; every lane affected by differences |
| RAS construction | exact `cav-7-macro-construction` tip versus durable M5 tip | `cav-7` row |
| Macro-local multiplicative semantics | exact `cav-8a-macro-local-multiplicative` tip versus `cav-7` | `cav-8a` row |
| Outer restricted-additive semantics | exact `cav-8b-outer-restricted-additive` tip versus `cav-8a` | `cav-8b` row |
| Serial RAS production integration | exact `cav-9a-cav-ras-fac-production` tip versus `cav-8b` | `cav-9a` row |
| **Working serial CAV-RAS production gate** | exact `cav-9b-cav-ras-parity-acceptance` tip after **M9** | **all six independently at the same SHA** |
| Distributed RAS | `[DISTRIBUTED_RAS_TIP]` versus frozen serial M9 tip | A, N, D, C, E; S whenever API/state/abstraction changes |

The M5 all-lane gate must pass before M6 retires integration scaffolding or M7 begins. The M9 all-lane gate certifies serial CAV-RAS only; it does not certify distributed execution.

## Re-audit

Every fix creates a new candidate SHA and a new or explicitly versioned suite. Rerun every affected lane. A and N are mandatory after semantic or cross-cutting changes. D and C rerun after ownership, storage, lifecycle, traversal, resource, or hot-path changes. E and S rerun after API, option, abstraction, state, backend, or review-boundary changes. Freeze new reports and rerun synthesis. Prior lane or milestone outcomes never transfer automatically.
