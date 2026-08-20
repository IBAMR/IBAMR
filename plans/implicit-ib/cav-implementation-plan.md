# Multiplicative CAV and Macro-Subdomain CAV-RAS Implementation Plan

**Status:** authoritative and active

**Created:** 2026-08-18

**Foundation:** `0c8447e3593113a8b62f94df9482c60d896f0d3b`

**Plan branch:** `codex/implicit-ib-cav-0-recovery-plan`

> **This document is the authoritative execution plan for the CAV implementation. Update it whenever implementation findings materially change the planned branch structure, algorithm contract, audit strategy, or acceptance criteria. Do not silently depart from it.**

## 1. Goals and non-goals

### Goals

1. Recover a correct current-master implementation of global multiplicative coupling-aware Vanka (CAV) from the cleaned IBAMR foundation.
2. Prove that implementation against the sandbox algorithmic oracle, layer by layer, using matrices, vectors, patch definitions, and operators exported from live IBAMR objects.
3. Reconstruct the complete capability path from historical branches 7 and 8 as small, independently testable, mergeable review units.
4. Add macro-subdomain CAV-RAS only after multiplicative CAV parity is green: multiplicative CAV within each macro-subdomain and restricted-additive composition across macro-subdomains.
5. Preserve reviewability and numerical correctness simultaneously: every semantic claim receives a focused native test and a matching parity check at the earliest useful layer.
6. Report original-system physical residuals, separated into momentum, divergence, and immersed-boundary (IB) components, for end-to-end Krylov acceptance.

### Non-goals

1. Do not rewrite, move, force-update, or otherwise modify any historical branch or tag.
2. Do not merge the disposable integration reconstruction wholesale into the durable review stack.
3. Do not wholesale cherry-pick historical branches 7 or 8; recover capabilities by intent after re-evaluating their proper review boundaries.
4. Do not modify `/Users/boyceg/code/implicit-ib-sandbox-cav-smoothing` or treat its working tree as an IBAMR implementation staging area.
5. Do not begin CAV-RAS implementation before current-master global multiplicative CAV reproduces the sandbox and the durable review decomposition is established.
6. Do not substitute external surrogate matrices, synthetic transfers, or reimplemented operators for live IBAMR objects in exports or parity tests.
7. Do not promote the initial `K = 1e8` case to an acceptance gate; it is diagnostic until conditioning and attainable accuracy are understood.
8. Do not push branches, open pull requests, or commit generated parity/build artifacts without explicit user approval.

## 2. Sources of truth and distinct roles

All identifiers and repository states in this table were checked on 2026-08-18 rather than assumed.

| Source | Verified state | Distinct role |
|---|---|---|
| Cleaned IBAMR foundation | Commit `0c8447e3593113a8b62f94df9482c60d896f0d3b` exists and is the tip of local and remote `codex/implicit-ib-6.4-stokes-eigen-schur-shell-backend-restacked`; the worktree began detached exactly at this commit. | Only base for new integration and durable review branches. It supplies the cleaned 1.1–6.4 solver, transfer, CAV patch-construction, shell-backend, and FAC foundation. |
| Historical branch 7 | Commit `4f83a028013e6ce8d716ea74e25b4c194c9d587f` exists at `codex/implicit-ib-7-shell-smoother-blas-lapack-restacked`. Its merge base with the cleaned foundation is `bd010f326eebdb6b23dbde5c2779789c6586de9f`, not the cleaned foundation. | Intent and implementation evidence for BLAS/LAPACK shell local-solve modes and surrounding historical recovery work. It is not a directly mergeable continuation of 6.4. |
| Historical branch 8 | Commit `d400a07ebe9809c581adbf0d3a3778239ae5d434` exists at `codex/implicit-ib-8-live-parity-audit-tests-restacked` and descends from historical branch 7. | Intent, fixtures, and evidence for live operator audits, CAV parity, residual sweeps, and conditioning-aware thresholds. It is an audit source, not a wholesale cherry-pick source. |
| Sandbox executable oracle | Source repository `/Users/boyceg/code/implicit-ib-sandbox-cav-smoothing`; exact pin `5b77344db6746269f8c77695c99e9043907ba74b`. The existing checkout is at `115aea99531c1bd3448f46ef5a4479ac172faf53` with unrelated dirty files, so neither implementation nor audit may use it as the pinned working tree. Every parity run uses a separate clean worktree or local clone detached at the exact pin and records an empty dirty state. | Executable evidence for global multiplicative CAV and macro-subdomain CAV-RAS. Agreement is necessary evidence but is not by itself the definition of mathematical correctness or IBAMR architecture. |
| Formal CAV/CAV-RAS algorithms | At sandbox pin `5b77344d...`: `docs/CAV_ALGORITHMS.tex`, blob `0900b817bf5abd40116d8164b641e4c732619375`, and `docs/CAV_RAS_PARALLEL_NOTES.md`, blob `fb2cc3417bff136b09d4c8cf6bf9e522a5e0b8bc`. | Formal pressure-cell-seeded CAV and inner-multiplicative/outer-restricted-additive CAV-RAS definitions, plus the intended distributed ownership/data contract. These must be derived and checked independently, not merely executed. |
| Authoritative manuscript context | Clean local `paper-coupling-aware-vanka` checkout at commit `9c523e0135b880758e773b373b9e6b787b7bb486`; `manuscript-body.tex`, blob `4955185c61cf02a0f1d3776b1f50b92d7b3d2b07`, especially “Coupling-aware Vanka patches” and “Multiplicative patch relaxation.” | Mathematical context for pressure-centered patch construction and multiplicative relaxation. It does not define later CAV-RAS production details and cannot silently override the pinned formal notes. |
| Historical review findings ledger | Canonical repository copy `plans/implicit-ib/cav-historical-review-findings-ledger.md`, incorporated from the original attachment whose SHA-256 is `2eb0d1b842a503c7850946d88f17bebf2ee822f6feceb0ce9b7ac39b5258a6b6`. | Durable authoritative record of cleanup findings, dispositions, and design intent. Audits use the committed canonical copy and never depend on the ephemeral attachment path. |
| Live IBAMR objects in this worktree/build | Objects created through IBAMR's actual assembly, patch hierarchy, DOF maps, transfers, smoother, FAC, and Krylov paths. | Sole source for exported IBAMR matrices/vectors/patches and all production parity evidence. |

## 3. Historical review findings ledger

The complete canonical ledger is `plans/implicit-ib/cav-historical-review-findings-ledger.md`; the table below is its execution summary. The dispositions are preserved from the handoff ledger. A recovery change must not undo the associated design intent, and an auditor requires no external attachment to recover it.

| # | Finding and code area | Disposition | Implementation evidence | Preserved design intent / future action |
|---|---|---|---|---|
| 1 | Add tests in branch 1.1 support APIs: `ibtk/include/ibtk/PETScMatUtilities.h` and related `IndexUtilities`, `IBImplicitStrategy`/`IBMethod`, and `StaggeredStokesPETScLevelSolver` support | **Intentionally deferred for this cleanup pass; not accidentally missed.** | None | Add only if explicitly requested, and then in the review unit that introduces the support APIs. |
| 2 | Factor shared Stokes-IB operator plumbing in nonlinear and Jacobian paths | **Implemented.** | `3fded45fc`, `d4a736372` | Keep public APIs unchanged; share private/common time, velocity-state, Stokes-apply, and boundary plumbing so both paths retain identical semantics. |
| 3 | Collapse duplicate FAC residual vector handling in `computeResidual()` | **Implemented.** | `d65dd80fe` | Centralize PETSc/SAMRAI vector lifetime and copy/synchronization handling while preserving distinct rediscretized and non-rediscretized residual formulas. |
| 4 | Reuse transfer-entry construction instead of two full passes | **Implemented.** | `52dd7d9e4` | Share per-row candidate construction between PETSc preallocation and insertion; retain separate PETSc phases and do not add row caching without a separate reason. |
| 5 | Deduplicate velocity-block embedding in `StaggeredStokesPETScLevelSolver` | **Implemented.** | `efe5bc430` | Use one local row-copy helper; keep matrix assembly at call sites. |
| 6 | Hide mutable coupling-aware map internals | **Implemented.** | `2e77def20` | Keep `PatchLevelCellClosureMapData` public and `clear()` public; keep internals private and expose only narrow test access. |
| 7 | Share Eigen shell smoother sweep logic | **Implemented.** | `b222a691f`, `75c612f55` | Shared Eigen shell infrastructure owns traversal, correction application, residual update, and postprocessing; the Stokes Schur backend owns only block algebra and local solve behavior. |
| 8 | Consolidate pressure/side transfer test helpers | **Implemented.** | `70f2fae95` | Share parsing and common checks while keeping geometry-specific profile filling separate. |

## 4. Algorithm contracts

### 4.0 Matrix/live interaction-kernel consistency

The original matrix-free Fortran kernels are the source of truth for every IB interaction kernel. Specifically, the authoritative production paths are `ibtk/src/lagrangian/fortran/lagrangian_interaction2d.f.m4`, `ibtk/src/lagrangian/fortran/lagrangian_interaction3d.f.m4`, and the scalar definitions they include from `ibtk/src/lagrangian/fortran/lagrangian_delta.f.m4`. For **every** regularized-delta kernel supported by the PETSc matrix-based interpolation/spreading builder, the C++ matrix builder must follow the applicable original Fortran path's fractional-coordinate convention, support/stencil convention, weight formulas, branch choices, and arithmetic order, in both 2D and 3D. This includes component/transverse kernel pairs used by composite kernels. When an optimized interaction routine expands a recurrence instead of calling a scalar delta function, that optimized live recurrence—not an algebraically equivalent auxiliary formula—is authoritative for the matrix port.

The production invariant is matrix/live agreement for the same coordinates and grid, not agreement obtained by changing geometry, copying a matrix, adding an alternate accepted operator, or normalizing away a mismatch in the comparator. Existing matrix-free Fortran implementations must not be edited merely to reconcile the C++ matrix builder. A live-kernel change is permitted only after a focused test independently reproduces an actual live defect and distinguishes it from a matrix-port defect. The 3D IB5 spread correction is the current sole exception: a matrix-transpose action check independently exposed the stale `ic0` loop index, and the corrected five-point x support replaced a fixture that encoded that bug.

The exhaustive implementation inventory defines the shared mandatory set as `PIECEWISE_LINEAR`, `BSPLINE_3`, `BSPLINE_4`, `BSPLINE_5`, `BSPLINE_6`, `COMPOSITE_BSPLINE_32`, `COMPOSITE_BSPLINE_43`, `COMPOSITE_BSPLINE_54`, `COMPOSITE_BSPLINE_65`, `IB_3`, `IB_4`, `IB_5`, and `IB_6`. The composite suffix gives the component-axis/transverse pair. Live-only `COMPOSITE_BSPLINE_23/34/45/56`, `PIECEWISE_CONSTANT`, `PIECEWISE_CUBIC`, `DISCONTINUOUS_LINEAR`, `IB_4_W8`, and `USER_DEFINED` have no current PETSc matrix-builder counterpart and must be reported as unsupported by matrix-based Jacobian construction, not silently counted as passing this contract.

This consistency work is a separate foundational review claim, `cav-1c-matrix-live-kernel-consistency`, before patch construction. Its focused native tests must exercise interpolation and spreading through live IBAMR objects and the matrix builder on tiny hand-checkable periodic cases for every supported kernel and dimension; compare stencil indices, ordered weights, action values, transpose/spreading behavior, constants, and boundary/half-cell cases; and include controlled formula/order mutations that the tests detect. The fixed per-entry matrix/live bound is `1024 epsilon_machine` times `max(1, |a|, |b|)`, approximately `2.28e-13` in binary64; this accommodates measured high-order B-spline cancellation without masking a coefficient, support, orientation, or ownership defect. The test inventory must identify unsupported or live-only kernels explicitly instead of silently treating them as covered. Production code must not retain a duplicate interpolation/spreading matrix or introduce repeated index translations to satisfy these tests.

#### Fortran reuse decision

Do not construct matrix rows by applying the full Fortran interpolation/spreading kernels to basis fields. That would require scratch patch arrays, repeated zeroing/copying, and one or more hot-kernel calls per marker/row merely to recover weights already known during setup. The included Fortran scalar functions are callable for piecewise linear, B-spline 3--6, and IB3--IB5, but they are not a universal single implementation: optimized IB4/IB5 interaction paths expand weight recurrences directly, and IB6 has no corresponding scalar function. Redirecting only the C++ builder to those helpers would therefore follow an auxiliary formula rather than eliminate duplication with the actual matrix-free path.

The recovery implementation consequently keeps the C++ setup-time weight callbacks as small, explicitly Fortran-authoritative ports guarded by the exhaustive native cross-path test. Introducing a shared Fortran weight-vector ABI and refactoring every optimized live kernel to call it is not part of CAV recovery: it changes the matrix-free hot path, requires separate 2D/3D ABI, inlining, performance, and lifecycle evidence, and conflicts with the rule against touching correct live kernels. It may be proposed later as its own behavior-preserving review unit only if measured maintenance or performance evidence outweighs that complexity. No such refactor is needed for multiplicative CAV or CAV-RAS parity.

Cross-language geometry generation is outside this claim. Coordinates that differ by a trigonometric rounding unit are distinct inputs; kernel consistency is assessed first with bitwise-identical coordinates. Oracle parity then separately measures the effect of the declared IBAMR-to-oracle geometry mapping without rewriting either geometry generator.

### 4.1 Base-patch construction: RELAXED and STRICT

Patch construction is deterministic data production, separate from patch application. Its output includes ordered seed traversal plus overlapping patch DOFs and, where required by composition, owned/restricted DOFs.

- The pinned formal/oracle construction is **pressure-cell seeded**. For pressure seed `i`, let `U_i` be its standard Vanka incident-velocity set and let `E_h = S E_L J` be the live Eulerian elasticity block. The construction graph is the numerical row-or-column nonzero graph of `E_h` only; the fluid mass/viscous/Stokes graph is not a patch-enlargement graph.
- In IBAMR, construct and retain patches in the full global velocity-pressure DOF index space. The live `E_h` must use that same full space; its pressure rows and columns must be structurally or numerically zero and are not part of the graph. Scan the live matrix once during setup into an undirected row-or-column velocity adjacency indexed directly by global DOF number. Do not extract or copy a velocity block and do not materialize a transpose. Compressed patch-local numbering begins only inside the local solution backend after global patch membership/order is fixed.
- Treat that full global numbering as the canonical persistent algebraic representation. Matrices, global vectors, seed identifiers, patch memberships, residuals, and corrections remain in it across construction, composition, FAC orchestration, exports, and parity mapping. Do not maintain parallel global and velocity-subdomain matrix representations. Translation to compact patch-local numbering is permitted only at the local solve boundary, where the backend needs a local factorization/solve, and the resulting correction is scattered directly back by the same patch index list.
- Reuse live PETSc matrices by borrowed/non-owning reference whenever the owner lifetime covers the operation. Do not duplicate `A`, `E_h`, a velocity block, a transpose, or a backend-format matrix merely for traversal, inspection, export, or parity. An unavoidable matrix conversion or copy must have a named consumer that cannot use the live representation, a documented lifetime/ownership boundary, a focused test, and a measured setup/apply-time and memory cost. Exporters must stream the live matrix and never become a second production representation.
- At the local-solve boundary, retain only the data required to apply the selected factorization. A dense local factor buffer and pivots are permitted when the factorization cannot consume the live sparse submatrix in place; retaining both an unfactored local matrix and an equivalent factor copy, copying the full level matrix, or caching duplicate residual-update blocks requires separate measured justification. Apply-time workspaces must be allocated during setup and reused rather than reconstructed or copied per patch application.
- `RELAXED` is the exact formal targeted-IB rule. Expand `U_i` through row-or-column `E_h` adjacency. If no velocity outside `U_i` is found, the patch is exactly `U_i union {p_i}`. Otherwise take every pressure cell incident to the expanded velocity set and the complete pressure-plus-incident-velocity closure of those cells, retaining every expanded velocity DOF. This is the multiplicative sandbox-parity target.
- `STRICT` is an explicit IBAMR extension on the same pressure seed and `E_h` graph. After the expansion, retain only candidate pressure cells whose complete incident-velocity closure is contained in the expanded velocity set; the seed cell remains present because `U_i` is the initial set. It may produce smaller patches and must keep its strict/classical native compatibility checks. The pinned sandbox does not independently define this option, so it must never be represented as executable-oracle evidence that the sandbox itself supplies.
- The cleaned foundation's existing **velocity-seeded** RELAXED/STRICT construction is retained as historical compatibility evidence during integration, not relabeled as the formal algorithm. It selects one velocity component, expands rows of the full level `A00`, then closes through adjacent cells. Static inspection confirms that production `A00` contains fluid/Stokes terms in addition to the IB augmentation; therefore this family is not equivalent to the pinned `E_h` pressure-seeded family by mapping alone.
- Pressure-seed stride and logical traversal order are explicit algorithm inputs. The sandbox's 2D pressure order is mapped to the declared IBAMR cell traversal; velocity seed axis applies only to the legacy construction. Patch ordering is observable for multiplicative methods and must be exported and compared, not reconstructed from unordered sets after the fact.
- A patch-definition test must compare mapped seed order, overlapping DOF order/content, restricted/owned DOF content where applicable, closure policy, numerical-zero tolerance, and coverage.

#### Pressure-cell seeds versus IBAMR velocity seeds

The M1/M2 inventory has now resolved the main architectural mismatch: cleaned IBAMR enumerates selected velocity-component DOFs and expands the full level `A00`; the formal sandbox starts from every pressure cell and expands only `E_h`. The integration reconstruction must add the pressure-seeded path and supply the live `E_h` separately to construction while continuing to extract local solve matrices from the original full saddle-point operator.

Milestone M2 must still verify:

1. exact extraction and DOF mapping of live `E_h` on every level, including the sign/scale irrelevance of its sparsity pattern;
2. row-or-column adjacency without assuming numerical symmetry;
3. exact pressure-cell logical order, periodic identification, boundary behavior, seed stride, and duplicate handling;
4. exact RELAXED agreement with the pinned targeted-IB family and native STRICT behavior under the declared extension;
5. separation of the `E_h` construction graph from the original full saddle-point local matrix and original-residual update graph.

Final smoother agreement cannot waive exact patch structure. The legacy velocity-seeded result and the pressure-seeded result are distinct named algorithms until structural evidence proves equality for a particular case.

### 4.2 Global multiplicative CAV

For ordered base patches `P_1, ..., P_m`, one sweep is sequential. Starting with the original residual `r_0 = b - A x_0`, patch `j` solves its local system using the current residual restricted to `P_j`, applies the local correction to the global iterate, and updates/recomputes the residual before patch `j+1`. Later patches therefore see all earlier corrections. Optional sweep reversal or symmetric forward/backward composition must be represented as explicit traversal policy, not hidden in container order.

The correctness oracle is the sandbox global multiplicative CAV map after exact IBAMR-to-sandbox DOF mapping. A flat sum of independently computed patch corrections is additive ASM and does not satisfy this contract.

### 4.3 Macro-subdomain CAV-RAS

Partition the ordered base patches into ordered, possibly overlapping macro-subdomains. Inside macro-subdomain `M_q`, execute CAV **multiplicatively**: every base-patch solve sees the residual after preceding base-patch corrections within `M_q`. Form a macro correction from that internal sweep, then restrict it to the macro-subdomain's owned/nonoverlap output region. Across macro-subdomains, compose these restricted macro corrections **additively** against the same outer residual/iterate state, followed by the defined global residual recomputation.

The two composition levels must stay distinct in the API and tests:

1. inner composition: ordered multiplicative base-patch CAV;
2. outer composition: restricted-additive macro-subdomain assembly.

**Existing flat additive-then-restrict patch processing is not sandbox CAV-RAS.** Restricting each independently computed base-patch correction and adding those corrections omits the required multiplicative residual evolution inside each macro-subdomain.

### 4.4 Original residual and pressure gauge

- The residual used to establish or validate a sweep is the residual of the original operator: `r = b - A x`. Any transformed, augmented, Schur, or local-backend residual is diagnostic only and cannot replace it.
- At every parity layer that updates a correction, compare both the correction and a fresh original-residual recomputation. Incremental residual updates may be used internally, but must be checked against fresh recomputation at defined audit boundaries.
- Pressure-nullspace handling is part of the operator contract. Gauge projection/normalization must use the same live IBAMR nullspace and weighting semantics on both sides of a comparison.
- Gauge-equivalent pressure vectors may be mapped to the same declared gauge before comparison; no gauge adjustment may alter velocity, momentum residual, divergence residual, or IB residual.
- End-to-end FGMRES reports the true original reduced-system residual and the exact live decomposition below. Preconditioned or internally reported Krylov norms are supplementary.

The pinned oracle and IBAMR velocity-path Jacobian solve a reduced Eulerian state `x=(u,p)`, so the original reduced residual has two row blocks, not three:

`r = b - (A_Stokes + A_IB)x`, `r_m = R_u r`, and `r_d = R_p r`.

Here `r_m` is the total momentum residual and `r_d` is the divergence residual. IB elasticity is a physical contribution inside the momentum equation; it must not be mislabeled as a fictitious third Eulerian row block. Report its live consistency defect by independently evaluating

`c_IB^mf = R_u[(A_live - A_Stokes)x]`, `c_IB^mat = R_u A_IB^FAC x`, and `r_IB = c_IB^mf - c_IB^mat`,

where `A_live` is the original matrix-free IBAMR Jacobian action and `A_IB^FAC` is the FAC-owned live `d_SAJ_mat` used by construction. Report both IB action norms and the defect norm. This is the mandatory **IB coupling residual** for the reduced implementation. If a later formulation carries an independent Lagrangian position/constraint unknown, report its kinematic residual separately; do not manufacture that state for the eliminated velocity path or substitute an external surrogate.

## 5. Four-layer software decomposition

| Layer | Responsibility | Must not own | Initial review/test seam |
|---|---|---|---|
| 1. Patch construction | RELAXED/STRICT closure maps, seed order, base patches, macro grouping, overlap and restricted/owned index sets | Dense/local factorization; correction application; FAC recursion | Native construction tests plus deterministic export/hash and mapped patch parity |
| 2. Local solution backend | Extract/factor/solve a patch matrix; PETSc, Eigen, and BLAS/LAPACK backend-specific algebra and conditioning/failure policy | Patch discovery/order; additive vs multiplicative semantics; grid transfers | Backend matrix/solve unit tests on live extracted matrices, including singular/ill-conditioned modes |
| 3. Correction composition | Global multiplicative traversal and residual evolution; inner-multiplicative/outer-RAS macro composition; relaxation and reversal policy | Patch discovery; backend-specific algebra; coarse-grid recursion | Exact small-system correction and residual-update tests plus sandbox sweep parity |
| 4. FAC/V-cycle orchestration | Pre/post smoothing, restriction/prolongation, coarse solve, level residual lifecycle, nullspace/gauge propagation | Redefining patch semantics or backend solve math | Transfer parity, smoother-in-FAC tests, V-cycle parity, and end-to-end FGMRES physical residuals |

This separation is an acceptance requirement, not merely a preferred class layout. Construction and application changes land separately unless an unavoidable compile seam is documented in the decision log.

Across all four layers, persistent algebraic state uses the full global velocity-pressure index space. The only representation boundary is the explicitly local patch solve. Shared live matrices are borrowed, not cloned; small derived graph/index metadata may be cached only when it removes repeated setup/apply work and its ownership, invalidation, and memory cost are explicit. No layer may introduce a copied matrix as an implicit convenience cache.

## 6. Branch and recovery strategy

### Immutable inputs

- Historical refs remain untouched.
- All new branches start from the cleaned foundation or the preceding durable review unit.
- Branch 7 and branch 8 are mined for intent, tests, and focused code fragments; they are not wholesale cherry-picked.

### Two green tips with different jobs

1. **Temporary integration reconstruction.** Create a disposable branch from the cleaned foundation and reconstruct enough of historical 7/8 plus current algorithm work to prove current-master end-to-end global multiplicative parity quickly. It may contain deliberately coarse commit boundaries, but it must stay buildable and auditable. It is never merged wholesale.
2. **Durable stacked review branches.** Build a separate stack of small branches, one conceptual claim per review unit, each independently green with native tests. Port only understood changes from the integration reconstruction.

The integration tip and durable-stack tip are both kept green. Once the durable stack reaches the same multiplicative parity, record the equivalence evidence and retire the integration branch without rewriting or deleting historical refs. CAV-RAS begins only after this point.

### M0-confirmed durable review units

The read-only inventory confirms that historical development boundaries are too broad to be presumed pull-request boundaries. The durable stack uses these smaller units; a lettered pair means two independently reviewable tips, not one PR with two commits:

1. `cav-0-recovery-plan`: authoritative plan, canonical ledger, inventory, and baseline commands/status; no solver semantics.
2. `cav-1a-live-operator-audit-schema`: narrow const access/snapshot interfaces and provenance for live operator, block, DOF, nullspace, and gauge data; no serializer or comparator.
3. `cav-1b-raw-export-comparator-contract`: raw serialization, declared mapping, comparator, hand-checkable live case, and comparator mutation tests; no solver behavior change.
4. `cav-1c-matrix-live-kernel-consistency`: one canonical coordinate/stencil/weight/arithmetic contract for every matrix-builder-supported regularized-delta kernel and its live Fortran interpolation/spreading path in 2D/3D, with exhaustive focused cross-path tests; no CAV construction or solver behavior change.
5. `cav-2-patch-construction-audit`: explicitly reconcile pressure-cell and velocity seed formulations; deterministic RELAXED/STRICT base-patch construction/order export and focused native tests; no correction application.
6. `cav-3a-blas-lapack-lu-backend`: current-foundation BLAS/LAPACK shell-backend integration with LU, workspace lifecycle, failure propagation, and focused local-matrix tests.
7. `cav-3b-blas-lapack-robust-modes`: SVD, QR, symmetric-indefinite, and any retained Cholesky compatibility/error policy, with rank/rcond and singular/near-singular tests; no traversal change.
8. `cav-4-multiplicative-composition`: production global multiplicative correction/original-residual semantics with exact per-patch state tests; backend-independent.
9. `cav-5-smoother-parity`: sandbox one-sweep parity on live IBAMR data, including mapping, gauge, original residual, and forward/reverse policy checks.
10. `cav-6a-fac-stage-observability`: narrow FAC/transfer/coarse-stage snapshots and lifecycle checks needed for first-mismatch localization; no FAC algorithm change.
11. `cav-6b-fac-vcycle-fgmres-parity`: production transfer/FAC/V-cycle integration, FGMRES parity, physical residual components, and equivalence to the disposable integration tip. This is the exact M5 durable production tip.
12. `cav-7-macro-construction`: macro geometry, ordered base-patch membership, solve/overlap partitions, and unique restricted/owned output sets only.
13. `cav-8a-macro-local-multiplicative`: a backend-independent local macro sweep in which each base patch sees the updated macro residual; no outer assembly.
14. `cav-8b-outer-restricted-additive`: restrict each completed macro correction and add across macro-subdomains from the common outer residual; no FAC integration.
15. `cav-9a-cav-ras-fac-production`: integrate serial CAV-RAS into the production smoother/FAC/V-cycle path and report original physical residual components.
16. `cav-9b-cav-ras-parity-acceptance`: live sandbox parity, required numerical matrix, production example configuration, and user documentation. This is the exact M9 serial production tip.

These boundaries separate schema from serialization, construction from application, backend integration from optional robust modes, observability from FAC behavior, inner macro traversal from outer assembly, and production RAS integration from acceptance tooling. A later finding may revise them only through a dated plan update.

### Milestone 0 cleaned-stack capability and test inventory

Every exact cleaned tip below is an ancestor of foundation `0c8447e3...`. The exact lettered tips—not similarly named aggregate or backup refs—define the architectural 1.1–6.4 stack.

| Boundary | Exact tip | Capability present at the foundation | Native coverage introduced or carried by that boundary | Recovery consequence |
|---|---|---|---|---|
| 1.1 | `b4c847e98d1fc5561588f068524031b54f14a2c3` | Implicit-IB support APIs in PETSc matrix/level-solver, index, and IB strategy surfaces | No focused promoted native test; ledger finding 1 preserves this as an intentional cleanup-pass deferral | Do not silently claim direct API coverage; add only alongside a review unit that needs a changed support API |
| 1.2 | `3fded45fce94db32e86e820dd4cf13e8d15f06f6` | Stokes-IB nonlinear operator plus shared time/operator plumbing | `implicit_stokes_ib_solver_components_01`, one-level case | Reuse private shared plumbing; do not duplicate it in parity hooks |
| 1.3 | `d4a7363724af6b8d2f9be2bad1253b5f036d87ba` | Stokes-IB Jacobian operator using the shared schedule/plumbing | Same focused one-level solver-components test, extended for Jacobian behavior | Audit/export interfaces must observe the live Jacobian rather than reconstruct it |
| 1.4 | `d65dd80feebc6f3461b42e6caaf45c4932dcb3e0` | Stokes-IB FAC level operator and centralized residual-vector lifetime | Same focused solver-components test | Preserve the ledger-resolved residual-vector deduplication |
| 1.5 | `f50e50ede2186f56cfc1ee7a9700ee62fb4c540a` | Stokes-IB Jacobian FAC preconditioner | `implicit_stokes_ib_solver_components_01` at one, two, and three levels | Supplies the live multilevel orchestration seam for `cav-6a/6b` |
| 2 | `adbbe69946cbfd74811f9f8d28784b642a7a22cb` | Integrator refactored onto solver/operator components and shared time schedule; implicit ex0 | `implicit_stokes_ib_01` plus the solver-components cases | End-to-end recovery must use this production operator path, not an external driver surrogate |
| 3.1 | `52dd7d9e44d9da76fb0c2b4349ea59a4f77a4301` | Linear cell-pressure prolongation and shared side/row construction | 2D/3D linear pressure prolongation parity, including nonlinear profile cases | Reuse the live transfer builder and its two-phase PETSc row construction |
| 3.2 | `1bfbfe3b3788fd86efafe5ae4f36b3531b558b18` | Linear cell-pressure restriction/coarsening | 2D/3D linear pressure restriction parity, including nonlinear profile cases | Supplies live restriction for multilevel first-mismatch checks |
| 4 | `70f2fae95527a3eea5c2bc5876f44b61a579c538` | Maintained transfer-parity surface and consolidated helpers | RIDP sums, RT0 prolongation/restriction equivalence, RT0 refine affine/divergence checks, pressure transfer, and coarse `S A J`, in 2D/3D where supplied | Keep parity helpers consolidated but geometry-specific data filling separate |
| 5 | `2e77def2046901ef628e6e634ee303c94826830d` | Coupling-aware map cache; velocity-seeded RELAXED/STRICT construction; seed order/stride controls; implicit-ex0 selection | Sparse closure, patch-level maps, reference patch parity, ASM/CAV modes, 2D strict/classical equivalence, 3D strict equivalence, and CAV shell semantics | Construction exists, but pressure-seed equivalence and exact ordered export remain open M2 obligations; mutable cache internals stay private |
| 6.1 | `efe5bc430dfccf9d4d6a9d714faed4b57b9f84cb` | PETSc shell-backend framework, parsed backend state, RHS boundary adjustment, option precedence, and deduplicated velocity-block embedding | PETSc/default shell semantics, option override, and RHS-adjustment tests | New backends plug into this seam; do not fork traversal or matrix assembly |
| 6.2 | `b222a691fd3bf38067f6b473326c9fd6f3eb8ff0` | Eigen reference/factorized local-solve framework, shared serial guard and sweep logic | Eigen reference, all-local-solvers exercise, option override, and additive-mode rejection | Canonical local/sweep reference for recovery and backend cross-checks |
| 6.3 | `aab6a143cf22406278eb6f021fe347db642f90d6` | Optimized Eigen factorization and pseudoinverse backends with shared dense dispatch | Eigen and Eigen-pseudoinverse shell tests | Preserve one shared composition path across dense solvers |
| 6.4 | `0c8447e3593113a8b62f94df9482c60d896f0d3b` | Stokes Eigen Schur-complement backend reusing the shared sweep; cleanup formatting | Eigen-Schur shell case and expected rank-2 unsupported-mode failure | Exact foundation and production comparison backend; still lacks current-oracle end-to-end parity evidence |

### Milestone 0 historical 7/8 inventory

| Historical commit | Capability/evidence | Native or parity coverage | Durable recovery unit(s) | Disposition |
|---|---|---|---|---|
| `4de0cd9fd` | Initial BLAS/LAPACK shell backend and additive/multiplicative/restrict application modes | Baseline BLAS/LAPACK shell semantics | `cav-3a`; composition intent only for `cav-4`/later RAS | Mine backend integration by intent; do not recover a flat restrict mode as CAV-RAS |
| `4f83a0280` | LU, SVD, QR, symmetric-indefinite, Cholesky path, rcond/rank policy, reusable workspaces, and active residual-update rows | LU/QR/SVD, near-singular QR/SVD, symmetric-indefinite, Cholesky behavior | `cav-3a`, then `cav-3b`; residual-row ideas reviewed separately in `cav-4` | Split core backend from robust modes and composition; do not cherry-pick the 979-line change wholesale |
| `ad3dcdb13` | Large live subdomain-construction and operator-chain audit tests plus production visibility hooks | `implicit_stokes_ib_cav_subdomain_construction_parity_01`; multilevel `implicit_stokes_ib_operator_chain_parity_01` | `cav-1a`, `cav-2`, `cav-6a` | Recover only narrow const observability and focused layer tests; avoid re-exposing mutable internals |
| `c9ce85d2f` | Trimmed the maintained audit surface | Audit-surface reduction | `cav-1a`, `cav-6a` | Evidence for deletion-oriented boundaries, not code to restore automatically |
| `6f9037942` | Live ex0 export/comparison workflow, FAC-stage hooks, and IBAMR-only preconditioner residual sweeps | Comparator/driver plus `K=1e2,1e4,1e6,1e8` residual cases | `cav-1b`, `cav-5`, `cav-6a`, `cav-6b` | Split raw export/comparator from smoother and FAC/FGMRES claims; rerun from current live objects |
| `1529bca5b` through `d400a07eb` | Strict-profile checks, restored fixtures, FAC assertion correction, tightened sign/gauge/order semantics, and conditioning-aware high-K handling | Historical layerwise and residual evidence | Same focused units as above | Preserve intent, but low-K fixed targets and current pinned oracle supersede old logs; `K=1e8` remains diagnostic |

Historical branch 8 contains no macro-subdomain inner-multiplicative/outer-restricted-additive implementation. RAS therefore comes only from the pinned formal/oracle sources after M5, not from relabeling its flat restricted patch path.

### Ref topology and exclusions recorded by M0

- The required historical local refs are exactly branch 7 `4f83a028...` and branch 8 `d400a07e...`. Their currently observed remote-tracking tips, `d54c4941...` and `a2ae04a1...`, differ and are not substitutes.
- The required historical refs share merge base `bd010f326eebdb6b23dbde5c2779789c6586de9f` with the foundation and are not its descendants.
- Aggregate local refs for branch 1 (`3228efc9...`), branch 2 parity (`89e16311...`), branch 3 (`626150c1...`), and branch 6 (`c6c661be...`) are not ancestors of `0c8447e3...`; the aggregate remote branch-3 (`d2ba7101...`) and branch-6 (`53cc911c...`) tips also differ. They are excluded as architectural boundaries. The exact 1.1–6.4 tips in the table above are the cleaned ancestry.
- Backup refs remain historical evidence only. No local, remote-tracking, backup, tag, or sandbox ref was changed by this inventory.

### Milestone 0 conclusions

1. The foundation already contains all four architectural layers—velocity-seeded patch construction, several local backends, one shared global multiplicative sweep mechanism, and FAC/V-cycle orchestration—but it does not contain current pinned-oracle end-to-end evidence. Existing native CAV/shell tests establish local semantics, not sandbox reproduction of the live ex0 operator chain.
2. BLAS/LAPACK is a backend recovery, not the definition of multiplicative composition. Core LU integration and robust factorization modes are independently testable and therefore remain `cav-3a` and `cav-3b`.
3. Historical branch 8's operator-chain test, example exporter, comparator, FAC hooks, and residual sweeps span several conceptual claims and thousands of lines. They must be recovered through `cav-1a/1b`, `cav-2`, `cav-5`, and `cav-6a/6b`, with narrow production visibility instead of restoring the entire audit surface.
4. The historical IBAMR-only stiffness residual sweep is a production regression guard, not cross-code numerical parity. Its old logs and thresholds cannot validate the current candidate.
5. No historical 7/8 capability supplies true macro-subdomain CAV-RAS. The durable RAS path begins only after the M5 all-lane gate and keeps construction (`cav-7`), inner multiplicative application (`cav-8a`), outer restricted-additive assembly (`cav-8b`), FAC integration (`cav-9a`), and acceptance (`cav-9b`) separate.

### Foundation baseline commands and status

The reproducible Debug baseline uses a build tree outside the repository:

```sh
/Users/boyceg/code/autoibamr/dbg/packages/cmake-3.30.6/bin/cmake \
  -S /Users/boyceg/.codex/worktrees/0760/IBAMR \
  -B /tmp/ibamr-cav-foundation-0c8447e-dbg \
  -DCMAKE_BUILD_TYPE=Debug \
  -DLIBMESH_ROOT=/Users/boyceg/code/autoibamr/dbg/packages/libmesh-1.7.8 \
  -DLIBMESH_METHOD=devel \
  -DPETSC_ROOT=/Users/boyceg/code/autoibamr/dbg/packages/petsc-3.23.3 \
  -DHYPRE_ROOT=/Users/boyceg/code/autoibamr/dbg/packages/petsc-3.23.3 \
  -DHDF5_ROOT=/Users/boyceg/code/autoibamr/dbg/packages/hdf5-1.12.2 \
  -DSAMRAI_ROOT=/Users/boyceg/sfw/samrai/IBSAMRAI2/darwin-clang-dbg \
  -DSILO_ROOT=/Users/boyceg/code/autoibamr/dbg/packages/silo-4.11-bsd \
  -DNUMDIFF_ROOT=/Users/boyceg/code/autoibamr/dbg/packages/numdiff-5.9.0
```

Build the named focused targets for solver components, implicit IB, pressure/RT0/RIDP transfers, CAV construction/modes, and PETSc/Eigen/Schur shell semantics:

```sh
/Users/boyceg/code/autoibamr/dbg/packages/cmake-3.30.6/bin/cmake \
  --build /tmp/ibamr-cav-foundation-0c8447e-dbg --target \
  tests-IB_implicit_stokes_ib_01 \
  tests-IB_implicit_stokes_ib_coarse_saj_transfer_parity_01 \
  tests-IB_implicit_stokes_ib_solver_components_01 \
  tests-navier_stokes_stokes_ib_cav_patch_level_maps_2d \
  tests-navier_stokes_stokes_ib_cav_reference_parity_2d \
  tests-navier_stokes_stokes_ib_cav_sparse_map_closure_2d \
  tests-navier_stokes_stokes_petsc_level_solver_asm_modes_2d \
  tests-navier_stokes_stokes_petsc_level_solver_asm_modes_3d \
  tests-navier_stokes_stokes_petsc_level_solver_cav_shell_semantics_2d \
  tests-navier_stokes_stokes_petsc_level_solver_shell_multiplicative_semantics_2d \
  tests-navier_stokes_stokes_petsc_mat_utilities_ops_2d \
  tests-navier_stokes_stokes_petsc_mat_utilities_pressure_transfer_2d \
  tests-navier_stokes_stokes_petsc_mat_utilities_pressure_transfer_3d \
  tests-navier_stokes_stokes_petsc_mat_utilities_prolongation_2d \
  tests-navier_stokes_stokes_petsc_mat_utilities_prolongation_3d \
  tests-navier_stokes_stokes_petsc_mat_utilities_restriction_2d \
  tests-navier_stokes_stokes_petsc_mat_utilities_restriction_3d \
  tests-navier_stokes_stokes_petsc_mat_utilities_ridp_sums_2d \
  tests-navier_stokes_stokes_petsc_mat_utilities_ridp_sums_3d \
  tests-refine_rt0_refine_01_2d \
  tests-refine_rt0_refine_01_3d \
  -j8
```

The legacy RT0-refine tests need their input/output symlinks beside the built executables before `attest` discovers them:

```sh
bash /Users/boyceg/.codex/worktrees/0760/IBAMR/tests/link-test-files.sh \
  /Users/boyceg/.codex/worktrees/0760/IBAMR/tests/refine \
  /tmp/ibamr-cav-foundation-0c8447e-dbg/tests/refine
```

Then run from `/tmp/ibamr-cav-foundation-0c8447e-dbg`:

```sh
/Users/boyceg/.codex/worktrees/0760/IBAMR/attest -j8 --verbose \
  -R '^(IB/implicit_stokes_ib_(01|solver_components_01|coarse_saj_transfer_parity_01)|navier_stokes/stokes_(ib_cav|petsc_level_solver|petsc_mat_utilities)|refine/rt0_refine_01)'
```

The first 24-test solver/CAV-focused subset passed 24/24. After building the additional 2D/3D pressure-transfer, prolongation/restriction, RIDP, RT0-refine, ASM-mode, and implicit-IB targets, the expanded exact-foundation selection passed **69/69**. The build emitted only the already-visible libMesh override and duplicate-library/rpath warnings; no target or test failed. These are implementation-task baseline results only and cannot satisfy an independent-audit gate.

### Durable `cav-1a` implementation-validation evidence

The uncommitted `cav-1a-live-operator-audit-schema` candidate is based on plan tip `89b6da10e259342f6383c0caa2d619a2e72f0a58` in `/tmp/ibamr-cav-durable-1a`, with Debug build `/tmp/ibamr-cav-durable-1a-dbg`. The focused executable target and the normal directory target both build successfully:

```sh
/Users/boyceg/code/autoibamr/dbg/packages/cmake-3.30.6/bin/cmake \
  --build /tmp/ibamr-cav-durable-1a-dbg \
  --target tests-navier_stokes_stokes_petsc_level_solver_shell_multiplicative_semantics_2d -j8
/Users/boyceg/code/autoibamr/dbg/packages/cmake-3.30.6/bin/cmake \
  --build /tmp/ibamr-cav-durable-1a-dbg --target tests-navier_stokes -j8
```

The directory target is required because its standard post-build step links the input/output fixtures into the build test tree. The build compiles the changed level solver in both 2D and 3D and emits only the already-known libMesh override and duplicate-library/rpath warnings.

The focused discovery selector finds exactly `navier_stokes/stokes_petsc_level_solver_shell_multiplicative_semantics_2d.live_operator_state.input`:

```sh
/tmp/ibamr-cav-durable-1a/attest \
  --mpi-executable /opt/homebrew/bin/mpiexec \
  --numdiff-executable /opt/homebrew/bin/numdiff --verbose -N \
  -R 'navier_stokes/stokes_petsc_level_solver_shell_multiplicative_semantics_2d\.live_operator_state\.input'
```

The same command without `-N` passes **1/1** twice in 1.96 s and 0.48 s, below the 15 s focused budget. Its same-unit input/output fixtures are `tests/navier_stokes/stokes_petsc_level_solver_shell_multiplicative_semantics_2d.live_operator_state.{input,output}`; the expected result is exit status zero with exact reported operator identity, global ownership/block partition, provenance, nullspace/gauge, and preinitialize/post-deallocate lifecycle checks. A disposable build-tree mutation changing only `verify_live_operator_state_view` from `TRUE` to `FALSE` leaves the solver execution path unchanged but is rejected by `attest` because all live-view result fields disappear. The unmutated rerun remains green. The mutation is retained only under `/tmp` and is not a candidate artifact.

The canonical broad discovery command finds these 18 cases:

- `IB/explicit_ex0.quartersquare3.input`
- `IB/explicit_ex0.expect_error=true.input`
- `IB/implicit_stokes_ib_coarse_saj_transfer_parity_01.multilevel_l2.input`
- `IB/implicit_stokes_ib_01.input`
- `IB/ib_body_force.input`
- `IB/explicit_ex0.quartersquare2.input`
- `IB/explicit_ex1.input`
- `IB/implicit_stokes_ib_solver_components_01.max_levels=2.input`
- `IB/ib_body_force_kirchhoff.redundant.input`
- `IB/implicit_stokes_ib_solver_components_01.max_levels=1.input`
- `IB/implicit_stokes_ib_solver_components_01.max_levels=3.input`
- `IB/explicit_ex0.input`
- `IB/explicit_ex0.quartersquare.input`
- `IB/ib_body_force_kirchhoff.standard.input`
- `IB/explicit_ex0.expect_error=true.mpirun=4.input`
- `IB/explicit_ex0.quartersquare.mpirun=4.input`
- `IB/explicit_ex0.restart=300.mpirun=4.input`
- `IB/explicit_ex0.quartersquare2.mpirun=4.input`

After `tests-IB` builds and links all required fixtures, the exact broad command is:

```sh
/tmp/ibamr-cav-durable-1a/attest \
  --mpi-executable /opt/homebrew/bin/mpiexec \
  --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R IB/
```

It runs with no exclusions in 82.93 s: all 14 serial cases pass, while all four `mpirun=4` cases fail before entering IBAMR because the managed sandbox rejects Open MPI's TCP listener. Thus the exact broad result is **14/18 with four explicitly environment-blocked rank-4 cases**, not a broad pass. The supplemental command with `-E 'mpirun=4'` confirms the serial scope **14/14** in 80.75 s; this labeled subset does not replace the broad result. Every discovered case expects overall `attest` success under its filename-declared normal or `expect_error=true` semantics. Rank-4 evidence must be rerun in a suitable environment before full review readiness. These results are implementation validation only, not independent audit evidence.

## 7. Reviewability rules

1. One conceptual claim per pull request/review unit.
2. Every durable review unit that introduces or materially changes production functionality includes its corresponding focused native regression tests in the same commit/PR; a production-code-now/tests-later boundary is forbidden.
3. Mechanical refactors, interface cleanup, and test-helper consolidation are separate from changes to patch definitions, solve behavior, traversal, residual updates, or FAC orchestration.
4. Patch construction is separate from patch application.
5. Local-backend selection is separate from correction-composition policy.
6. A review unit must be understandable from its diff and tests without requiring the disposable integration branch.
7. Preserve the historical ledger's implemented refactors and deferred disposition.
8. Commit no build directories, generated matrices/vectors, run logs, or comparison output unless explicitly requested as a curated fixture.
9. Each review unit records its native command, result, and parity layer in this plan before advancing.
10. No final numerical claim relies solely on a test's exit status; record the compared quantity, tolerance, operator/gauge semantics, and original physical residual components.

### Native CI and comment/style contract

The same-unit native tests are deterministic, require no network or external repository, use no large generated fixture, are discoverable and runnable through the normal repository-root `attest` workflow with their required input/output fixtures and expected-result/exit semantics, and are fast enough for normal IBAMR CI. CMake may build the executable and arrange fixtures according to neighboring repository patterns, but a CMake target alone is not CI-integration evidence. External oracle/parity runs remain required where assigned but are heavier reproducibility evidence and never substitute for a native `attest` regression. Reuse a common test executable with separate input/output or expected-result fixtures when cases exercise the same linkage and initialization path. A separate executable requires a documented dimensional, linkage, initialization, expected-failure, or other technical seam. Audit hooks stay narrow and may not turn production APIs into a testing framework. This rule applies prospectively to durable CAV units and does not reopen ledger Finding 1 unless a unit actually introduces or materially changes those support APIs.

Use the checkout's actual repository-root runner from the configured build directory. The canonical broad IB regression form is:

```sh
/path/to/attest \
  --mpi-executable /opt/homebrew/bin/mpiexec \
  --numdiff-executable /opt/homebrew/bin/numdiff \
  --verbose -R IB/
```

`-R IB/` selects all discovered examples/tests under the `IB/` subtree. Development first uses the narrowest `-R '<feature-case-expression>'` that selects the new cases, but a focused pass never substitutes for the broad `-R IB/` regression required before review-unit readiness. The broad candidate-checkout regression has a provisional 30-minute reference Debug budget; any unavailable case, environmental exclusion, expected-failure behavior, or budget exceedance is recorded explicitly and cannot be silently dropped. Preserve the environment's explicit `--mpi-executable` and `--numdiff-executable` paths in every reproducibility command where required.

Before a unit is committed, its recovery-manifest entry records both focused and broad scope: selector, checkout-root `attest` path, exact command, discovered native case names, required input/output fixtures, expected result/exit semantics, expected budget, measured local Debug wall time, deterministic rerun, sensitivity evidence, result, and exclusions. The provisional added-test budgets below are per review unit on the reference local Debug build; a unit that exceeds its focused budget must be split or document why the runtime is technically necessary. Cross-repository parity/audit time is reported separately.

| Review unit | Native CI form | Added-test budget |
|---|---|---:|
| `cav-1a` | Reused level-solver executable; live-view lifecycle fixture | 15 s |
| `cav-1b` | Reused exporter/comparator test driver; tiny raw/mutation fixtures | 30 s |
| `cav-1c` | One kernel-consistency source built in 2D/3D | 120 s total |
| `cav-2` | Reused patch-construction executable; RELAXED/STRICT and order fixtures | 90 s |
| `cav-3a` | Reused shell-backend executable; LU/lifecycle/failure fixtures | 45 s |
| `cav-3b` | Same backend executable; robust-mode fixtures | 90 s |
| `cav-4` | Reused shell-composition executable; traversal/residual fixtures | 60 s |
| `cav-5` | Reused live-smoother driver; forward/reverse/profile fixtures | 90 s |
| `cav-6a` | Reused FAC executable; stage/lifecycle/regrid fixtures | 120 s |
| `cav-6b` | Reused implicit-IB executable; V-cycle/FGMRES/residual fixtures | 180 s |
| `cav-7` | Reused construction executable; macro ownership/order fixtures | 90 s |
| `cav-8a` | Reused composition executable; macro-local multiplicative fixtures | 90 s |
| `cav-8b` | Same composition executable; outer restriction/addition fixtures | 90 s |
| `cav-9a` | Reused implicit-IB/FAC executable; serial production fixtures | 180 s |
| `cav-9b` | Same production executable plus compact acceptance fixtures | 300 s |

New and materially changed code follows the surrounding IBAMR Doxygen and inline-comment style. Comments explain only non-obvious mathematical intent, invariants, ordering, gauge/nullspace behavior, ownership/lifetime/invalidation, supported scope, failure behavior, and performance-sensitive representation choices. They explain why a choice is required rather than narrating code; avoid generated-looking essays, redundant commentary, speculative frameworks, and text that can drift from implementation. Each new file/class/function is checked against nearby IBAMR examples, and public configuration/failure behavior is documented consistently.

The pre-commit review-unit checklist is therefore: one claim; feature and focused test together; focused `attest` discovery/execution with fixtures and expected status; deterministic sensitivity; appropriate broad candidate-checkout `-R IB/` regression within its budget or explicit blocking/exclusion record; no focused-for-broad substitution; no unnecessary executable/helper/fixture; comments and public API style checked against neighboring code; exact selectors/commands/results/runtimes recorded; no generated artifact committed.

## 8. Verification and parity ladder

Advance only when the preceding layer is structurally and numerically green.

| Layer | Required comparison/evidence | Gate |
|---|---|---|
| 1. Problem/operator | Grid, BCs, coefficients, time/scale factors, DOF maps, `A`, `b`, nullspaces/gauge, IB coupling blocks from live objects | Exact metadata/structure after mapping; raw cross-runtime numerical differences reported and localized; fixed numerical targets evaluated in the declared common-arithmetic replay |
| 2. Patch definitions/order | Seed order; RELAXED/STRICT base patches; overlap and restriction/ownership sets; macro membership/order | Exact mapped structural agreement; deterministic repeatability |
| 3. Local matrices | Extracted local matrices/RHS, row/column order, backend solution, factorization/failure mode | Exact structure/order; common-arithmetic solution parity within the fixed target; live backend backward error and forward differences reported separately |
| 4. Corrections/residual updates | Per-patch correction, iterate after each stage, incremental residual, fresh `b-Ax`, gauge projection | Common-arithmetic stagewise low-stiffness target; live incremental/fresh residual consistency and first-mismatch trace |
| 5. Transfers | Live restriction, prolongation, redistribution/synchronization, pressure transfer and gauge behavior | Existing native transfer tests plus exact mapped structure and low-stiffness numerical target |
| 6. Smoother | One forward sweep and each supported traversal/composition mode | Common-arithmetic original-residual and correction parity; live production trace retained; no comparison based only on internal norm |
| 7. V-cycle | Pre/post sweeps, coarse RHS/solve/correction, level residuals, final correction | Common-arithmetic stagewise parity plus live original-residual reduction |
| 8. FGMRES | Original operator solve with production preconditioner | Common-arithmetic final-output limits plus separate live momentum/divergence/IB residual breakdown, convergence status, and iteration history |

Audit exports must include provenance: Git commit, build identity, dimension, MPI size, input/configuration, matrix/vector dimensions, DOF mapping version, gauge declaration, stiffness `K`, and floating-point precision.

## 9. Acceptance criteria

### Structural acceptance

- After a declared IBAMR-to-sandbox mapping, problem metadata, operator sparsity, DOF ownership/order, patch definitions/order, local row/column order, transfer structure, macro membership, and restriction/ownership sets agree exactly.
- Missing, duplicated, reordered, or silently dropped patches/DOFs are failures even if a final norm happens to be small.

### Numerical acceptance

Numerical acceptance has two mandatory tracks with different purposes. Neither may substitute for the other.

1. **Common-arithmetic algorithm parity.** Independently exported candidate and oracle metadata, DOF maps, operator structure, patch/macro definitions and order, transfer definitions, and traversal policy are mapped into one canonical ordering. A single documented replay implementation then executes both declared algorithms with the same mapped operator/RHS, local-solve routine, sparse/dense multiplication order, gauge operation, and floating-point runtime. Run that paired replay twice: once on the mapped live-IBAMR operator/RHS and once on the independently exported oracle operator/RHS. The candidate replay must consume the candidate's exported ordered patch/macro/transfer artifacts rather than silently rebuilding them with oracle construction code. The oracle replay consumes the independently exported oracle artifacts. The fixed `K=1,1e2,1e4` table below applies to corrections, residual updates, smoother/V-cycle stages, solution, and final original residual in both common-input runs.
2. **Live production behavior.** The actual IBAMR matrix-free/operator, recovered local backend, multiplicative shell, FAC/V-cycle, and FGMRES path must be executed. Report convergence reason, iteration count/history, fresh original residual, local backward errors, incremental-versus-fresh residual consistency, and the absolute/scaled momentum, divergence, and IB coupling residual components defined below. For the mandatory low-`K` cases the independently recomputed live original relative residual must be at most the declared `1e-10` solve tolerance; the scalar vector `(rho, ||r_m||_2, ||r_d||_2)` and iteration history must be compared with the oracle run under the fixed per-`K` table. The live-only coupling consistency check must report both IB action norms and satisfy `rho_IB <= 1e-10`; it is not compared to a fictitious oracle third row, since the oracle assembles that contribution directly. Raw entrywise live-versus-MATLAB residual-vector differences remain mandatory diagnostics with first-mismatch localization, but vendor/runtime-dependent entrywise identity is not the algorithm-parity gate.

The common replay is not permission to replace live construction with a surrogate. Exact mapped structure and order, Fortran-authoritative matrix/live kernel consistency, raw independently constructed operator exports, and comparator/mapping mutation controls remain mandatory. Any raw numerical difference must be reported with its earliest source; it cannot be deleted, normalized away, or described as common-arithmetic evidence. Conditioning-aware thresholds remain forbidden at `K <= 1e4` in the common-arithmetic track.

After exact structural mapping, common-arithmetic vector and canonicalized-matrix-entry parity uses

`E_inf(x,y) = ||x-y||_inf / max(1.0, ||x||_inf, ||y||_inf)`.

The denominator floor is exactly `1.0`, so sub-unit quantities are checked by absolute max error. Always report the unscaled absolute max error too. Local solve backward error is

`eta_i = ||A_i delta_i-r_i||_inf / (||A_i||_inf ||delta_i||_inf + ||r_i||_inf + 1e-30)`,

and incremental/fresh residual consistency is

`E_r = ||r_inc-(b-Ax)||_inf / max(1.0, ||r_inc||_inf, ||b-Ax||_inf)`.

Original residual reporting uses `rho = ||b-Ax||_2/max(||b||_2,1e-30)`, `rho_m = ||r_m||_2/max(||b_m||_2,sqrt(eps)||b||_2,1e-30)`, and `rho_d = ||r_d||_2/max(||b_d||_2,sqrt(eps)||b||_2,1e-30)`. The IB coupling consistency metric is `rho_IB = ||r_IB||_2/max(||c_IB^mf||_2,||c_IB^mat||_2,sqrt(eps)||b||_2,1e-30)`. Absolute total/component/action/defect norms are mandatory so denominator floors cannot hide a defect.

| Stiffness | Fixed target | Gate rule |
|---:|---:|---|
| `K=1` | `E_inf <= 1e-13` | Near-roundoff expectation; a larger reproduced discrepancy blocks acceptance. |
| `K=1e2` | `E_inf <= 1e-12` | Fixed low-stiffness target. |
| `K=1e4` | `E_inf <= 2e-11` | Fixed `1e-11`-scale target; every final-output comparison also has the `1e-10` hard ceiling. |
| `K=1e6` | Conditioning-aware | Allowed only with measured condition estimate and backward error; never relax structure, gauge, or original-residual reporting. |
| `K=1e8` | Diagnostic | Record behavior and failure mode; no gate until a dated decision adds one. |

The `1e-10` final ceiling is not permission to waive a per-K common-arithmetic target. A common-arithmetic target miss must be localized and corrected before an implementation milestone is declared green. A raw cross-runtime entrywise miss is instead a mandatory localized diagnostic; it blocks only if it reveals a structural/formula/semantic defect, invalid mapping/comparator, failed live convergence, or a violation of the separate live-production criteria.

### Mandatory numerical matrix

The exact matrix, including modes, traversal, backends, dimensions, supported MPI scope, formulas, and outcomes, is normative in `plans/implicit-ib/cav-audit-numerical-prompt.md`. Its required rows are summarized here:

- `N0`: live `4x4` 2D rank-1 hand mapping/comparator case at `K=1`;
- `N1`: single-level `8x8`, RELAXED/STRICT, forward/reverse, all low K, Eigen reference;
- `N2`: two-level `8x8 -> 16x16`, forward/reverse/palindrome, all low K, through FGMRES;
- `N3`: three-level `8x8 -> 16x16 -> 32x32` production multiplicative path, all low K;
- `N4`: backend cross-check on `8x8` at `K=1,1e4` for Eigen reference, production Eigen/Schur, and recovered BLAS/LAPACK LU/SVD;
- `R0`–`R2`: serial CAV-RAS from single-level macro traces through the three-level M9 production path, with fixed macro blocks and ghost widths;
- `P2`: native 2D rank-2 ownership/lifecycle/rank-consistency scope;
- `P3`: native 3D Debug/Release rank-1 build/smoke/lifecycle scope.

Full sandbox parity is required in supported 2D rank-1 scope. Rank-2 and 3D checks are native production-robustness evidence unless broader parity is explicitly advertised. The M9 gate certifies serial CAV-RAS only; distributed RAS has a later gate.

### Export/mapping/comparator acceptance

Parity cannot pass until raw pre-mapping exports are retained and the mapping/comparator are independently checked on a live `4x4` case. Required controls mutate disposable raw copies with DOF permutations, velocity/pressure sign changes, legal pressure-only gauge shifts, illegal velocity shifts, omissions, and patch reorderings. Declared maps and legal pressure gauge normalization must pass; undeclared permutations, sign errors, omissions, illegal shifts, and meaningful reorderings must be detected. Record raw/mapped/comparator hashes and an unmutated control.

## 10. Milestone checklist

- [x] **M0 — Inventory and review decomposition:** audit cleaned 1.1–6.4, historical 7/8, relevant native tests, build topology, and establish baseline commands/status. Implementation validation: exact-foundation Debug build and focused selection passed 69/69; this is not independent-audit evidence.
- [ ] **M1 — Live parity harness foundation:** establish all-kernel matrix/live interpolation-and-spreading consistency, then export/map problem and operator data from live IBAMR objects.
- [ ] **M2 — Base-patch parity:** exact RELAXED/STRICT definitions and order after mapping.
- [ ] **M3 — Local backend recovery:** robust local-solve capability with focused native tests and live local-matrix parity.
- [ ] **M4 — Global multiplicative CAV:** exact correction/residual semantics and one-sweep sandbox parity.
- [ ] **M5 — Multiplicative production path:** exact durable `cav-6b-fac-vcycle-fgmres-parity` tip has the two-source common-arithmetic transfer/FAC/V-cycle/FGMRES matrix green at fixed low-`K` tolerances, has the separate live production convergence and momentum/divergence/IB residual checks green, and matches the disposable integration reconstruction. The disposable reconstruction now satisfies this numerical contract; the durable tip and its independent gates remain open.
- [ ] **M6 — Retire integration scaffolding:** only after the M5/cav-6b all-six independent gate passes, record two-tip equivalence and stop using the disposable integration branch.
- [ ] **M7 — Macro-subdomain construction:** deterministic grouping and restriction/ownership tests.
- [ ] **M8 — Macro-subdomain CAV-RAS:** inner multiplicative/outer restricted-additive semantics and sandbox parity.
- [ ] **M9 — Production acceptance:** exact `cav-9b-cav-ras-parity-acceptance` tip has the serial CAV-RAS low-stiffness matrix green, high-stiffness diagnostics recorded, documentation/example finalized, and generated artifacts excluded.
- [ ] **Audit-suite readiness (audit-coordinator/auditor-owned):** all six lane sensitivity controls have activated mutations or structured challenges, unmutated controls, immutable evidence, and the required detectors.
- [ ] **Working multiplicative-CAV audit gate (auditor/synthesis-owned):** all six independent lane reports plus synthesis are frozen and hashed for the exact M5 `cav-6b` durable production tip.
- [ ] **Working serial-CAV-RAS audit gate (auditor/synthesis-owned):** all six independent lane reports plus synthesis are frozen and hashed for the exact M9 `cav-9b` production tip.

The implementation task may prepare candidates and evidence but must not mark the three audit-owned checklist items complete. The separate audit coordinator owns coverage/custody, lane auditors own their outcomes, and the synthesizer owns only the frozen-suite milestone recommendation.

## 11. Recovery manifest

| Capability | Historical source | New review unit | Native tests | Parity layer | Status |
|---|---|---|---|---|---|
| Authoritative recovery control, canonical ledger, and baseline | Cleaned 6.4 plus original handoff ledger | `cav-0-recovery-plan` | Documentation cross-check; exact-foundation Debug focused build/test baseline 69/69 | Provenance/pre-layer | **M0 complete as implementation validation; independent audit not launched** |
| Nonowning live operator/block/DOF/gauge view | Historical 8 `ad3dcdb13` and trimmed intent at `c9ce85d2f` | `cav-1a-live-operator-audit-schema` | Reused 2D level-solver executable with same-unit live-view fixture; 2D/3D library build; focused `attest` 1/1 twice; hook-disabled sensitivity detected; broad `-R IB/` 14 serial passes plus four rank-4 environment blocks, with supplemental serial 14/14 | Problem/operator | Candidate implementation validated locally; exact base `89b6da10e`; production commit pending; rank-4 broad rerun required for full review readiness; independent audit not launched |
| Raw exporter, declared mapping, and comparator contract | Historical 8 `6f9037942` through `553454567`, recovered by intent | `cav-1b-raw-export-comparator-contract` | Live `4x4` hand case, raw round-trip, exact control, and permutation/sign/gauge/omission/order mutations | Problem/operator; mapping | Proposed after `cav-1a` |
| All-kernel matrix/live interpolation and spreading consistency | Original live Fortran interaction kernels in `lagrangian_interaction2d/3d.f.m4` and included `lagrangian_delta.f.m4`; cleaned PETSc matrix builder follows them | `cav-1c-matrix-live-kernel-consistency` | Tiny periodic 2D/3D cases for all 13 shared scalar or component/transverse choices: matrix-row basis interpolation, transpose-scaled spreading, constants, generic, half-cell, lower-boundary, and upper-boundary coordinates; live periodic representatives are folded into canonical global DOFs; captured-row formula and weight-order mutations must be detected | Live operator factors and actions | Disposable integration complete and green: baseline exposed incorrect matrix IB5/IB6 weights (`1.24e-1`/`9.54e-2` in 2D) and a missing inner-loop x-index update in live 3D IB5 spreading (`1.64e-1`). Matrix formulas were ported from the authoritative Fortran paths; the sole live edit fixes the independently reproduced 3D IB5 stale-index bug. Formula/orientation/index fixes pass in both dimensions, including actual wrapped stencils, with maximum shared-factor differences `1.30e-13` (2D) and `7.23e-14` (3D); every kernel detects the formula and weight-order controls. The complete pre-existing serial interpolation/spreading suite is green **126/126** after replacing the one 3D IB5 spreading fixture that encoded the corrected index bug; the replacement is independently justified by the matrix-transpose check. Direct full-kernel basis probing and a shared weight-ABI refactor were rejected for this recovery because they add setup copies/calls or modify the correct matrix-free hot path. Focused implicit-IB construction/smoother regressions are green **7/7**. Fresh `K=1`, `1e2`, and `1e4` end-to-end exports are byte-for-byte unchanged from the pre-correction IB4 operator, smoother, V-cycle, solution, and original-residual vectors. Durable port remains open |
| Pressure-seeded `E_h` RELAXED/STRICT base patches and explicit legacy velocity-seed compatibility | Cleaned branch 5 plus pinned formal documents and manuscript | `cav-2-patch-construction-audit` | Existing sparse/map/ASM/strict-equivalence tests plus exact ordered pressure-seed mapped patch test and legacy/non-equivalence discriminator | Patch definitions/order | Disposable integration: serial RELAXED constructor consumes live full-space `d_SAJ_mat` directly; one setup scan builds global-DOF row-or-column adjacency without velocity-block/transpose copies; zero/row/column tests, explicit production selection, and live two-level FAC path are green. STRICT extension, exact mapped oracle parity, and durable recovery remain open |
| BLAS/LAPACK LU shell backend | Historical 7 `4de0cd9fd` and LU/workspace portions of `4f83a0280` | `cav-3a-blas-lapack-lu-backend` | Live extracted matrix LU solve/backward error, lifecycle, failure propagation, shell semantics | Local matrices/solutions | Disposable minimal LU recovery is implemented with canonical global numbering outside the local factor boundary, one overwritten factor buffer, pivots, and reusable workspaces. Selected-only observer callbacks reconstruct transient unfactored local matrices on demand and retain no matrix copy. The positive observer/reinitialize fixture and normalized singular-factorization expected-failure fixture pass **2/2**. At controlled `K=1e4`, first coupled patch 236 has an exact mapped matrix, `7.11e-15` RHS difference, `2.30e6` 2-norm condition number, `3.14e-11` LU forward difference, and `2.31e-18` IBAMR backward error. PETSc SVD gives a comparable `3.25e-11` forward difference, so this is not an LU-only semantic defect. Durable review unit remains unimplemented |
| Robust BLAS/LAPACK modes | Historical 7 `4f83a0280` | `cav-3b-blas-lapack-robust-modes` | SVD/QR/symmetric-indefinite; near-singular rank/rcond; explicit Cholesky policy | Local matrices/solutions | Proposed after `cav-3a` |
| Global multiplicative composition | Cleaned 6.2–6.4 shared Eigen sweep plus historical 7 active-row intent and sandbox | `cav-4-multiplicative-composition` | Exact live per-patch residual/correction/state and fresh-original-residual tests | Corrections/residual updates | Proposed after M3 |
| One-sweep sandbox parity | Historical 8 live parity workflow, recovered by intent | `cav-5-smoother-parity` | Mandatory N0/N1 live export comparison; forward/reverse and RELAXED/STRICT | Smoother | Proposed after M4 |
| FAC/transfer/coarse observability | Historical 8 `ad3dcdb13`, `6f9037942`; cleaned 1.5 and 3.1–4 | `cav-6a-fac-stage-observability` | Focused stage snapshot, transfer/gauge, init/deinit/reinit, and regrid invalidation tests | Transfers; FAC stages | Disposable integration: passive const observer on the single production recursion exports fine pre/post and coarse RHS/correction stages through temporary level vectors only; no matrix copy or alternate recursion. Focused native lifecycle/regrid tests and durable recovery remain open |
| Multiplicative production V-cycle/FGMRES path | Cleaned FAC/transfers plus historical 8 parity/residual-sweep intent | `cav-6b-fac-vcycle-fgmres-parity` | N2–N4 two-source common-arithmetic replay; separate live V-cycle, FGMRES, momentum/divergence/IB residuals; integration-tip equivalence | V-cycle; FGMRES | Disposable reconstruction green at `K=1`, `1e2`, and `1e4`: candidate and oracle outputs are exactly equal in one MATLAB R2025b arithmetic runtime for one multiplicative sweep, two-level V-cycle correction/fresh residual, FGMRES solution/fresh residual, and relative iteration history on both independently exported live and oracle operators while consuming each side's own ordered patches. Separate live IBAMR runs converge in 9/8/9 iterations with original relative residuals `3.056e-11`, `9.178e-11`, and `1.245e-11`; live matrix-free-versus-`SAJ` IB-action relative defects are `2.976e-15`, `3.953e-15`, and `4.081e-15`. Raw common-replay-versus-exported-oracle final-residual differences remain named non-gating diagnostics (`1.78e-14`, `1.34e-12`, `1.13e-10`). Exact M5 durable tip and independent audits remain open |
| Macro-subdomain geometry/grouping/ownership | Sandbox/formal pin `5b77344d`; no IBAMR historical equivalent | `cav-7-macro-construction` | Deterministic macro order/membership, solve/owned partitions, overlap/restriction coverage | Patch definitions/order | Blocked by M5 all-lane gate |
| Macro-local multiplicative CAV | Sandbox/formal pin `5b77344d` | `cav-8a-macro-local-multiplicative` | Exact inner patch residual/correction/state sequence and completed macro correction | Corrections/residual updates | Blocked by M7 |
| Outer restricted-additive macro composition | Sandbox/formal pin `5b77344d` | `cav-8b-outer-restricted-additive` | Ownership uniqueness, completed-correction restriction, common-outer-state additive assembly | Corrections/residual updates; smoother | Proposed after `cav-8a` |
| Serial CAV-RAS FAC/V-cycle production integration | Sandbox oracle plus durable multiplicative production path | `cav-9a-cav-ras-fac-production` | R0/R1 plus native FAC/V-cycle/physical-residual tests | Smoother; V-cycle; FGMRES | Proposed after M8 |
| Serial CAV-RAS parity/acceptance and user surface | Pinned oracle/formal sources plus `cav-9a` | `cav-9b-cav-ras-parity-acceptance` | R2, full low-K matrix, high-K diagnostics, example/docs checks | Full ladder | Exact M9 serial production tip; proposed |

## 12. Six-lane independent audit architecture

The audit suite treats six concerns as distinct failure modes with independent ownership. Passing one lane never supplies evidence for another.

| Lane | Failure mode | Controlling prompt |
|---|---|---|
| A | Algorithm semantics | `plans/implicit-ib/cav-audit-algorithm-prompt.md` |
| N | Numerical reproduction/parity | `plans/implicit-ib/cav-audit-numerical-prompt.md` |
| D | Data-structure efficiency | `plans/implicit-ib/cav-audit-data-structures-prompt.md` |
| C | Cache/data-layout efficiency | `plans/implicit-ib/cav-audit-cache-layout-prompt.md` |
| E | Concrete extensibility | `plans/implicit-ib/cav-audit-extensibility-prompt.md` |
| S | Simplicity/minimality | `plans/implicit-ib/cav-audit-simplicity-prompt.md` |

`plans/implicit-ib/cav-independent-audit-prompt.md` version `cav-audit-v3` is the common committed protocol and launch index. Each auditor receives it plus exactly one versioned lane prompt. Reviewability is not a seventh technical lane; it is adjudicated after reports are frozen.

### 12.1 Permanent separation of duties, isolation, and custody

The implementation task must not perform, approve, coordinate, synthesize, or certify its own independent audit.

| Implementation task may | Implementation task must not |
|---|---|
| Maintain committed audit instructions and canonical ledger. | Execute any lane, synthesis, or audit-coordinator role. |
| Assemble provenance-bearing reproducibility artifacts from live IBAMR objects. | Issue, approve, downgrade, resolve, close, or certify findings. |
| Record implementation validation explicitly as implementation evidence. | Issue a lane/milestone outcome or mark an independent gate complete. |
| Address frozen findings in a new candidate commit. | Edit frozen reports or treat a prior outcome as applying to the new SHA. |

Auditors and the synthesizer are permanently read-only with respect to the candidate, oracle, comparator, pinned refs, prompts, and frozen reports. They may build and run tests and write only their own report/evidence in a separate audit-results workspace. They cannot receive later authorization inside the same task to edit or fix the candidate. Fixes always return to implementation and require newly launched affected-lane audits at the new SHA.

Every lane runs in its own fresh clean IBAMR worktree at the same exact candidate SHA and uses a second separate clean sandbox-oracle worktree detached at `5b77344d...`. The existing dirty/mismatched sandbox checkout is never audit evidence. Auditors receive committed instructions, pinned refs, candidate diffs, and reproducibility artifacts—not implementation conversational reasoning or another lane's report.

A separate audit coordinator, neither implementation nor a lane auditor in that suite, owns only launch coverage and immutable custody. It assigns a suite ID, records Git blob IDs and SHA-256 prompt versions, withholds reports until every lane submits, and freezes them under `audit-results/implicit-ib/cav/<candidate-sha>/<suite-id>/` in a separate append-only checkout with `MANIFEST.sha256`. It cannot edit findings or implementation.

### 12.2 Common exact-SHA, authority, provenance, and reporting protocol

Every lane must:

1. resolve and record candidate/parent/foundation/stack SHAs, local and remote refs, submodules, worktree paths, and complete dirty states;
2. use a separate clean detached oracle worktree at `5b77344d...`, or return `INCOMPLETE`;
3. read applicable `AGENTS.md`, this plan, the committed canonical ledger, recovery manifest, common protocol, and its lane prompt;
4. use the pinned formal documents and manuscript context in Section 2 while treating executable sandbox agreement as evidence rather than the definition of correctness;
5. rerun evidence from builds tied to active clean checkouts; old logs and implementation claims are context only;
6. use live IBAMR objects and retain raw pre-mapping exports with full provenance;
7. leave candidate, oracle, comparator, refs, and prompts unchanged;
8. report each finding with classification, severity, blocking status, exact SHA/file/line or immutable artifact/hash, reproduction, expected/actual, impact, and recommended disposition;
9. conclude with exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.

The numerical lane independently validates the exporter/mapping/comparator with raw exports, a live hand-checkable `4x4` case, exact-match controls, and activated mutations for permutations, signs, legal/illegal gauges, omissions, and reorderings. Failure to validate the comparison machinery makes that lane `INCOMPLETE`.

### 12.3 Lane contracts

#### A — Algorithm semantics

Independently derive global multiplicative CAV and macro-subdomain CAV-RAS from the formal authorities, then trace production code and tests. Reconcile the pressure-cell-seeded row-or-column `E_h` formulation with IBAMR's velocity-seeded row-`A00` RELAXED/STRICT construction without assuming equivalence. Require original-residual evolution, pressure gauge, ordered traversal, inner multiplicative macro sweeps, completed-correction restriction, unique ownership, and outer additive assembly. Flat inner additive-restrict fails if labeled CAV-RAS.

#### N — Numerical reproduction/parity

Use the mandatory matrix and formulas in Section 9 and the lane prompt. Localize the first mismatch through configuration; live operator/blocks; mapping/sign/gauge; seed reconciliation and patch order; local solves; per-patch states; macro restriction/assembly; transfers; smoother; V-cycle; FGMRES; and original momentum/divergence/IB residuals. Own numerical production robustness in supported dimension/MPI scope and never use conditioning to excuse `K <= 1e4`.

#### D — Data-structure efficiency

Separate setup from apply work. Audit scaling, ownership, allocation lifetime, lookup, conversion, recomputation, synchronization, serial/distributed responsibilities, initialize/deallocate/reinitialize, hierarchy regridding, cache invalidation, stale-source rejection, and retained resources. Findings require measured or bounded impact, not container taste.

#### C — Cache/data-layout efficiency

Focus on the hot path using the fixed Release workloads in the common protocol. Audit locality, indirection, temporaries, traversal, cached factors/mappings, memory footprint, and invalidation. Use paired hardware/build provenance, warmups/repetitions, median/MAD normalization, baselines, and the >10% plus noise-floor materiality rule. Do not demand speculative optimization.

#### E — Concrete extensibility

Consider only RELAXED/STRICT; multiplicative/CAV-RAS; local-solver backends; forward/reverse/palindromic traversal; and eventual serial/distributed execution. Require narrow tested seams and structured next-change exercises without hypothetical frameworks.

#### S — Simplicity/minimality

Be deletion-oriented and audit every production unit that introduces complexity. Challenge abstractions, state, duplicate paths, generic mechanisms, options, and compatibility layers. Ask whether roughly 30% could be deleted while preserving requirements, tests, measured performance, and concrete axes. Retained complexity needs specific evidence.

### 12.4 Audit-lane-to-review-unit matrix

Do not run all six lanes on every review unit. `R` means required; blank means not routinely assigned unless a finding or changed scope expands coverage. Simplicity is required where production complexity is introduced.

For every production unit, lane N is the named owner for configuring/building the active checkout, verifying `attest` discovery, running the exact fast native cases through repository-root `attest`, and recording command, discovered names, fixtures, expected status/exit semantics, budget, measured runtime, deterministic rerun, and meaningful failure sensitivity. Other assigned lanes verify the accuracy of comments for their own invariants; synthesis owns the cross-unit same-commit/PR-boundary, `attest` discoverability, executable/fixture reuse, and CI-suitability check, not test execution.

| Review unit / gate | A | N | D | C | E | S |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| `cav-0-recovery-plan` |  |  |  |  | R | R |
| `cav-1a-live-operator-audit-schema` |  | R | R |  |  | R |
| `cav-1b-raw-export-comparator-contract` |  | R | R |  | R | R |
| `cav-1c-matrix-live-kernel-consistency` | R | R | R |  |  | R |
| `cav-2-patch-construction-audit` | R | R | R |  | R | R |
| `cav-3a-blas-lapack-lu-backend` |  | R | R | R | R | R |
| `cav-3b-blas-lapack-robust-modes` |  | R | R | R | R | R |
| `cav-4-multiplicative-composition` | R | R |  | R | R | R |
| `cav-5-smoother-parity` | R | R |  |  |  |  |
| `cav-6a-fac-stage-observability` |  | R | R | R |  | R |
| `cav-6b-fac-vcycle-fgmres-parity` | R | R | R | R |  | R |
| `cav-7-macro-construction` | R | R | R |  | R | R |
| `cav-8a-macro-local-multiplicative` | R | R | R | R | R | R |
| `cav-8b-outer-restricted-additive` | R | R | R | R | R | R |
| `cav-9a-cav-ras-fac-production` | R | R | R | R | R | R |
| `cav-9b-cav-ras-parity-acceptance` | R | R | R | R | R | R |
| **M5 working multiplicative gate at exact durable cav-6b tip** | **R** | **R** | **R** | **R** | **R** | **R** |
| **M9 working serial-CAV-RAS gate at exact cav-9b tip** | **R** | **R** | **R** | **R** | **R** | **R** |

The M5 all-lane gate must pass before M6 retires integration scaffolding and before M7 starts. The M9 all-lane gate certifies serial rank-1 CAV-RAS only. Foundation restack, each lettered review tip, review-stack equivalence, RAS construction, and future distributed RAS remain separate entry points.

### 12.5 Production-robustness ownership

Production robustness is explicitly assigned within the six lanes:

| Concern | Primary | Supporting |
|---|---|---|
| Debug/Release 2D and 3D build/smoke | N | D |
| supported MPI ranks and rank consistency | N | A, D |
| initialization/deallocation/reinitialization | D | N, C |
| hierarchy regrid/cache invalidation/stale-source rejection | D | N, C |
| local-solver singular/ill-conditioned/failure behavior | N | A, D |
| NaN/Inf and sanitizer/resource/leak checks | N | D, C |
| tiny/empty/boundary/periodic/no-IB/ownership-boundary cases | N | A, D |

The audit coordinator verifies that every row has evidence and a named report owner. Unsupported scope must reject cleanly or be explicitly excluded. It cannot silently pass.

### 12.6 Executable performance contract

C and D run the common protocol's fixed workloads P1 (single-level `64x64` multiplicative sweep), P2 (three-level `16/32/64` V-cycle), P3 at M9 (same fine workload with `8x8` RAS macro blocks, standard ghost 0, IB ghost 2), and ten setup/deallocate cycles. Candidate and exact parent use identical Release toolchains/hardware/binding. Record warmups, repetitions, medians, MADs, normalized time, allocations, bytes, peak RSS, and scaling where claimed.

A material regression is a normalized time or peak-memory increase above 10% whose delta also exceeds three times the larger MAD/timer resolution, a new unnecessary per-apply allocation/recomputation scaling with problem size, or a reproducibly superlinear path contradicting a linear claim. Missing comparable evidence yields `INCOMPLETE`; intuition alone cannot fail the lane.

### 12.7 Audit-suite sensitivity

Fault exercises mutate only disposable harness state/raw copies or use purpose-built fixtures; candidate and oracle remain untouched. Every exercise records an ID/hash, activation proof, unmutated control, expected detector, and actual detector.

| Lane | Required sensitivity |
|---|---|
| A | stale residual, patch order, RAS ownership/restriction, overlap, gauge, inner additive substitution |
| N | A defects plus mapping permutation/sign/gauge/omission/reordering and first-mismatch localization |
| D | forced per-apply rebuild, duplicated full-domain ownership metadata, or equivalent lifecycle/recomputation challenge |
| C | cache-disable/per-patch allocation or locality-degrading traversal with measured material effect |
| E | fixed concrete-axis change exercises plus a seeded cross-layer-entanglement example |
| S | per-unit/30% deletion challenge plus a seeded redundant abstraction/state/option example |

An inactive mutation makes the exercise `INCOMPLETE`. Failure to detect an activated required defect is a blocking audit-suite-readiness failure, not evidence that production necessarily contains that defect.

### 12.8 Severity, outcomes, frozen reports, synthesis, and reviewability

Findings use independent classification (`confirmed defect`, `untested risk`, `optional cleanup`) and severity:

| Severity | Meaning and blocking rule |
|---|---|
| `S0 critical` | wrong algorithm, corrupted state/result, or unsafe resource behavior; always blocking |
| `S1 high` | required correctness, reproducibility, supported robustness, or material performance failure; always blocking |
| `S2 medium` | material scoped weakness/evidence gap; blocking when it touches a milestone requirement |
| `S3 low` | local cleanup/clarity issue without acceptance impact; normally nonblocking |

Lane outcomes are `PASS` (complete evidence, no blocker), `FAIL` (confirmed blocker or mandatory failure), and `INCOMPLETE` (missing/invalid environment, evidence, scope, or sensitivity control). `INCOMPLETE` blocks a gate.

After coordinator custody freezes all reports for one exact SHA, a separate permanently read-only synthesis task verifies hashes and coverage. It preserves every finding, assesses reviewability, reconciles cache duplication versus simplicity and concrete extensibility versus overengineering, retains dissent, and may recommend PR-boundary changes. It does not edit code, findings, or reports and cannot close anything.

Reviewability means one coherent claim, focused same-unit tests, normal `attest` discoverability/execution with required fixtures and expected status, recorded CI budget/runtime and sensitivity, no hidden later dependency, refactor separate from behavior, construction separate from application, comments conforming to surrounding IBAMR practice, and no canonical-ledger regression. Missing, undiscoverable, nondeterministic, nonsensitive, or unjustifiably slow required native tests are blocking or make the responsible lane `INCOMPLETE`. Materially misleading comments are defects; optional comment-density preferences are not.

### 12.9 Re-audit rules

Every fix produces a new candidate SHA and new/versioned suite. Prior outcomes do not transfer.

- Rerun every lane affected by changed files, behavior, evidence contract, or prior finding.
- A and N are mandatory after semantic or cross-cutting changes.
- D and C rerun after ownership, storage, lifecycle, traversal, resource, or hot-path changes.
- E and S rerun after API, option, abstraction, state, backend, or branch-boundary changes.
- Newly launched auditors rerun evidence independently; the implementation task cannot resolve/close its own findings.
- Freeze new reports and rerun synthesis before the new exact-SHA gate can pass.

## 13. Decision log

| Date | Decision | Reason/evidence |
|---|---|---|
| 2026-08-18 | Use `0c8447e3593113a8b62f94df9482c60d896f0d3b` as the only recovery base. | Verified branch 7 diverges before the cleaned stack (`merge-base bd010f326...`), so historical continuation assumptions are unsafe. |
| 2026-08-18 | Maintain disposable integration reconstruction and durable review stack as separate green tips. | End-to-end recovery and small-PR reviewability are different jobs; neither is allowed to substitute for the other. |
| 2026-08-18 | Recover historical 7/8 by intent, not wholesale cherry-pick. | Their commit topology includes development-era boundaries and a divergent pre-cleanup foundation. |
| 2026-08-18 | Gate all CAV-RAS work on durable current-master multiplicative parity. | Prevents conflating a base smoother recovery error with macro grouping/composition behavior. |
| 2026-08-18 | Treat patch order, gauge, and original-residual recomputation as algorithm semantics. | Multiplicative corrections depend on traversal and current residual; saddle-point comparisons can otherwise appear to disagree or agree spuriously. |
| 2026-08-18 | Treat existing flat additive-then-restrict patch processing as non-equivalent to sandbox CAV-RAS. | It lacks multiplicative residual evolution inside each macro-subdomain. |
| 2026-08-18 | Establish permanently read-only independent auditor and synthesizer tasks, with fixes reserved to implementation and fresh audits at every fixed candidate SHA. | Reviewability and integrated numerical correctness require fresh evidence from exact clean tips; implementation authority must not leak into audit, synthesis, or finding closure. |
| 2026-08-18 | Pin the sandbox algorithmic oracle at `5b77344db6746269f8c77695c99e9043907ba74b` for implementation and independent audit. | The commit is verified and introduces the required global multiplicative and macro-subdomain CAV-RAS processing; later sandbox commits or dirty working-tree changes must not drift the oracle silently. |
| 2026-08-18 | Assign independent-audit findings, dispositions, and milestone gates exclusively to a separate auditor task pinned to an exact candidate SHA. | The implementation task may prepare evidence and address findings but cannot audit, approve, close, or certify its own work; every new candidate requires a fresh independent re-audit. |
| 2026-08-18 | Replace a generic technical audit with six independent lanes: algorithm, numerical, data structures, cache/layout, concrete extensibility, and simplicity. | These are distinct failure modes; scoped lane assignment avoids ritual all-lane review on every unit while both working-algorithm gates require all six. |
| 2026-08-18 | Place reviewability and cross-lane tension resolution in a post-freeze synthesis task, not a seventh technical lane. | Synthesis can reconcile cache duplication versus simplicity and concrete extensibility versus overengineering without hiding technical findings or editing code. |
| 2026-08-18 | Require controlled audit-suite fault injection and exact-SHA re-audit of every affected lane. | An audit suite must demonstrate sensitivity to representative semantic defects; prior passes cannot certify a changed candidate. |
| 2026-08-18 | Canonicalize the historical ledger in-repository and require a separate clean oracle worktree detached at `5b77344d...`. | Audit inputs must not depend on an ephemeral attachment or the currently dirty/mismatched sandbox checkout. |
| 2026-08-18 | Pin the formal algorithm documents and manuscript context, while leaving pressure-cell-seed versus velocity-seed equivalence as an explicit M2 proof obligation. | Executable agreement is evidence rather than a mathematical definition, and the two constructions are not equivalent by assumption. |
| 2026-08-18 | Fix the low-stiffness matrix at `K=1,1e2,1e4`, with `2e-11` as the `K=1e4` `E_inf` target and `1e-10` as the hard final-output ceiling. | This retains a fixed `1e-11`-scale criterion while accommodating the historical STRICT discrepancy of approximately `1.289e-11`; conditioning-aware criteria remain forbidden at low K. |
| 2026-08-18 | Assign immutable suite custody to a separate audit coordinator and production robustness across N/D/C, with lane outcomes `PASS`, `FAIL`, or `INCOMPLETE`. | Exact prompt/report hashes, fixed performance workloads, lifecycle/build/MPI evidence, and explicit evidence gaps are required for reproducible gates. |
| 2026-08-18 | Replace broad `cav-1`, `cav-3`, `cav-6`, `cav-8`, and `cav-9` units with the lettered M0-confirmed decomposition. | Historical 7/8 inventory shows independently reviewable seams between schema/export, core/robust backends, observability/behavior, inner/outer RAS composition, and production/acceptance. |
| 2026-08-18 | Mark M0 implementation inventory/validation complete at the exact foundation after the focused Debug selection passes 69/69. | The exact 1.1–6.4 ancestry, divergent historical/aggregate refs, historical 7/8 capabilities/tests, build commands, and durable review boundaries are now recorded; no independent audit has run. |
| 2026-08-18 | Treat pressure-seeded `E_h` construction as the formal RELAXED parity target and retain the foundation's full-`A00` velocity-seeded RELAXED/STRICT path only as explicitly named legacy compatibility behavior. | Source inspection at the pinned oracle shows pressure seeds, standard Vanka `U_i`, and row-or-column Eulerian-elasticity adjacency; foundation production instead expands a selected velocity seed through the full Stokes-plus-IB `A00`, so equivalence by index mapping is false in general. |
| 2026-08-19 | Stage the first pressure-seeded constructor as serial, RELAXED-only disposable integration code, with distributed ownership and the declared STRICT extension left as separate explicit obligations. | The formal/sandbox RELAXED rule can be tested independently through zero-`E_h`, row-edge, and column-edge cases; silently claiming partial distributed or STRICT semantics would weaken both correctness and reviewability. |
| 2026-08-19 | Give the level solver an explicit `PRESSURE_CELL` versus `VELOCITY_COMPONENT` seed-family choice and a separate non-owning construction-matrix channel. | FAC can now supply its live `d_SAJ_mat` only to patch construction while the shell backend continues to extract local matrices and update residuals with the full original saddle operator; the legacy family remains the default until production inputs opt in. |
| 2026-08-19 | Extract the velocity-space construction block from FAC's live `d_SAJ_mat` at the level-solver seam instead of assuming `d_SAJ_mat` is already velocity-sized. | The first two-level live run measured 192 rows in `d_SAJ_mat` versus 128 velocity DOFs; the constructor's matrix-space guard caught the full-coupled representation. Patch adjacency must use only its live velocity submatrix, while local solves continue to use the separate full saddle operator. |
| 2026-08-19 | Supersede the proposed `d_SAJ_mat` velocity-submatrix extraction: use the existing full global velocity-pressure index space for patch construction and delay local numbering until factorization/application. | The full-space FAC representation already shares indices with pressure seeds, Vanka closures, final local saddle extraction, and residuals. Treating only its velocity-to-velocity nonzeros as the `E_h` graph removes a redundant mapping/extraction seam and is simpler to validate; nonzero pressure adjacency is an error. |
| 2026-08-19 | Make the full global construction representation exclusive for the recovered pressure-seeded path and build row-or-column adjacency with one direct setup scan. | This removes the temporary velocity-submatrix allocation, its index translation, and the materialized transpose. The compact adjacency cache stores only the graph needed for repeated patch queries, in the same global DOF numbering; focused construction, shell, and live two-level FAC tests pass 5/5. |
| 2026-08-19 | Make one full velocity-pressure global index space and borrowed live matrices the production representation contract. | Seed, patch, residual, correction, export, and parity objects can share one numbering; only local factorization needs compact patch numbering. Avoiding duplicate `A`/`E_h`/velocity-block/transpose/backend matrices removes translations, ownership hazards, setup cost, and memory growth. Any unavoidable conversion now requires an explicit consumer, lifetime, test, and measured cost. |
| 2026-08-19 | Add the first disposable `ex0` construction bundle as a streaming view of live full-space objects, with generated output confined to `/tmp`. | The bundle records configuration, global DOF maps, pressure-seed/patch order, and MatrixMarket views of each live level `A` and `E_h` without extracting or converting matrices. This supplies the first-mismatch problem/construction inputs while keeping the exporter out of production algebra and commits. |
| 2026-08-19 | Require explicit 17-digit row streaming for raw live MatrixMarket exports and declare the pressure-equation row-sign map separately from the DOF permutation. | The first pinned-oracle comparison found exact mapped sparsity and patch structure but PETSc's default ASCII MatrixMarket view rounded live entries to roughly six digits, creating a false `1e-6`-scale `E_h` discrepancy. It also localized the genuine convention difference to divergence rows (`IBAMR=-div`, sandbox `=+div`). Export precision and equation-sign mapping must be visible comparator inputs, never silent normalizations. |
| 2026-08-19 | Accept the first implementation-task `K=1e2` problem/operator/base-patch parity layer after full-precision export, while leaving independent audit and the complete numerical matrix open. | The declared logical-DOF bijection has zero coordinate error; metadata and stiffness mapping are exact; both level `A`/`E_h` numerical sparsity patterns are exact; all 1024 fine pressure seeds, overlap memberships, and first-owner subsets agree; and mapped `E_inf` is `2.05e-15` coarse and `4.92e-15` fine after the explicit pressure-row sign map. This is implementation evidence only. |
| 2026-08-19 | Accept the first implementation-task `K=1e2` complete finest-level multiplicative sweep comparison through the active PETSc shell backend. | A deterministic live `PCApply` probe uses the actual 1024 patch matrices/local solves, current-residual recomputation, correction assembly, and pressure postprocessing. Against the pinned sandbox, mapped `E_inf` is `4.25e-14` for the correction and `7.97e-15` for a fresh original residual. No comparison matrix or alternate solver path was introduced; direct per-patch trace evidence and independent audit remain open. |
| 2026-08-19 | Recover branch-8 FAC diagnostics by intent as a passive observer on the single production recursion, not as a copied diagnostic recursion with retained hierarchy vectors. | Historical branch 8 duplicated V-/mu-cycle control flow and copied multiple full hierarchy vectors. A const stage observer can stream requested pre/post-smoother and coarse-solve states without changing the cycle, owning algebra, or copying matrices; this preserves one semantic path and fits `cav-1a-native-observability-schema`. |
| 2026-08-19 | Permit FAC parity exports to allocate only one temporary PETSc vector for the requested level while continuing to forbid copied matrices and retained diagnostic hierarchy state. | The stage observer receives const live hierarchy vectors and streams each selected level immediately. This is an explicit, bounded serialization boundary; it neither creates a second algebraic representation in production nor changes the observed cycle. |
| 2026-08-19 | Make RT0 velocity refinement/coarsening and linear/conservative pressure transfer choices explicit in the multiplicative parity configuration. | The first-mismatch trace showed matching fine pre-smoothing but a `2.23e-1` coarse-RHS error under the generic constant/conservative velocity defaults. Selecting `RT0_REFINE`/`RT0_COARSEN` makes the coarse RHS agree to `5.49e-15`, the coarse correction to `3.40e-14`, and the prolongated fine state to `3.06e-14`; silent transfer drift is therefore an algorithmic mismatch. |
| 2026-08-19 | Use one-step Richardson, not `preonly`, to invoke the multiplicative shell as a FAC level smoother when the incoming error is nonzero. | `preonly` overwrote the prolongated post-smoothing state and produced a `3.77e-1` final-correction discrepancy. One Richardson step applies the shell to the current residual and adds the correction, matching the sandbox post-smoother and final V-cycle map without adding a copied operator or alternate apply path. |
| 2026-08-19 | Accept the first implementation-task `K=1e2` two-level V-cycle comparison after stagewise first-mismatch localization. | Fine RHS and pre-sweep agree at `1.25e-14` and `1.40e-14`; coarse RHS/correction and prolongated post input agree at `5.49e-15`, `3.40e-14`, and `3.06e-14`; the final correction agrees at `2.39e-14`, and its fresh residual has `7.94e-13` normalized error (`2.44e-12` absolute). The external trace is asserted against unmodified pinned `v_cycle_precomp`; this remains implementation evidence, not an independent audit. |
| 2026-08-19 | Accept the implementation-task mandatory low-stiffness two-level V-cycle matrix at `K=1`, `1e2`, and `1e4`. | At `K=1`, final correction error is `1.63e-15` and fresh-residual error is `2.73e-14`; at `K=1e2`, they are `2.39e-14` and `7.94e-13`; at `K=1e4`, final correction error is `9.43e-13` normalized and `4.72e-12` absolute, while the fresh-residual error is `6.51e-12` normalized. All final corrections satisfy the fixed low-`K` target and `1e-10` hard output limit without conditioning allowances. |
| 2026-08-19 | Keep the `K=1e4` entrywise fresh FGMRES residual as an unresolved low-stiffness parity failure even though the final solution, iteration history, and physical residual norms agree. | At `K=1` and `1e2`, mapped final-solution errors are `3.42e-16` and `1.34e-14`, and mapped fresh-residual errors are `1.95e-14` and `3.27e-12`. At `K=1e4`, the final solution is green at `1.39e-12`, all 9 residual-history entries agree to `1.72e-15`, and physical total/momentum/divergence norms differ by at most `8.04e-13`, but the mapped fresh residual has absolute/`E_inf` error `2.51e-10`, above the fixed `2e-11` target and `1e-10` ceiling. Do not waive this as conditioning; localize local-solve/backend arithmetic first. |
| 2026-08-19 | Recover a minimal local LU backend by intent before adding historical robust factorization modes or cached residual-update matrices. | Passive traces of the actual live backend show exact mapped patch membership and sparsity. Representative IB-coupled local matrices agree with the oracle to at most `2.64e-15` normalized, but solving each exact same IBAMR local matrix/RHS with MATLAB LU changes the PETSc-SVD result by `2.08e-13` on the first 118-DOF coupled patch and `1.47e-13` on the 221-DOF maximum patch; LU backward errors are approximately `2.05e-18` and `9.00e-19`. This makes LU recovery a correctness-motivated experiment. The implementation must borrow the full global operator, store only necessary patch-local factors/pivots and reusable vectors, and must not recover branch 7 wholesale or copy full/update matrices without later measured justification. |
| 2026-08-19 | Keep the IBAMR-to-sandbox sign reconciliation as a pressure-equation row map only: IBAMR uses `-Div`, while the sandbox uses `+Div`. | Source inspection shows IBAMR assembles `[H,G;-D,0]` and the sandbox `[H,G;D,0]`; their gradient/pressure-column convention agrees. Before mapping, the finest operator has a `64` maximum discrepancy (`-32` versus `+32`) and the smoother RHS has a `97.02` discrepancy. Multiplying only mapped IBAMR pressure-equation rows/RHS/residual entries by `-1` reduces the operator and RHS comparisons to `4.63e-15` and `9.28e-15`. Corrections and pressure solution values receive no sign change. The remaining `K=1e4` FGMRES maximum residual discrepancy is at velocity DOF 489, where this map is the identity, so it is not a Div-sign artifact. |
| 2026-08-19 | Do not retain copied residual-update matrices in the minimal LU backend based on the present evidence. | The disposable backend borrows the full global operator, overwrites one necessary dense local buffer with LU factors, and retains factors/pivots plus reusable global residual/correction/action vectors and local RHS workspaces. An incremental live-matrix update gives green low-`K` smoother/V-cycle results, while a same-live-operator MATLAB replay still differs from the C++ sweep by about `5.5e-13`; therefore copying historical dense update blocks is not established as either necessary or sufficient. Any later cache must have measured setup/apply/memory evidence and a demonstrated numerical or performance consumer. |
| 2026-08-19 | Keep the `K=1e4` FGMRES residual failure assigned to both final solver sensitivity and the live matrix-free/assembled operator seam, not to CAV construction or the pressure sign map. | With minimal LU, `K=1` and `1e2` pass through FGMRES; `K=1e4` passes patch structure, smoother, V-cycle correction (`8.32e-13`), final solution (`1.49e-12`), and 9-iteration agreement, but its mandatory matrix-free residual error is `2.69e-10`. Re-evaluating the same solution with the FAC-assembled Stokes-plus-SAJ operator reduces the error only to `1.93e-10`; the independently measured matrix-free/SAJ IB-action defect is `7.95e-10`. PETSc modified Gram-Schmidt worsened the entrywise result, and always-refined classical Gram-Schmidt was bit-for-bit unchanged, so both experiments were removed. The existing matrix-based SAJ Jacobian path rejects two-level vectors and cannot be used to bypass the production residual contract. |
| 2026-08-19 | Treat the `K=1e4` outer-FGMRES trace as a cancellation-sensitive numerical seam, not evidence of a new CAV semantic mismatch. | A passive PETSc iterate observer and an external instrumented copy that exactly reproduces pinned `fgmres_opt` show an initial mapped residual difference of `1.27e-14`, the first post-V-cycle iterate difference at iteration 1 (`4.39e-12`), and the first iterate/residual crossing of `1e-11` at iteration 2. At iteration 3 the residual difference grows to `2.91e-10` while the iterate difference falls to `2.16e-12`; the trace therefore requires a same-operator/same-RHS decomposition before changing the solver, tolerances, or matrix ownership policy. |
| 2026-08-19 | Move the active `K=1e4` first mismatch to live RHS formation, ahead of outer FGMRES. | The ten-state residual identity closes to at most `2.22e-16`. The mapped live RHS values already differ by `3.13e-10` at the identical zero iterate. At later iterates, effects from the iterate under the common oracle matrix (`2.8e-10`--`8.5e-10`), assembled operator/RHS representation (`2.0e-10`--`2.4e-10`), and IBAMR matrix-free versus assembled application (`1.5e-10`--`3.1e-10`) cancel to produce the observed final residual difference. Do not tune FGMRES or relax low-`K` criteria until the identical live seed separates matrix-free RHS formation from assembled-operator representation. |
| 2026-08-19 | Keep both assembled representation and matrix-free application seams active after exact live-seed reconciliation. | The mapped `K=1e4` RHS seeds, pressure means, and oracle recomputation agree exactly. Applying the exact seed separates the `3.13e-10` live RHS difference into a `3.71e-10` mapped assembled-operator effect and a `2.29e-10` IBAMR matrix-free-versus-assembled effect that partially cancel. Neither changing FGMRES nor copying a matrix addresses this first mismatch; split both effects into Stokes and IB-elasticity actions next. |
| 2026-08-19 | Assign the dominant exact-seed arithmetic seam to IB-elasticity composition, not the ordinary Stokes operator. | After forming the Stokes remainder at the matrix-entry level instead of subtracting already-applied large actions, its mapped representation difference is at most `1.42e-14` at `K=1`, `1e2`, and `1e4`. The assembled IB-elasticity effect grows from `2.71e-14` to `3.87e-12` to `3.71e-10`; the IBAMR matrix-free/assembled IB effect grows from `2.04e-14` to `1.99e-12` to `2.00e-10`. Explicit `A=(A-E_h)+E_h` recomposition contributes `5.82e-11` at `K=1e4`. This is strong stiffness-scaling evidence but does not waive the fixed low-`K` ceiling. |
| 2026-08-19 | Classify the remaining exact-seed IB seam as floating-point composition order, not a seed, coefficient, ownership, sign, or duplicated-force-Jacobian defect. | A disposable trace applies the same live objects four ways without retaining or copying another production matrix. The FAC and matrix-free Lagrangian force Jacobians are exactly equal. The live hierarchy interpolation and `MatMult(J,seed)` differ by only `2.60e-18` at every stiffness; after force application the differences are `1.39e-16`, `1.42e-14`, and `1.59e-12` for `K=1`, `1e2`, and `1e4`. FAC `MatPtAP` versus staged `J^T A J` contributes `1.64e-14`, `2.08e-12`, and `2.16e-10`; hierarchy spread versus `MatMultTranspose(J,force)` contributes `2.66e-15`, `2.27e-13`, and `2.18e-11`. The final direct hierarchy-versus-SAJ defects are `2.04e-14`, `1.99e-12`, and `2.00e-10`. Seed order/membership, scales, signs, matrices, and direct hierarchy recomposition are exact or roundoff-level; the fixed `K=1e4` failure remains open. |
| 2026-08-19 | Reject an extra left-associated production matrix product as the remedy for `K=1e4`; localize the factor-level difference to IB4 weights. | Replaying both `((J^T A)J)` and `J^T(AJ)` in MATLAB from raw live IBAMR factors does not improve the pinned-oracle matrix or exact-seed action relative to FAC `MatPtAP`, so an intermediate setup matrix has no demonstrated correctness benefit. The scaled 80-by-80 Lagrangian force Jacobian is exactly equal to the oracle factor. The mapped 80-by-3072 interpolation matrix has exact numerical sparsity, but corresponding IB4 weights differ by up to `8.05e-16`; stored-entry counts differ only because IBAMR retains explicit zeros. The previous “no coefficient defect” classification therefore means no semantic or scale defect, not bitwise-equal interpolation weights. Trace the formulas/coordinate evaluation before considering any production change. |
| 2026-08-19 | Make matrix/live formula identity a production invariant for every regularized-delta kernel supported by the matrix builder, as a separate `cav-1c` foundation unit. | The IB4 trace exposed two independently meaningful effects: 12 of 80 circle-coordinate values differ from MATLAB by one binary64 unit because C++ and MATLAB evaluate trigonometry separately, while the matrix and live paths duplicate coordinate/weight formulas. Cross-language geometry remains a mapping concern; within IBAMR, every matrix/live interpolation and spreading path must use the same coordinate, stencil, weight, and arithmetic contract in 2D/3D. Focused exhaustive tests precede production formula edits, and no matrix copy or alternate operator is permitted as a substitute. |
| 2026-08-19 | Correct matrix IB5/IB6 semantics and live 3D IB5 spreading inside the isolated `cav-1c` claim, and use one exhaustive row/action test rather than per-kernel surrogate matrices. | A live 2D baseline found IB5 and IB6 matrix-row errors of `1.2444e-1` and `9.5368e-2`; 3D found `6.7290e-2` and `5.7549e-2`, plus a `1.6410e-1` IB5 spreading error. IB5 had been evaluated as an incorrect pointwise scalar instead of the live centered five-weight recurrence, IB6 used the reversed fractional coordinate (`r-2` instead of the live-equivalent `3-r`), and the 3D IB5 spreading innermost loop omitted `ic0=ic_lower(0)+i0`. After aligning all matrix callback formulas with the live formulas and fixing that loop index, all 13 shared choices pass direct matrix-row versus live interpolation and volume-scaled spreading in 2D/3D. The test introduces no retained matrix, no alternate operator, and no persistent index translation. |
| 2026-08-19 | Retain the corrected `spread_01_3d.ib_5.output` fixture as native regression evidence only after an independent operator-action check identified the old fixture as encoding the missing `ic0` update. | The old expected output restricted x support to indices 4--6 because the live 3D IB5 spreading loop reused a stale x index. The corrected implementation has the five-point IB5 x support at indices 2--8 for this case, agrees with the independently constructed matrix transpose, passes the focused rerun, and leaves the complete serial interpolation/spreading suite green **126/126**. This is an investigated semantic correction, not an unexplained golden-output refresh. |
| 2026-08-19 | Keep the all-kernel correction even though the production IB4 CAV parity case is byte-for-byte unchanged, and retain the existing `K=1e4` failure classification. | Fresh IBAMR and pinned-oracle exports at `K=1`, `1e2`, and `1e4` reproduce the prior IB4 operator, smoother, V-cycle, FGMRES solution, and original-residual vectors exactly. The mapped smoother/V-cycle/original-residual errors remain `7.77e-16`/`1.64e-15`/`2.58e-14`, `3.86e-14`/`2.35e-14`/`2.39e-12`, and `6.80e-13`/`8.32e-13`/`2.69e-10`, respectively. Thus the fixes are necessary for IB5/IB6 consistency but neither improve nor regress IB4; `K=1e4` remains outside the fixed low-stiffness final-output limit and is not promoted. |
| 2026-08-19 | Treat missing BLAS/LAPACK-backend local-solve observer callbacks as an open evidence-plumbing defect, not as a failed numerical comparison or a reason to weaken parity. | Enabling the disposable `CAV_EXPORT_LOCAL_SOLVE_TRACE` path writes the selected-patch manifest but no local matrix/RHS/solution files because the observer is currently invoked only by the PETSc shell backend. Symmetric FGMRES evidence was regenerated without requesting that incomplete trace. Observer support belongs with the recovered local-solver backend review unit and must be implemented and tested before the local-solve evidence gate; current native LU semantics tests remain valid implementation evidence. |
| 2026-08-19 | Close the disposable all-kernel focused-test coverage only after actual periodic wrapping and sensitivity controls pass for every shared kernel. | Two additional points place every coordinate inside the lower and upper boundary cells of the one-patch periodic grid. Basis interpolation populates every periodic representative of one canonical global DOF; spreading folds ghost-region contributions back by the live DOF map. All 13 shared choices exercise a stencil containing both low- and high-side canonical indices and pass in 2D/3D. Captured matrix rows also prove that a coefficient perturbation and an unequal-weight column swap are rejected for every choice. Maximum observed differences remain below the unchanged `1024 epsilon_machine` bound. |
| 2026-08-19 | Make exact exported live marker coordinates the controlled-input variant for separating implementation parity from cross-language geometry generation, but do not treat it as a waiver of fixed low-`K` limits. | At `K=1e4`, an external wrapper feeds the 80 raw IBAMR coordinate values to the unmodified pinned oracle functions. The mapped `J` structure is exact with maximum entry difference `2.78e-17`, and the scaled force Jacobian is exact. The final residual error improves from `2.69e-10` to `2.40e-10`, while the smoother, V-cycle, and solution errors are `4.92e-13`, `1.00e-11`, and `1.12e-12`; the hard final-residual ceiling is still missed. Exact input therefore removes one problem-generation effect but leaves matrix-composition and live-versus-assembled action order as the first mismatch. The oracle checkout remains clean at its exact pin. |
| 2026-08-19 | Keep the `K=1e2` fixed residual miss assigned to implementation arithmetic after completing the exact-shared-coordinate matrix. | Feeding each case's raw live IBAMR marker coordinates to the external pinned-oracle wrapper gives smoother/V-cycle/solution/original-residual errors of `6.66e-16`/`1.59e-15`/`2.78e-16`/`2.49e-14` at `K=1`, `4.20e-14`/`2.25e-14`/`1.45e-14`/`2.47e-12` at `K=1e2`, and `4.92e-13`/`1.00e-11`/`1.12e-12`/`2.40e-10` at `K=1e4`, with exact iteration counts `9`, `8`, and `9`. The `K=1e2` original-residual error therefore remains above its fixed `1e-12` target even after removing independent geometry generation. Do not attribute it to input geometry or relax the contract; continue at the factor/composition first mismatch. |
| 2026-08-19 | Preserve the live IBAMR formulas for every shared kernel and do not emulate MATLAB's binary64 elementary-function rounding. | The maximum exact-coordinate IB4 factor mismatch is row `60`, mapped velocity column `522`, at marker `(0.49999999999999994, 0.3)` and side DOF `(axis 0, [16,10])`: IBAMR stores `0.14788689868556618` and MATLAB stores `0.14788689868556620`. Stencil indices and the fractional coordinates `r_x=0.9999999999999982`, `r_y=0.09999999999999964` agree; the discrepancy is one binary64 unit in the square-root evaluation/tensor product, not a formula, sign, order, or support defect. Replacing only this one coefficient in a fixed-order external `32 J^T(-F)J` replay changes an `E_h` entry by `9.09e-13` at `K=1e2` and `2.33e-10` at `K=1e4`; replacing all observed `J` coefficients changes an entry by `3.64e-12` and `4.66e-10`. The stiffness-scaled seam is therefore sufficient to explain the current magnitude. The production invariant remains one common coordinate/stencil/weight/arithmetic contract between the C++ matrix builder and live Fortran for all 13 shared scalar/composite choices. A new near-gridline live basis/spreading point passes that invariant in 2D and 3D; no production formula or matrix copy is added. |
| 2026-08-20 | Use an operator-controlled sandbox replay only to separate multiplicative-application arithmetic from independently reconstructed-operator arithmetic; retain both comparisons. | An external wrapper maps the raw live IBAMR `A`, `E_h`, RHS/seed, and probe vectors into canonical sandbox order and then invokes the unmodified pinned multiplicative functions. The mapped operators are exact by construction and are not independent construction evidence. At `K=1`, smoother/V-cycle/solution/live-original-residual errors are `1.22e-15`/`1.66e-15`/`2.78e-16`/`2.31e-14`; at `K=1e2`, `4.03e-14`/`2.32e-14`/`1.24e-14`/`2.36e-12`; iteration counts remain exactly `9` and `8`. At `K=1e4`, however, the same-operator first mismatch is already the multiplicative sweep: smoother correction `4.06e-11`, V-cycle correction `7.61e-11`, solution `2.58e-12`, live original residual `2.80e-10`, and exact 9 iterations. Thus operator reconstruction is not the only `K=1e4` seam; cross-runtime local factor/solve and accumulated residual-update arithmetic must be traced. The independently reconstructed-operator comparison stays red and remains mandatory. |
| 2026-08-20 | Close the missing local-solve-observer evidence gap without retaining an unfactored local matrix or broadening production indexing. | The backend-neutral observer accepts a selection predicate. The LU backend checks it before diagnostics, reconstructs only the selected local submatrix transiently from the borrowed global operator, wraps the already-copied selected RHS and solution workspace in transient sequential vectors, calls the observer, and destroys the diagnostic objects. Its native fixture verifies exactly one selected callback, dimensions, finite/tiny backward error, deallocation, and reinitialization. A normalized expected-failure fixture verifies that an exactly singular local factorization reports subdomain 0 and LAPACK `info=3`; the two focused fixtures pass **2/2**. This is implementation validation, not an audit result. |
| 2026-08-20 | Classify the first controlled `K=1e4` local mismatch as conditioning-sensitive cross-library forward error, not a matrix, patch-order, sign, or LU-specific semantic defect. | Patch 0 is exact through its local solve. Patch 236, the first IB-coupled patch, has exact mapped membership and matrix, `7.11e-15` RHS difference, and `3.14e-11` LU solution difference. Replaying the IBAMR RHS with the oracle's dense LU changes the solution by only `2.22e-16`; both solves have tiny backward error, while the exact mapped matrix has 2-norm condition number `2.30e6` and smallest singular value `8.46e-1`. On the same patch PETSc SVD differs from the oracle by `3.25e-11`, so switching the production backend does not remove the seam. On later patch 364, accumulated RHS drift is `1.89e-9`; its `8.53e-12` LU solution difference separates into `3.56e-12` RHS effect and `7.05e-12` same-RHS backend effect. No acceptance criterion is relaxed and no solver is changed merely to tune cross-runtime roundoff. |
| 2026-08-19 | Reject PETSc PtAP algorithm selection as a remedy for the remaining same-geometry seam. | The supported serial `rap` path produces byte-for-byte identical `A`, `E_h`, smoother, V-cycle, solution, and original-residual exports to the default path. The installed serial implementation exposes `scalable`, `rap`, and `hypre`; `nonscalable` is not a valid choice. No option change or additional matrix is retained. |
| 2026-08-19 | Accept the live FAC construction seam only after a two-level implicit-IB native run exercises nonzero `d_SAJ_mat` with `PRESSURE_CELL` and multiplicative shell application. | The corrected path passed through the native runner with `test_failures = 0`, three Krylov iterations, and finite residual `3.67879e-06`; this establishes live plumbing/application, not sandbox numerical parity or an independent audit result. |
| 2026-08-19 | Define the reduced-system IB residual as an independently evaluated live coupling-consistency defect, not a nonexistent third Eulerian row block. | The formal oracle and IBAMR velocity-path Jacobian both solve `x=(u,p)`: total momentum and divergence are the only reduced residual rows. Comparing the matrix-free IB action with FAC-owned live `d_SAJ_mat x` provides a real IB coupling defect while preserving original-operator semantics and avoiding a surrogate or fabricated Lagrangian state. |
| 2026-08-20 | Complete the selected-patch residual-update trace before considering any change to multiplicative composition. | The passive observer now also lends the current global pre-update residual; the exporter writes that vector and a transient embedded patch correction for selected patches and their immediate successor. Restrictions to every selected local RHS are exact. At `K=1e2`, first coupled patch 236 has condition number `5.39e3`, `2.53e-14` local-solution difference, and a `2.50e-12` exact-matrix correction effect; PETSc-versus-MATLAB multiplication order contributes `9.09e-13`, and actual patch-237 pre-residual difference is `2.05e-12`. At `K=1e4`, the corresponding values are `2.30e6`, `3.14e-11`, `1.36e-8`, `6.80e-10`, and `1.36e-8`. Patch 237 sees exactly the actual updated residual in both implementations. The mismatch is therefore fully localized to backward-stable local-solve forward sensitivity plus sparse-multiply accumulation order, not stale residuals, patch ordering, restriction, sign/gauge handling, or a hidden correction-composition formula. Fixed low-`K` acceptance remains red and is not waived. |
| 2026-08-20 | Apply the fixed low-`K` parity table to a two-source common-arithmetic algorithm replay and make live production residuals a separate mandatory gate. | The completed first-mismatch trace proves that the remaining entrywise cross-runtime drift comes from backward-stable vendor/runtime local solves and sparse-multiply accumulation order, not a semantic defect. Both candidate and oracle algorithms must therefore be replayed in one arithmetic runtime on each independently exported operator/RHS, while consuming their own exported ordered construction artifacts. Raw live-versus-MATLAB differences remain reported diagnostics. The actual IBAMR path must separately converge at the declared `1e-10` original-relative-residual tolerance and report/compare iteration history plus momentum, divergence, and IB consistency components. Fixed low-`K` tolerances are not made conditioning-aware or waived. |
| 2026-08-20 | Make the original Fortran interaction kernels authoritative for every matrix-built IB kernel, without refactoring correct matrix-free paths into a new shared ABI during CAV recovery. | The existing scalar Fortran functions are callable for several kernels, but optimized IB4/IB5 routines expand their own recurrences and IB6 has no scalar helper. Calling full kernels to recover rows would require scratch arrays and repeated setup calls; adding a shared weight-vector ABI would modify the matrix-free hot path and require separate inlining/performance evidence. The C++ matrix callbacks remain small faithful ports locked by exhaustive 2D/3D interpolation/spreading tests. The only live edit is the independently reproduced 3D IB5 stale-index defect. |
| 2026-08-20 | Record the disposable multiplicative reconstruction as green under the selected two-track contract, without promoting the durable M5 or any independent-audit gate. | Schema `cav-common-arithmetic-multiplicative-v3` replays candidate and oracle semantics in MATLAB R2025b on both independently exported operator/RHS sources. Ordered patch structure is exact, and candidate-versus-oracle correction/residual/history results are bitwise equal through one sweep, the two-level V-cycle, and FGMRES for `K=1`, `1e2`, and `1e4`. Separate live IBAMR runs meet the `1e-10` original-relative-residual and IB-consistency criteria and match oracle iteration counts and normalized histories within the fixed table. Cross-runtime replay-versus-export differences remain explicit non-gating diagnostics; at `K=1e4` the raw final-residual difference is `1.13e-10`, consistent with the already-localized runtime arithmetic seam. Harness/report hashes are recorded in the current stopping point. |
| 2026-08-20 | Make same-unit native CI tests and IBAMR comment/style conformance binding production and audit criteria for every prospective durable CAV feature unit. | Each feature commit/PR must carry deterministic, `attest`-discoverable, CI-suitable focused tests; shared executables with separate fixtures are preferred, runtime budgets and measurements are recorded, and external parity cannot replace native regression coverage. Lane N owns actually running and timing the fast set, assigned technical lanes verify their documented invariants, simplicity challenges harness/API/comment excess, and synthesis checks the feature/test PR boundary. This strengthens prospective CAV units without reopening historical-ledger Finding 1. |
| 2026-08-20 | Correct the native CI authority and focused-versus-broad scope to IBAMR's repository-root `attest` workflow. | The user's authoritative broad example is `/path/to/attest --mpi-executable /opt/homebrew/bin/mpiexec --numdiff-executable /opt/homebrew/bin/numdiff --verbose -R IB/`, with `-R IB/` selecting all IB examples/tests. CMake builds executables and may arrange fixtures, but no other runner or target-only result establishes CI integration. Focused feature cases and the broad regression are labeled and recorded independently; focused success cannot substitute for a required broad run. Lane N records exact paths/selectors/cases/fixtures/status/runtime/exclusions, and synthesis verifies both scopes and same-unit boundaries. |
| 2026-08-20 | Recover the first durable production unit as a borrowed live-state view in full coupled velocity-pressure global numbering, with no new matrix or algebraic-object copy. | `cav-1a` exposes the solver-owned live PETSc operator handle, existing global-numbered locally owned velocity/pressure DOF vectors, ownership counts, operator provenance, nullspace declarations/attachment, and pressure-gauge behavior through one invalidatable nonowning view. It does not extract blocks, create a second operator representation, translate to a local subdomain numbering, serialize data, construct patches, or change solver behavior. The focused same-unit fixture proves identity, partition, lifecycle, and gauge semantics; exporter/mapping work remains `cav-1b`. |
| 2026-08-20 | Record the first durable unit's focused and broad `attest` scopes without treating an environmental limitation as a pass or silently substituting the focused case. | The focused live-view case passes twice and detects a disabled-hook mutation. The exact broad `-R IB/` command discovers 18 cases and yields 14 serial passes plus four rank-4 failures at Open MPI TCP-listener setup in the managed sandbox; the separately labeled serial subset passes 14/14. This permits a local feature commit with an explicit evidence limitation, but full review readiness retains the rank-4 rerun obligation. |

## 14. Superseded stopping-point snapshot

- **Active branch/worktrees:** authoritative plan on `codex/implicit-ib-cav-0-recovery-plan` in `/Users/boyceg/.codex/worktrees/0760/IBAMR`; disposable implementation on `codex/implicit-ib-cav-integration-reconstruction` in `/tmp/ibamr-cav-integration-reconstruction`. Both started at exact foundation `0c8447e3593113a8b62f94df9482c60d896f0d3b`. The clean oracle clone is detached at `5b77344db6746269f8c77695c99e9043907ba74b` in `/tmp/implicit-ib-sandbox-cav-oracle-5b77344`.
- **Last green validation:** disposable integration Debug build in `/tmp/ibamr-cav-integration-build-dbg`: the complete pre-existing serial interpolation/spreading suite passes **126/126**, including the corrected 3D IB5 spreading fixture; the focused implicit-IB construction, shell-smoother, local-LU, and live two-level FAC selection passes **7/7**. The exhaustive matrix/live kernel test passes **1/1 in 2D and 1/1 in 3D** for all 13 shared choices, comparing every stored matrix-row weight with live Fortran interpolation basis probes and periodic-DOF-folded, volume-scaled live spreading at generic, center, half-cell, near-gridline, lower-boundary, and upper-boundary points. Every choice proves that a wrapped stencil was exercised and that the captured-row formula and unequal-weight order mutations are detected. Maximum differences are `1.30e-13` in 2D and `7.23e-14` in 3D, both within the unchanged `1024 epsilon_machine` factor bound; IB4 is within `2.78e-16` in 2D and `1.53e-16` in 3D. Fresh passive-trace live parity at `K=1`, `1e2`, and `1e4` confirms exact configuration, mapping, patch seed/order/membership, and iteration counts. With independently generated but nominally identical geometry, the mapped smoother/V-cycle/original-physical-residual errors are respectively `7.77e-16`/`1.64e-15`/`2.58e-14`, `3.86e-14`/`2.35e-14`/`2.39e-12`, and `6.80e-13`/`8.32e-13`/**`2.69e-10`**. The gauge-aligned final-solution errors are `2.88e-16`, `1.07e-14`, and `1.49e-12`; iteration counts are exactly 9, 8, and 9. Original IBAMR physical residual components `(momentum, divergence, independently evaluated IB consistency)` are `(2.389e-10, 3.925e-10, 5.799e-14)`, `(3.414e-9, 4.857e-9, 7.704e-12)`, and `(4.459e-8, 6.135e-8, 7.952e-10)`. With exact exported live coordinates, the smoother/V-cycle/solution/original-residual errors are `6.66e-16`/`1.59e-15`/`2.78e-16`/`2.49e-14`, `4.20e-14`/`2.25e-14`/`1.45e-14`/`2.47e-12`, and `4.92e-13`/`1.00e-11`/`1.12e-12`/`2.40e-10`, with exact iteration counts `9`, `8`, and `9`. Thus exact shared geometry leaves the `K=1e2` fixed residual target and `K=1e4` fixed residual ceiling red. The new IB4 operator, smoother, V-cycle, solution, and original-residual export files are byte-for-byte identical to the pre-all-kernel candidate at all three stiffnesses. Both `IBTK2d` and `IBTK3d` and the affected IBAMR libraries build cleanly; the exact-foundation expanded selection remains **69/69**. Raw bundles and generated reports remain under `/tmp`. These are implementation validations, not independent audit results.
- **Current work:** M0 is complete and the mandatory `K=1`, `1e2`, `1e4` problem/operator/base-patch, finest-level smoother, transfer/coarse, complete two-level V-cycle, FGMRES, and original physical residual comparisons have been freshly rerun against the all-kernel candidate. `K=1` is green. Under the plan's current per-`K` table, `K=1e2` has green smoother/V-cycle/solution results but its exact-shared-coordinate `2.47e-12` entrywise original-residual error exceeds the fixed `1e-12` target; it is a cross-runtime factor/composition mismatch, not an independent-geometry or internal matrix/live formula mismatch, and must not be described as a completed numerical gate. `K=1e4` remains unpromoted because its mandatory matrix-free entrywise residual misses both the `2e-11` target and `1e-10` final ceiling, including with exact shared coordinates. The maximum mapped IB4 coefficient differs from MATLAB by one binary64 unit at identical fractional coordinates; a controlled same-order replay shows that this one coefficient alone can produce `9.09e-13` and `2.33e-10` `E_h` entry changes at `K=1e2` and `1e4`, while all coefficient differences can produce `3.64e-12` and `4.66e-10`. An external operator-controlled replay using exact mapped live `A`, `E_h`, RHS/seed, and probe vectors is green through the full multiplicative path at `K=1`; at `K=1e2` its smoother/V-cycle/solution errors are green but its `2.36e-12` final live residual still misses the fixed `1e-12` target. At `K=1e4`, its first same-operator mismatch is the finest multiplicative sweep (`4.06e-11` correction), followed by `7.61e-11` V-cycle correction and `2.80e-10` live final residual with exact 9 iterations. The disposable minimal real-scalar serial LU backend uses canonical full global numbering outside the local factor boundary: it borrows the full operator; stores only one overwritten dense factor buffer, pivots, patch references, and reusable vector workspaces per patch; and performs incremental residual evolution with reusable global vectors and the live matrix. It does not copy the full operator, a velocity block, a transpose, an unfactored local matrix, or residual-update matrices. Its positive native semantics and lifecycle case is green; a stable native singular-factorization failure regression and BLAS/LAPACK-backend local-solve observer callbacks remain open. The same-geometry `K=1e4` comparison makes the construction first mismatch more precise: mapped `J` structure is exact with `2.78e-17` maximum entry error, the force Jacobian is exact, and the mapped `E_h` structure is exact, but double-precision composition produces `2.00e-10` mapped assembled-IB action and `2.29e-10` live-versus-assembled action effects on the exact seed. PETSc's supported `rap` PtAP path is byte-for-byte identical to the default; explicit sparse product associations did not improve oracle agreement and would add an unjustified temporary matrix. The all-kernel `cav-1c` disposable implementation is green for identical-coordinate 2D/3D interpolation and spreading, actual periodic wrapping, controlled formula/order sensitivity, the full existing serial suite, and unchanged end-to-end IB4 behavior. Its contract covers every shared matrix/live kernel—piecewise linear, B-splines 3--6, composite B-splines 32/43/54/65, and IB kernels 3--6—with the live coordinate, stencil, one-dimensional weight, tensor-product, and traversal formulas treated as the production source of truth. Only its durable-stack port remains open. No acceptance criterion is waived, and no copied production matrix or alternate accepted operator has been introduced. Observer lifecycle/regrid tests, durable-stack decomposition/port, and independent audit remain open. The constructor/exporter and LU backend remain explicitly serial-only; construction makes no sandbox-STRICT claim. No auditor, audit coordinator, or synthesizer has been launched, and no independent gate is represented as complete.
- **Blockers:** none.
- **Historical next action (completed):** Add backend-neutral local-solve observer callbacks to the disposable BLAS/LAPACK LU path, with no retained unfactored matrix, and use them to export the same selected `K=1e4` patch matrices/RHS/solutions and per-patch sweep states as the oracle so the same-operator multiplicative first mismatch can be localized before changing solve or residual-update arithmetic.

## 15. Superseded stopping-point snapshot

- **Active branch/worktrees:** the authoritative plan remains on `codex/implicit-ib-cav-0-recovery-plan` in `/Users/boyceg/.codex/worktrees/0760/IBAMR`; the disposable implementation remains on `codex/implicit-ib-cav-integration-reconstruction` in `/tmp/ibamr-cav-integration-reconstruction`. Both started at exact foundation `0c8447e3593113a8b62f94df9482c60d896f0d3b`. The clean oracle worktree remains detached and unmodified at `5b77344db6746269f8c77695c99e9043907ba74b` in `/tmp/implicit-ib-sandbox-cav-oracle-5b77344`.
- **Last green validation:** the exhaustive all-kernel matrix/live consistency cases pass in 2D and 3D for all 13 shared choices, with maximum differences `1.30e-13` and `7.23e-14`; all formula and unequal-weight-order sensitivity controls activate, and the complete existing serial interpolation/spreading suite remains green **126/126**. The previously selected implicit-IB/FAC validation remains **7/7** and the exact-foundation selection remains **69/69**. After the local observer addition, the positive LU observer/reinitialize fixture and the normalized singular-factorization expected-failure fixture pass together **2/2**. Both affected libraries and the focused native executable build cleanly. Generated traces and reports remain under `/tmp`. These are implementation validations, not independent audit results.
- **Current work:** all matrix-builder-supported kernels now use the same coordinate, stencil, one-dimensional weight, tensor-product, and traversal formulas as their live 2D/3D Fortran paths; no matrix copy or persistent index translation was added. The selected-only LU observer reconstructs an unfactored matrix only transiently after its predicate selects a patch. In the exact-operator `K=1e4` trace, patch 236 is the first coupled patch and has exact mapped membership and matrix, `7.11e-15` RHS drift, `2.30e6` 2-norm condition number, `3.14e-11` LU forward difference, and tiny backward errors. PETSc SVD has a comparable `3.25e-11` forward difference, so the mismatch is conditioning-sensitive cross-library arithmetic rather than an LU-only semantic defect. The fixed low-`K` criteria are unchanged: `K=1e2` still misses its `1e-12` entrywise final-residual target, and `K=1e4` still misses the `1e-10` hard final-output ceiling. The durable `cav-1c`/backend ports, FAC observer lifecycle/regrid tests, and all independent audits remain open. No auditor or synthesizer has been launched.
- **Blockers:** none.
- **Historical next action (completed):** add a selected-patch, read-only multiplicative-state observation point in the disposable integration path and export the global pre-update residual plus embedded patch correction at patch 236, then replay that one update externally with the oracle local solution to determine whether any residual-update semantic mismatch remains after the now-explained conditioning-sensitive local solve difference.

## 16. Superseded stopping-point snapshot

- **Active branch/worktrees:** the authoritative plan is on `codex/implicit-ib-cav-0-recovery-plan` in `/Users/boyceg/.codex/worktrees/0760/IBAMR`; the first durable feature unit is on `codex/implicit-ib-cav-1a-live-operator-audit-schema` in `/tmp/ibamr-cav-durable-1a`; and the disposable implementation remains on `codex/implicit-ib-cav-integration-reconstruction` in `/tmp/ibamr-cav-integration-reconstruction`. All derive from exact foundation `0c8447e3593113a8b62f94df9482c60d896f0d3b`; `cav-1a` was created from committed plan tip `b8a69bf23` and will be fast-forwarded to the current forward-only plan update before its feature commit. The clean oracle worktree remains detached and unmodified at `5b77344db6746269f8c77695c99e9043907ba74b` in `/tmp/implicit-ib-sandbox-cav-oracle-5b77344`.
- **Last green validation:** every matrix-builder-supported kernel passes the shared matrix/live 2D/3D formula contract, including wrapped and near-gridline cases and per-kernel formula/order sensitivity controls; maxima remain `1.30e-13` and `7.23e-14`, and the existing serial interpolation/spreading suite remains **126/126**. The latest passive observer change builds cleanly and its positive observer/reinitialize and singular-factorization expected-failure fixtures pass **2/2**. The earlier focused implicit-IB/FAC selection remains **7/7** and the exact-foundation selection remains **69/69**. Generated traces and reports remain under `/tmp`. These are implementation validations, not independent audit results.
- **Current work:** the all-kernel contract, global velocity-pressure indexing, borrowed-operator ownership, no-copy local-factor design, local failure path, and selected-patch state observation are implemented in the disposable reconstruction. Exact-operator traces at `K=1e2` and `K=1e4` prove exact patch order/membership, exact local restriction, exact local matrices, correct correction embedding, and correct successor-patch state. At the first coupled patch, cross-library backward-stable solve differences are `2.53e-14` and `3.14e-11`, amplified by the operator to `2.50e-12` and `1.36e-8`; PETSc-versus-MATLAB multiply order contributes `9.09e-13` and `6.80e-10`. No stale-residual, sign, gauge, patch-order, restriction, or hidden composition-formula defect remains at this first mismatch. Switching to PETSc SVD does not remove it. The fixed `K=1e2` entrywise residual target and `K=1e4` hard final-output ceiling remain red and are not waived. RAS has not started, and no auditor or synthesizer has been launched.
- **Blockers:** the mandatory low-`K` entrywise cross-runtime gate is still red after the full first-mismatch ladder localized the difference to vendor/runtime floating-point solve and sparse-multiply order, with no production-justified semantic, backend, index, or matrix-copy change remaining. Current instructions forbid conditioning-aware relief at these stiffnesses, and no commits/durable-stack branch construction are authorized.
- **Historical next action (completed):** present the acceptance-contract decision with this reproduced evidence—retain exact cross-runtime entrywise arithmetic as a hard blocker, or explicitly authorize a revised common-arithmetic algorithm-parity criterion while keeping live physical residuals separate—without changing code, starting RAS, or marking the multiplicative gate complete.

## 17. Superseded stopping-point snapshot

- **Active branch/worktrees:** the authoritative plan remains on `codex/implicit-ib-cav-0-recovery-plan` in `/Users/boyceg/.codex/worktrees/0760/IBAMR`; the disposable implementation remains on `codex/implicit-ib-cav-integration-reconstruction` in `/tmp/ibamr-cav-integration-reconstruction`. Both derive from exact foundation `0c8447e3593113a8b62f94df9482c60d896f0d3b`. The clean oracle worktree remains detached and unmodified at `5b77344db6746269f8c77695c99e9043907ba74b` in `/tmp/implicit-ib-sandbox-cav-oracle-5b77344`.
- **Last green validation:** every matrix-builder-supported kernel passes the Fortran-authoritative 2D/3D matrix/live contract, including wrapped and near-gridline cases and per-kernel formula/order sensitivity controls; maxima remain `1.30e-13` and `7.23e-14`, and the existing serial interpolation/spreading suite remains **126/126**. The positive observer/reinitialize and singular-factorization expected-failure fixtures pass **2/2**; the earlier focused implicit-IB/FAC selection remains **7/7** and the exact-foundation selection remains **69/69**. A disposable symbol/formula probe confirmed that callable scalar Fortran IB4/IB5 helpers agree with their optimized recurrence weights within `2.22e-16`, but also confirmed the architectural inventory: optimized live routines do not universally call those helpers and IB6 exposes no scalar helper. Generated probes, traces, and reports remain under `/tmp`. These are implementation validations, not independent audit results.
- **Current work:** the user selected the common-arithmetic interpretation of low-`K` parity. The plan and numerical-audit contract now require both candidate and oracle algorithms to be replayed in one documented arithmetic runtime on both independently exported operator/RHS sources, using each side's own exported patch/macro/transfer definitions. The live IBAMR path is a separate mandatory production check with original relative residual, iteration history, and momentum/divergence/IB consistency components. Raw cross-runtime entrywise differences remain visible diagnostics but are no longer mistaken for semantic parity. The original Fortran interaction kernels are now explicit authority for all matrix callbacks; correct matrix-free paths are not refactored or edited merely for sharing. Direct full-kernel matrix construction and a new shared Fortran weight ABI are rejected for this recovery because they add scratch/copy/setup cost or change the hot path without need. RAS has not started, and no auditor or synthesizer has been launched.
- **Blockers:** no algorithmic decision is blocked. Durable-stack branch construction and commits remain unauthorized, so current implementation work stays in the disposable reconstruction and generated evidence stays under `/tmp`.
- **Historical next action (completed):** implement the two-source common-arithmetic multiplicative replay from the raw live candidate/oracle exports and their independently exported ordered patch artifacts, then run the fixed `K=1`, `1e2`, and `1e4` comparison without changing production solver arithmetic or starting RAS.

## 18. Current stopping point

- **Active branch/worktrees:** this forward-only evidence update on `codex/implicit-ib-cav-0-recovery-plan` in `/Users/boyceg/.codex/worktrees/0760/IBAMR` is based on `89b6da10e259342f6383c0caa2d619a2e72f0a58`; the first durable unit is an uncommitted candidate based on that same tip on `codex/implicit-ib-cav-1a-live-operator-audit-schema` in `/tmp/ibamr-cav-durable-1a` and will be fast-forwarded to the evidence update before commit; and the disposable implementation remains on `codex/implicit-ib-cav-integration-reconstruction` in `/tmp/ibamr-cav-integration-reconstruction`. All derive from exact foundation `0c8447e3593113a8b62f94df9482c60d896f0d3b`. The clean oracle worktree remains detached and unmodified at `5b77344db6746269f8c77695c99e9043907ba74b` in `/tmp/implicit-ib-sandbox-cav-oracle-5b77344`.
- **Last green validation:** the uncommitted durable `cav-1a` candidate compiles through both 2D and 3D library paths. Its exact focused `attest` case is discovered and passes twice in 1.96 s and 0.48 s; disabling only the new verification hook in a disposable fixture is detected. The exact broad `-R IB/` command completes in 82.93 s with all 14 serial cases passing and four rank-4 cases blocked before IBAMR by the managed sandbox's Open MPI TCP-listener restriction; the explicitly labeled serial subset reruns green **14/14** in 80.75 s. Separately, disposable schema `cav-common-arithmetic-multiplicative-v3` passes at `K=1`, `1e2`, and `1e4`. Candidate and oracle results are exactly equal in the common runtime for ordered patch structure, one multiplicative sweep correction/fresh residual, two-level V-cycle correction/fresh residual, FGMRES solution/fresh residual, and relative iteration history on both independently exported operator/RHS inputs. Separate live IBAMR runs converge in exactly 9, 8, and 9 iterations with original relative residuals `3.0556e-11`, `9.1779e-11`, and `1.2450e-11`; their live matrix-free-versus-`SAJ` IB-action relative defects are `2.9759e-15`, `3.9531e-15`, and `4.0806e-15`. The live scalar-summary errors against the oracle are `1.97e-16`, `2.12e-14`, and `2.30e-12`, and normalized iteration-history errors are `2.26e-17`, `4.41e-16`, and `1.14e-15`, all within the fixed table. Raw cross-runtime final-residual differences remain explicit diagnostics (`1.78e-14`, `1.34e-12`, `1.13e-10`) and do not replace either mandatory track. Fresh focused native validation is green for both all-kernel programs (`test_failures = 0`; maxima `1.30e-13` in 2D and `7.23e-14` in 3D) and for all 17 applicable serial implicit-IB/local-backend cases. The excluded two-rank expected-error fixture cannot reach its intended assertion because this sandbox rejects the MPI TCP listener. Harness SHA-256 is `9d6b1800b53e278786ed3430bbcfa2190de848fa5e503deab7302d82dbeeb49f`; report hashes for increasing `K` are `4912a3ad04045a77617468ce8e5c9af436040f700a69da216564af601841e2b7`, `826147134cf4100815fe0f47f4dd44b8c9cc3169df56707eed331df24c1b25f2`, and `7f6785c1b8a5f8af73223006ee3fd6bf7da1e636e9fc2146800e825278f12e4a`. Generated evidence remains under `/tmp`. These are implementation validations, not independent audit results.
- **Current work:** `cav-1a` now contains only a public, nonowning, lifecycle-bounded view of the solver-owned live operator state and a focused same-unit case in the existing level-solver executable. It preserves full coupled global numbering, borrows the exact KSP operator, references existing ownership/field vectors, adds no matrix or block copy, and does not serialize, construct patches, add a backend/FAC observer, or alter solver behavior. Both dimensions compile; focused discovery, deterministic execution, and sensitivity are green. The broad result and environment limitation are recorded separately. RAS has not started, and no auditor, audit coordinator, or synthesizer has been launched.
- **Blockers:** no blocker to the authorized local `cav-1a` commit. Full review readiness still requires the four discovered rank-4 broad cases to be rerun outside this TCP-restricted sandbox; the focused case and all 14 serial broad cases are green.
- **Exactly one next action:** commit this forward-only `cav-0` evidence update, fast-forward `cav-1a` to it, and commit the already-validated same-unit live-view production/test diff without pushing or launching an audit.
