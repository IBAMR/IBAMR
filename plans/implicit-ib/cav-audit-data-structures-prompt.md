# Audit Lane D: Data-Structure Efficiency

**Prompt version:** `cav-audit-v8/D` (2026-08-21)

Use with `cav-independent-audit-prompt.md`. This lane is independent and permanently read-only with respect to the candidate, oracle, comparator, and pinned refs.

Use the common protocol and implementation-plan Section 7 for the canonical focused-versus-broad `attest` evidence contract; do not redefine or substitute those scopes in this lane.

## Question

Are the candidate's data structures, ownership models, and lifecycle rules efficient and scalable for their declared setup-time and apply-time responsibilities?

## Required evidence

- Inventory every material structure: creator, lifetime, owner, cardinality, asymptotic storage, consumers, mutation/invalidation, and serial/distributed meaning.
- Separate construction/regrid/setup work from per-sweep, per-level, per-macro, and per-Krylov-iteration work.
- Derive deterministic counts or tight bounds for allocations, copies, lookups, conversions, sorting/set work, recomputation, synchronization, and communication-sensitive ownership operations. Relate each count/bound to patches, macros, levels, DOFs, or ranks.
- Audit patch/macro overlap and owned-set representation, local matrix/factor storage, DOF maps, cached residual-update data, and raw-export state.
- Apply the common protocol's W1/W2 and, at M9, W3 structural scenarios to make setup/apply counts and ownership/storage growth concrete. These are inspection scenarios, not timing benchmarks.
- Exercise initialize/deallocate/reinitialize cycles, hierarchy regridding, cache invalidation, stale-source rejection, and retained-memory/resource behavior.
- Verify rank-local versus replicated data and advertised MPI ownership. For RAS, show that owned sets partition the global supported scope and that distributed design does not require repeated full-domain metadata or reconstruction.
- Cover empty/smallest-legal/boundary/no-IB and ownership-boundary cases where data cardinality or ownership changes. In 2D all SAMRAI hierarchy patches have extent at least `4x4`; use the actual neighboring dimensional rule in 3D.
- Inspect specifically for padding records, synthetic splits or hierarchy levels, alternate initialization, compatibility shims, or duplicate ownership whose purpose is to admit sub-minimum hierarchy patches. Those structures are outside supported scope and should be removed; invalid requested geometry belongs to normal SAMRAI/IBAMR validation. Do not confuse a small algebraic CAV patch with a hierarchy patch.
- Activate a disposable forced per-apply rebuild, duplicated full-domain ownership record, or equivalent recomputation/ownership sensitivity challenge and show that the audit detects it relative to an unmutated control.
- Verify that ownership, borrowing, lifetime, invalidation, regrid, supported-scope, and representation comments match actual state transitions and do not conceal mutable or duplicated internals.

## Exclusions

- Do not report container taste or naming as an efficiency defect.
- Do not conflate setup cost with apply cost.
- Do not demand distributed machinery before the distributed milestone.
- Do not require data structures for unsupported sub-minimum hierarchy patches.
- Do not propose speculative micro-optimization without scale or repeated-work evidence.
- Do not run or require production timing, quote speedups/throughput, or treat elapsed time gathered during concurrent audits as evidence. Focus on code structure, asymptotic work, and deterministic counts/bounds.
- Do not edit or fix the candidate.

## Pass/fail

**PASS** requires clear ownership/lifetime/invalidation, phase-appropriate work, complete structural-scenario analysis, supported-scope lifecycle coverage, and successful sensitivity detection, with no material scaling/allocation/lookup/recomputation/resource defect established by code and deterministic counts/bounds. A demonstrated material repeated cost, ambiguous ownership, stale cache, resource leak, or per-apply reconstruction that belongs in setup is **FAIL**. Missing structure/count/bound evidence, unsupported-scope declarations, lifecycle cases, or inactive sensitivity controls is **INCOMPLETE**. Missing timing is never a gap.

## Lane-specific report schema

After the common schema, include: structure/lifecycle table; setup-versus-apply operation counts/bounds; W1/W2/W3 structural analysis as applicable; storage-growth evidence; ownership/MPI analysis; regrid/reinitialize/invalidation/resource results; sensitivity control; concrete improvement candidates with structural justification; distributed-risk boundary; and exact-SHA `PASS`, `FAIL`, or `INCOMPLETE`.
