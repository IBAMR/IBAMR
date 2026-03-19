# Subdomain Index Cache Audit

This note records the current indexing/cache structures used in Stokes
subdomain setup and highlights places where more contiguous representations may
be preferable.

## Current hotspot

The main structure of interest is
[`StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData`](/Users/boyceg/code/IBAMR/include/ibamr/StaggeredStokesPETScMatUtilities.h),
which currently stores:

- `velocity_dof_to_adjacent_cell_dofs` as `std::map<int, std::set<int>>`
- `cell_dof_to_closure_dofs` as `std::map<int, std::set<int>>`
- `velocity_dof_to_component_axis` as `std::map<int, int>`
- `velocity_dof_to_paired_seed_velocity_dofs` as `std::map<int, std::set<int>>`

These choices are convenient for correctness and deduplication, but they are
not ideal for data locality.

## Audit findings

### 1. Some maps are good candidates for `vector<>`

The best candidates are maps built for all locally relevant keys in a
contiguous range:

- local ownership counts and per-local-DOF flags already use vectors in
  several places
- `construct_nonoverlap_subdomains_from_overlap()` currently uses
  `std::map<int, std::size_t> dof_owner`, but this is built over the locally
  owned DOF range and could use a dense vector indexed by local offset

This is the lowest-risk contiguous rewrite.

### 2. `PatchLevelCellClosureMapData` likely needs an explicit local indexing layer

The current closure maps use global DOF indices as keys and values. That makes
the data easy to reason about, but a direct replacement of `map<>` with
`vector<>` is not clean unless we first define local indexing.

A better representation would likely be:

- one dense local index space for locally relevant velocity DOFs
- one dense local index space for locally relevant cell/pressure DOFs
- translation tables between local ordinals and global DOF ids
- vector-backed adjacency/closure lists keyed by local ordinal

That would let us store:

- velocity component axis as `vector<int>`
- adjacent cell DOFs as `vector<vector<int>>`
- cell closures as `vector<vector<int>>`
- paired STRICT seed groups as `vector<vector<int>>`

while still exposing global DOF ids at the algorithm boundaries.

### 3. `set<int>` is often acting as “sort + deduplicate”

Several temporary sets appear to be used mostly for:

- duplicate elimination
- sorted iteration order

Examples:

- `axis_velocity_dofs`
- `seed_velocity_dofs`
- `initial_velocity_dofs`
- `closure_dofs`

For many of these, a more locality-friendly pattern would be:

1. append into `vector<int>`
2. `std::sort(...)`
3. `erase(std::unique(...), ...)`

This is likely a good fit for build-once caches, especially if insertions are
currently dense and bulk.

### 4. Some sets are still semantically useful

Not every `std::set<int>` should automatically become a vector.

Places that still need scrutiny:

- structures that are repeatedly unioned from multiple sources
- places where the algorithm depends on set membership tests during
  construction, not just after finalization
- overlap/nonoverlap subdomain outputs, which currently match the surrounding
  API expectations

So the right approach is selective replacement, not a blanket rewrite.

## Recommended phases

### Phase 1: Low-risk contiguous rewrites

- replace `dof_owner` in
  `construct_nonoverlap_subdomains_from_overlap()` with a dense local-owner
  vector
- replace temporary build-only `set<int>` uses with `vector<int> + sort/unique`
  where the code only needs final sorted unique results

### Phase 2: Closure-map representation redesign

Introduce a new local-indexed closure cache type with:

- `first_local_velocity_dof`, `first_local_cell_dof` or equivalent translation
  metadata
- dense ordinal-to-global and global-to-ordinal lookup support
- vector-backed adjacency/closure storage

This should be done as a new representation first, with side-by-side
validation, before deleting the old `map/set` implementation.

### Phase 3: Broader audit of subdomain setup caches

Once the closure maps are settled, audit related subdomain setup code for the
same patterns:

- repeated global-to-owner lookups
- build-only ordered sets
- sparse-key maps that are actually dense over a local contiguous range

## Immediate takeaway

Yes, we should adopt these patterns more broadly, but the best first step is
not to rewrite `PatchLevelCellClosureMapData` in place. The best first step is:

1. convert the obvious local-range maps and temporary build sets
2. then redesign the closure-map cache around explicit local ordinals

That gives us early locality wins without making the first change too risky.
