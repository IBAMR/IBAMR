# SAMRAI Compatability Migration Plan

## Overall Plan

### Goal

Port IBAMR from IBSAMRAI2-era SAMRAI usage toward SAMRAI 4.5.0 with minimal churn, while keeping the tree buildable at each step.

### Strategy

1. Establish a complete compatibility surface first.
2. Keep include and type migrations largely mechanical.
3. Defer behavioral/API redesign work until after alias/interposition rollout.
4. Land work in buildable batches.

### Current Architecture

- Canonical compatibility environment:
  - `ibtk/include/samrai_compatibility/samrai_compatibility_environment.h`
- Forwarding headers for existing include sites:
  - `ibtk/include/ibtk/samrai_compatibility_names.h`
  - `ibtk/include/ibtk/samrai_compatibility_legacy_aliases.h`
- Alias headers:
  - `ibtk/include/samrai_compatibility/<module>/SAMRAI*.h`
- Legacy include wrappers:
  - `ibtk/include/samrai_compatibility/legacy/...`

## What Is Done

### Compatibility Surface (Complete for One-to-One Matches)

- One-to-one SAMRAI alias headers are now in place for all classes/types that are consistent between IBSAMRAI2 and SAMRAI 4.5.0.
- Legacy wrapper coverage is complete for the alias set (no missing alias-to-legacy wrapper pairs).
- Compatibility include paths are exported in CMake so flat `#include <SAMRAI*.h>` works in build/install contexts.

### Namespace and Pointer Compatibility

- `Pointer` compatibility behavior has been established via `SAMRAIPointer` and bridge aliases where needed.
- Targeted legacy namespace bridges are centralized in the compatibility environment header.

### Tooling for Rollout

- Added per-file rewrite tool for non-legacy migration style:
  - `scripts/maintenance/rewrite_samrai_compatibility_file.py`
- This tool rewrites SAMRAI includes and type spellings to compatibility aliases and is intentionally heuristic (manual fix-up expected for ambiguous symbols).

### Process/Validation Progress

- The compatibility infrastructure has been committed and validated incrementally.
- Pilot updates (e.g. SnapshotCache) have been used to validate migration mechanics without mixing functional changes.

## What Happens Next

### Phase 1: Alias Adoption in Source Files

- Apply compatibility alias includes/types across IBTK/IBAMR in scoped batches.
- Prefer scripted rewrites first, then manual review/fixup for collisions and edge cases.
- Keep each batch buildable and focused on non-functional migration.

### Phase 2: Resolve Non One-to-One API Gaps

After alias adoption is broad, address true SAMRAI 2 -> 4 API differences that are not simple aliases (renames, removals, semantic shifts).

### Phase 3: Functional Validation

- Rebuild and run selected tests per subsystem during rollout.
- Run wider regression/testing once major migration batches are in.

## Working Rules

- Do not introduce umbrella headers that include all SAMRAI compatibility headers per translation unit.
- Keep migration commits separate from behavior changes.
- Prefer minimal per-file includes and explicit aliases.
- Treat scripted changes as first pass; always expect manual cleanup where names are ambiguous.
