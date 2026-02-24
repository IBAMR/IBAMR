# SAMRAI Compatibility Migration Plan

For a short adoption guide (with a simple before/after example), see
`doc/samrai_compatibility_adoption.md`.

## Purpose

Track the migration from IBSAMRAI2-era SAMRAI usage to the SAMRAI compatibility layer, and capture what is complete versus what is still pending.

## Current State

### Compatibility Infrastructure

- Canonical alias headers exist under `ibtk/include/samrai_compatibility/<module>/SAMRAI*.h`.
- Legacy wrapper coverage exists under `ibtk/include/samrai_compatibility/legacy/...`.
- Common entry headers are in:
  - `ibtk/include/ibtk/samrai_compatibility_names.h`
  - `ibtk/include/ibtk/samrai_compatibility_layer.h`
  - `ibtk/include/ibtk/samrai_compatibility_legacy_aliases.h`
- Compatibility include/install paths are exported by CMake.

### Adoption Status

One-to-one alias adoption is complete for the main code regions:

- IBTK library sources and headers.
- IBAMR library sources and headers.
- IBTK and IBAMR examples.
- Testsuite sources.

### Tooling

- Rewrite tool: `scripts/maintenance/rewrite_samrai_compatibility_file.py`.
- Legacy wrapper generator: `scripts/maintenance/generate_samrai_legacy_wrappers.py`.

## Remaining Work

1. Handle non one-to-one SAMRAI API differences between IBSAMRAI2 and SAMRAI 4.5.0.
2. Resolve any build/runtime issues that appear when switching builds from IBSAMRAI2 to SAMRAI 4.5.0.
3. Remove temporary migration-specific workarounds once SAMRAI 4 behavior is fully validated.

## Active Rules

- Keep migration changes behavioral-neutral unless a SAMRAI API difference requires functional updates.
- Keep compatibility alias headers centralized; avoid introducing ad hoc alias definitions in unrelated files.
- Prefer minimal, explicit includes over broad umbrella inclusion in ordinary translation units.
- Use script + manual fixup workflow for large edits; validate with incremental builds.
