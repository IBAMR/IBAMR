The Stefan tests also provide focused checks for the phase-change library APIs:

- `amr.restart=2` moves the refined region on each timestep, checks every stored liquid-fraction gradient component against the live fraction field, and compares every cell and side value after restart with the uninterrupted trajectory. The divergence source is included in that comparison.
- `open_boundary.restart=2` adds nonzero flow and physical open boundaries to the moving-mesh restart check. New fine patches reach those boundaries before the first restored timestep, so the test detects boundary wrappers that still use constructor defaults.
- `source_transfer` assigns a constant divergence source, moves the refined region, and checks all old and newly created patch cells. This rejects the former scratch-only registration, which left the source undefined on new patches.
- `tags` checks inclusive liquid-fraction thresholds at initialization, positive and negative gradient components later, and preservation of existing tags.
- `extrapolate0` and `extrapolate1` check extension from constant PCM liquid fraction into gas initialized on either side of that value. The PCM value, fraction bounds, and constant temperature must be preserved.
- `profile_*` checks all seven Allen-Cahn interpolation profiles through the integrator's source evaluation. Hybrid transition points are computed independently from their matching equation; fractions on both sides and temperatures below, at, and above melting are checked.
- `laser` checks initial zero forcing, integrated power, motion of the source centroid, callback scaling, two smoothing kernels, and three time discretizations.

All fixtures use actual SAMRAI patch data on initialized IBAMR hierarchies. The small cases supplement the original Stefan reference outputs. Run them from a Debug build with `attest -R phase_change` and the configured MPI and numdiff executables.

Numerical checks report compact values or error norms through `plog` for `numdiff`.
The restart cases retain per-rank uninterrupted field data only to compare the
same cell and side values after restart; the reported maximum scaled difference
uses `max(1, abs(reference))` as its scale. This detects local errors that separate
field norms can hide. Integer tagging and mesh-change checks use direct errors.
Verify discovery with `attest -N -R phase_change` from the configured build root.
