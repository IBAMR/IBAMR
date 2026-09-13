These tests cover the main phase-change capabilities and the defects repaired
with them. They use representative configurations rather than an exhaustive
matrix of input options. Solver, hierarchy, and boundary settings remain in
ordinary input files.

- The original Stefan cases retain their physical mass references for both
  Enthalpy and Allen-Cahn, including the existing MPI aliases.
- Allen-Cahn `amr.restart=2` and Enthalpy `open_boundary.restart=2` exercise
  hierarchy changes and restart through subsequent timesteps. The latter
  includes nonzero flow and physical boundaries in 2D and 3D. These cases also
  check the stored fraction gradient and divergence source.
- The 2D Enthalpy `source_transfer` case assigns a nonzero divergence source,
  moves the refined region, and checks the transferred values. One case covers
  this shared base-class operation.
- The 2D `extrapolate0.amr.restart=2` case extends the PCM liquid fraction into
  gas through regridding and restart, including the extrapolated field in the
  trajectory comparison. The 3D `extrapolate1` case approaches the PCM value
  from the other side and exercises the separate dimensional implementation.
  Both check PCM preservation, temperature, and extension into adjacent gas.
- The two `laser` cases combine laser forcing, phase-dependent material
  properties, and tagging. Laser power, location, and smoothing width are
  checked on live data. Distinct gas, solid, and liquid coefficients make
  incorrect material mixing observable; tagging assertions run in the same
  initialized hierarchy.

The separate Allen-Cahn profile sweep was removed: those formulas largely
predate this PR, and representative profiles are exercised by the retained
physical cases. The surface-tension, thermocapillary, regrid-projection, and
Fortran gradient regressions remain under `tests/multiphase_flow`.

All checks use actual SAMRAI patch data on initialized IBAMR hierarchies.
Numerical results are reported through `plog` to `output` for `numdiff`. The
restart cases compare each patch's cell and side values, including separate
copies of shared sides, against the uninterrupted trajectory. Differences are
scaled by `max(1, abs(reference))`; these checks do not rely on aggregate field
norms alone.

Run `attest -R phase_change` from the Debug build root. Use `attest -N` with the
same selection to inspect discovery.
