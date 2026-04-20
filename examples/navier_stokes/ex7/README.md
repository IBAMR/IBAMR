This example demonstrates turbulence-statistics support for `INSStaggeredHierarchyIntegrator` in a live periodic-box Navier-Stokes run.

The solver advances an exact, divergence-free periodic velocity field driven by a matching body force. The field contains

1. a spatially varying steady mode, which produces a nontrivial mean velocity field, and
2. several time-periodic harmonic modes, which produce nonzero Reynolds stresses and turbulent kinetic energy.

The `PeriodicFlow` input database also supports seeded harmonic perturbations through

1. `noise_amplitude`
2. `noise_num_modes`
3. `noise_seed`
4. `max_spatial_harmonic`

so the example can be made richer without losing reproducibility.

Two statistics configurations are provided:

1. Running-mean mode (`input*.cell`, `input*.node`)
   The example advances for several full periods and compares the recovered mean velocity, Reynolds stresses, and TKE against the exact time average.
2. Phase-average mode (`input*.periodic.cell`, `input*.periodic.node`)
   The example stores one-period phase snapshots and compares them against the exact periodic solution. Since the forcing is deterministic, the phase-locked Reynolds stresses and TKE should vanish up to numerical error.

The periodic inputs are set up to discard an initial transient before the statistics are accumulated and to use multiple timesteps per stored phase snapshot. The most useful control knobs are

1. `NUM_PERIODS`, which increases the number of stored cycles,
2. `STATISTICS_START_TIME`, which delays accumulation until after a warm-up interval,
3. `NUM_SNAPSHOTS`, which controls the number of phase points per period,
4. `TIMESTEPS_PER_SNAPSHOT` and `DT_MAX`, which control the timestep size, and
5. `STEADY_STATE_THRESHOLD`, which controls the periodic-steady-state convergence check.

Unlike the branch tests, this example is solver-driven: it does not inject manufactured side data directly into the statistics object. Exact manufactured verification remains in `tests/navier_stokes/turbulence_statistics.cpp`.

Typical runs are:

```bash
./main2d input2d.cell
./main2d input2d.periodic.node
./main3d input3d.cell
./main3d input3d.periodic.node
```

The run reports diagnostic errors for the recovered statistics. Those errors should decrease with refinement and with smaller time steps.
