This example set is for turbulence-oriented Navier-Stokes examples.

It reuses the generic Navier-Stokes driver from `ex0` and currently provides 3D benchmark-style setups:

1. `input3d.taylor_green_les`
   a periodic 3D Taylor-Green vortex with Smagorinsky LES enabled; this is a standard transition-to-turbulence benchmark and is a good first case for visualizing coherent structure breakdown
2. `input3d.decaying_hit`
   a periodic 3D decaying homogeneous-isotropic-turbulence-style case initialized from a smooth divergence-free multimode field and evolved with Smagorinsky LES
3. `input3d.channel_flow_les`
   a wall-bounded channel-flow LES scaffold with periodic streamwise/spanwise directions, no-slip walls, and constant streamwise forcing; this is the right starting point for mean-profile and Reynolds-stress studies, even though the default run is intentionally much lighter than a reference channel benchmark
4. `input3d.channel_flow_retau180_spinup`
   the spin-up phase for a more serious channel-flow verification case targeting approximately `Re_tau = 180`
5. `input3d.channel_flow_retau180_average`
   the averaging phase for that same case, intended to be run from restart after the spin-up phase so that turbulence statistics are accumulated only after transients have decayed

The defaults are chosen to be illustrative while still remaining light enough to run as examples:

1. `input3d.taylor_green_les` runs to `t = 0.4` on a `32 x 32 x 32` periodic grid
2. `input3d.decaying_hit` runs to `t = 0.5` on a `32 x 32 x 32` periodic grid
3. `input3d.channel_flow_les` runs to `t = 0.6` on a `48 x 32 x 32` grid with periodic `x`/`z` and wall-normal no-slip boundaries
4. the `Re_tau = 180` benchmark pair uses a larger `128 x 96 x 96` grid and a two-stage workflow instead of a single short run

Suggested visualization fields are:

1. velocity magnitude
2. vorticity
3. Q-criterion or vortex-magnitude surrogates when available
4. pressure

Typical runs are:

```bash
./main3d input3d.taylor_green_les
./main3d input3d.decaying_hit
./main3d input3d.channel_flow_les
```

The `Re_tau = 180` benchmark workflow is:

```bash
./main3d input3d.channel_flow_retau180_spinup
./main3d input3d.channel_flow_retau180_average restart_channel_flow_retau180_spinup <restart step>
```

The intended interpretation is:

1. `input3d.channel_flow_retau180_spinup` drives the flow to a turbulent state and writes restart files
2. `input3d.channel_flow_retau180_average` restarts from one of those files and accumulates turbulence statistics over a longer averaging window

This is the closest thing in `ex9` to a serious DNS/LES verification case. It is still limited by the current generic driver:

1. it uses constant body-force driving rather than a constant-mass-flux controller
2. the statistics are the current mean and second-moment turbulence statistics already implemented in IBAMR

The new `statistics_start_time` option in `TurbulenceStatistics` can also be used to delay accumulation in a single continuous run. In `input3d.channel_flow_retau180_average` it is set to `20.0`, so statistics only begin at that simulation time.

All cases write VisIt output by default.
