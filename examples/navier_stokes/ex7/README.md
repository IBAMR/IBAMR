This example demonstrates turbulence-statistics support for `INSStaggeredHierarchyIntegrator` using a manufactured velocity field instead of a full Navier-Stokes solve.

The main idea is simple: we prescribe a staggered velocity field, feed that field into the turbulence-statistics object, and compare the recovered statistics against exact values. This lets us verify the statistics machinery in isolation.

Two kinds of quantities are checked:

1. The mean velocity, `<U>`
2. The second-order statistics derived from `<UU>`, namely the Reynolds stresses
   `R_ij = <u_i u_j> - <u_i><u_j>`
   and the turbulent kinetic energy
   `k = 0.5 R_ii`

The example supports two analysis centerings:

1. `CELL`: statistics are accumulated and verified on cell centers
2. `NODE`: statistics are accumulated and verified on nodes

This is useful because solver data naturally lives on sides in the staggered formulation, while postprocessing and visualization are often easier on cells or nodes.

The repository includes two sampling modes:

1. Running-mean mode (`input*.cell`, `input*.node`)
   The prescribed field has a steady mean plus a periodic fluctuation. The example verifies that the accumulated mean, Reynolds stresses, and TKE match the exact time-averaged values.
2. Periodic-phase mode (`input*.periodic.cell`, `input*.periodic.node`)
   The statistics object stores snapshots over one period and updates them over repeated cycles. The example verifies that the recovered phase snapshots match the exact periodic solution and that the phase-locked Reynolds stresses and TKE are zero.

In all cases, the manufactured field is spatially uniform. That keeps the verification focused on the statistics logic itself rather than on interpolation error between side-, cell-, and node-centered data.

Typical runs are:

```bash
./main2d input2d.cell
./main2d input2d.periodic.node
./main3d input3d.cell
./main3d input3d.periodic.node
```

A successful run reports that the manufactured turbulence-statistics verification passed, with errors at roundoff level.
