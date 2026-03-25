This example verifies the Smagorinsky SGS implementation for `INSStaggeredHierarchyIntegrator` with manufactured velocity fields.

The example does not advance the Navier-Stokes equations. Instead, it:

1. builds a hierarchy and a staggered INS integrator with `TurbulenceModel { model_type = "SMAGORINSKY" }`
2. prescribes an exact velocity field on the staggered grid
3. verifies the native-location strain-rate reconstruction used by the SGS machinery
4. evaluates the SGS force and compares it to an independently assembled exact discrete force

The manufactured cases are:

1. `UNIFORM`
   checks that the reconstructed strain and SGS force are zero
2. `RIGID_ROTATION`
   checks that the symmetric strain vanishes even though the velocity gradient is nonzero
3. `SIMPLE_SHEAR`
   checks nonzero native shear strain with zero SGS-force divergence
4. `TRIGONOMETRIC`
   checks a fully space-varying case with nontrivial cell, node/edge, and side-centered SGS data

The verification is done at the same locations used by the implementation:

1. cell centers for the diagonal strain and normal stresses
2. nodes in 2D or edges in 3D for the shear strain and shear stresses
3. side centers for the SGS force

Typical runs are:

```bash
./main2d input2d.uniform
./main2d input2d.trigonometric
./main3d input3d.simple_shear
./main3d input3d.trigonometric
```

A successful run reports the maximum cell-strain, native-shear-strain, and SGS-force errors and ends with `verification = PASS`.
