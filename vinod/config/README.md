# Configuration Files

This directory contains parameter files and input configurations for IBAMR simulations.

## Files

### default_params.input
Comprehensive template with all common IBAMR parameters documented.

**Usage**:
```bash
cp config/default_params.input my_simulation.input
# Edit my_simulation.input for your specific case
./main2d my_simulation.input
```

**Key Sections**:
1. **Physical Parameters** - Fluid properties (viscosity, density)
2. **Grid Parameters** - Domain size and resolution
3. **Temporal Parameters** - Time stepping configuration
4. **AMR Parameters** - Adaptive mesh refinement settings
5. **IB Method Parameters** - Immersed boundary configuration
6. **Solver Parameters** - Linear/nonlinear solver options
7. **Output Parameters** - Visualization and data output
8. **Initial/Boundary Conditions** - Problem setup

## Parameter Categories

### Critical Parameters for Stability

**CFL Condition**:
```
CFL_MAX = 0.3              // Decrease if simulation is unstable
DT_MAX = 0.01              // Maximum timestep
```

**IB Kernel**:
```
IB_DELTA_FUNCTION = "IB_4" // Options: IB_2, IB_3, IB_4, IB_6
                           // Higher order = smoother but more expensive
```

**Grid Resolution**:
```
N = 128                    // Increase for better accuracy
IB_POINT_DENSITY = 2.0     // Lagrangian points per grid spacing
```

### Performance Parameters

**AMR Settings**:
```
MAX_LEVELS = 2             // More levels = finer resolution near IB
REGRID_INTERVAL = 4        // More frequent = more accurate but slower
```

**Parallel Performance**:
```
LoadBalancer {
   bin_pack_method = "SPATIAL"  // or "GREEDY"
   max_workload_factor = 1      // Load imbalance tolerance
}
```

### Output Control

**Visualization Frequency**:
```
VIZ_DUMP_INTERVAL = 100    // Write every N timesteps
                           // Smaller = more output, larger files
```

**Data to Output**:
```
INSStaggeredHierarchyIntegrator {
   output_U = TRUE         // Velocity field
   output_P = TRUE         // Pressure field
   output_F = FALSE        // IB force density
   output_Omega = TRUE     // Vorticity
   output_Div_U = FALSE    // Velocity divergence
}
```

## Common Configurations

### High Reynolds Number Flow
```
MU = 0.001                 // Lower viscosity
CFL_MAX = 0.2              // More restrictive for stability
DT_MAX = 0.005             // Smaller timestep
N = 256                    // Finer grid
```

### Low Reynolds Number / Stokes Flow
```
MU = 0.1                   // Higher viscosity
CFL_MAX = 0.5              // Can be less restrictive
DT_MAX = 0.02              // Larger timestep allowed
N = 64                     // Coarser grid may suffice
```

### 3D Simulation
```
NDIM = 3                   // Set at compile time
N = 64                     // Typically need coarser due to cost
MAX_LEVELS = 2             // Limit levels for 3D
```

### AMR-Heavy Problem
```
MAX_LEVELS = 3             // More refinement levels
REGRID_INTERVAL = 2        // Frequent regridding
USE_VORTICITY_TAGGING = TRUE
VORTICITY_REL_THRESH = 0.1 // More aggressive tagging
```

## Boundary Condition Examples

### Uniform Inflow (x-direction)
```
VelocityBcCoefs_0 {
   gcoef_function_0 = "1.0"    // x_lo: u = 1.0
   gcoef_function_1 = "0.0"    // x_hi: u = 0.0
}
```

### No-Slip Walls
```
VelocityBcCoefs_0 {
   acoef_function_2 = "1.0"    // y_lo: Dirichlet
   bcoef_function_2 = "0.0"
   gcoef_function_2 = "0.0"    // u = 0
}
VelocityBcCoefs_1 {
   acoef_function_2 = "1.0"
   bcoef_function_2 = "0.0"
   gcoef_function_2 = "0.0"    // v = 0
}
```

### Periodic Boundaries
```
CartesianGeometry {
   periodic_dimension = 1, 1   // Periodic in both x and y
}
```

### Outflow (Neumann)
```
VelocityBcCoefs_0 {
   acoef_function_1 = "0.0"    // x_hi: Neumann
   bcoef_function_1 = "1.0"
   gcoef_function_1 = "0.0"    // du/dn = 0
}
```

## Solver Options

### PETSc Krylov Solver
```
SOLVER_TYPE = "PETSC_KRYLOV_SOLVER"
MAX_ITERATIONS = 100
ABS_RESIDUAL_TOL = 1.0e-8
REL_RESIDUAL_TOL = 1.0e-2
```

### Hypre Preconditioner
```
PC_TYPE = "SHELL"
PC_LEVEL_TYPE = "HYPRE_LEVEL_SOLVER"
PC_HYPRE_TYPE = "PFMG"         // or "SMG", "BoomerAMG"
```

## Tips

1. **Start Simple**: Begin with coarse grid, low resolution
2. **Convergence Study**: Double resolution and check if results change
3. **Monitor Logs**: Watch for CFL violations, solver failures
4. **Visualization**: Check output frequently during development
5. **Restart Files**: Enable for long runs in case of crashes

## Troubleshooting

| Problem | Likely Cause | Solution |
|---------|-------------|----------|
| Simulation crashes | CFL violation | Reduce `DT_MAX` or `CFL_MAX` |
| Poor IB resolution | Grid too coarse | Increase `N` or `IB_POINT_DENSITY` |
| Solver not converging | Bad preconditioner | Try different `PC_HYPRE_TYPE` |
| Slow performance | Too much output | Increase `VIZ_DUMP_INTERVAL` |
| Load imbalance | Poor partitioning | Adjust `LoadBalancer` settings |

## References

- IBAMR Input Database Documentation: http://ibamr.github.io
- SAMRAI User Guide: https://computing.llnl.gov/projects/samrai
- PETSc Documentation: https://petsc.org/release/
