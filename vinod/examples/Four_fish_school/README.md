# Four Fish School Simulation - IBAMR Implementation

## Overview

This project implements a 2D fluid-structure interaction simulation of four fish swimming in a rectangular school formation using **IBAMR** (Immersed Boundary Adaptive Mesh Refinement). The simulation captures the complex hydrodynamic interactions between multiple self-propelled swimmers and the surrounding fluid.

## What is IBAMR?

**IBAMR** is an open-source software framework for solving fluid-structure interaction (FSI) problems using the immersed boundary (IB) method with adaptive mesh refinement (AMR). It enables high-resolution simulations of deformable and rigid bodies interacting with incompressible viscous fluids.

- **Website**: https://ibamr.github.io/
- **Documentation**: https://ibamr.github.io/docs
- **Repository**: https://github.com/IBAMR/IBAMR

### Key Features Used in This Project

- **Immersed Boundary Method**: Fluid-structure coupling without body-fitted meshes
- **Adaptive Mesh Refinement**: Automatic grid refinement near fish bodies
- **ConstraintIBMethod**: Kinematically prescribed body motions
- **INSStaggeredHierarchyIntegrator**: Incompressible Navier-Stokes solver
- **IBHydrodynamicForceEvaluator**: Computes forces and torques on swimmers

## Physics and Configuration

### Fish School Formation

The simulation models **4 fish in a rectangular formation**:

```
Fish-3 (top-left)      Fish-4 (top-right)
    x = 0.0, y = +0.2      x = 2.0, y = +0.2

Fish-1 (bottom-left)   Fish-2 (bottom-right)
    x = 0.0, y = -0.2      x = 2.0, y = -0.2
```

- **Axial spacing**: dx = 2.0L (head-to-head distance)
- **Lateral spacing**: dy = 0.4L (centerline-to-centerline)
- **Swimming mode**: In-phase undulation (all fish synchronized)

### Swimming Kinematics

Each fish undergoes traveling wave deformation based on the Bhalla et al. (2013) eel model:

**Body shape**: `h(s,t) = 0.125 * ((s + 0.03125)/1.03125) * sin(2πs - 6.28t)`

**Deformation velocity**: `-0.785 * ((s + 0.03125)/1.03125) * cos(2πs - 6.28t)`

Where:
- `s` is the arc length along the fish body (0 to 1)
- `t` is time
- Fish length L = 1.0

### Flow Parameters

- **Reynolds number**: Re = 5609 (based on fish length and swimming speed)
- **Fluid density**: ρ = 1.0
- **Dynamic viscosity**: μ = 0.785/Re ≈ 0.00014
- **Domain size**: 12.0 × 6.0 (length units)
- **Periodic boundary conditions**: All directions

### Computational Grid

- **Base grid**: 192 × 96 cells (N=64, domain 3N × 1.5N)
- **Refinement levels**: 3 levels with refinement ratio 4:1
- **Finest resolution**: Near fish bodies (automatically refined)
- **Time step**: Adaptive, CFL = 0.3, max Δt = 0.0001

## Project Structure

```
Four_fish_school/
├── example.cpp                    # Main simulation program (NS + IB + Odor)
├── IBEELKinematics.cpp            # Fish kinematics implementation
├── IBEELKinematics.h              # Kinematics class header
├── input2d                        # IBAMR input configuration
├── CMakeLists.txt                 # Build configuration
├── eel2d_1.vertex                 # Fish-1 geometry (bottom-left)
├── eel2d_2.vertex                 # Fish-2 geometry (bottom-right)
├── eel2d_3.vertex                 # Fish-3 geometry (top-left)
├── eel2d_4.vertex                 # Fish-4 geometry (top-right)
├── generate_4fish_vertices.py     # Script to generate fish positions
├── plot_eel_only.py               # Visualization: fish only
├── plot_combined_fluid_eel.py     # Visualization: fish + fluid
├── plot_odor_concentration.py     # Visualization: odor concentration field
├── analyze_odor_plumes.py         # Analysis: odor plume statistics
├── test_odor_transport_vortex_dynamics.py  # Validation script
├── README_ODOR_DYNAMICS.md        # Detailed odor transport documentation
└── *.pdf                          # Research papers (see References)
```

## Building the Simulation

### Prerequisites

1. **IBAMR** installed with dependencies:
   - PETSc (parallel linear algebra)
   - SAMRAI (structured AMR framework)
   - HDF5 (data output)
   - MPI (parallel computing)

2. **CMake** version ≥ 3.15

### Compilation Steps

```bash
# Create build directory
mkdir build
cd build

# Configure with CMake
cmake .. --debug-output

# Compile (using 4 parallel jobs)
make -j4

# Return to main directory
cd ..
```

This creates the executable `build/main2d`.

## Running the Simulation

### Basic Execution

```bash
# Run with 8 MPI processes
mpirun -np 8 ./build/main2d input2d
```

### Monitoring Performance

```bash
# Run with memory monitoring
mpirun -np 8 /usr/bin/time -v ./build/main2d input2d 2>&1 | grep -E "Maximum resident|Memory"

# Check available system memory
free -h
```

### Adjusting Resolution

If memory is limited, reduce the grid resolution in `input2d`:

```bash
# Reduce from N=64 to N=32 (8x less memory)
sed -i 's/^N = 64$/N = 32/' input2d
mpirun -np 4 ./build/main2d input2d
```

## Output Files

The simulation generates several output directories:

### Visualization Data

- **viz_eel2d_Str/**: VisIt/Silo visualization files
  - Written every 40 time steps
  - View with VisIt: `visit -o viz_eel2d_Str/dumps.visit`

### Restart Data

- **restart_IB2dStrDiv/**: Checkpoint files for restarting
  - Written every 150 time steps
  - Restart: `mpirun -np 8 ./build/main2d input2d restart_IB2dStrDiv <number>`

### Post-Processing Data

- **Eel2dStr/**: Hydrodynamic forces and fish kinematics
  - `Eel2d.curve`: Time series of forces, torques, velocities, COM positions
  - Updated every time step

### Log Files

- **IB2dEelStr.log**: Main simulation log
- **IB.log**: IBAMR initialization log

## Configuration Options

Key parameters in `input2d`:

### Physical Parameters

```
Re = 5609.0              # Reynolds number
MU = 0.785/Re            # Dynamic viscosity
RHO = 1.0                # Fluid density
```

### Solver Parameters

```
MAX_LEVELS = 3           # AMR levels
REF_RATIO = 4            # Refinement ratio between levels
N = 64                   # Base grid resolution
CFL_MAX = 0.3            # CFL number
DT_MAX = 0.0001          # Maximum time step
VORTICITY_TAGGING = TRUE # Refine based on vorticity
```

### Simulation Duration

```
START_TIME = 0.0
END_TIME = 30.0          # Simulation end time (≈ 5 tail beats)
```

### Output Intervals

```
viz_dump_interval = 40        # Visualization output
restart_dump_interval = 150   # Restart file output
timer_dump_interval = 100     # Performance timing output
```

## Fish Kinematics Configuration

Each fish has independent kinematics in the `ConstraintIBKinematics` section of `input2d`:

```
eel2d_1 {
    structure_names = "eel2d_1"
    structure_levels = MAX_LEVELS - 1
    calculate_translational_momentum = 1,1,0
    calculate_rotational_momentum = 0,0,1
    lag_position_update_method = "CONSTRAINT_POSITION"
    body_shape_equation = "0.125*((X_0+0.03125)/1.03125)*sin(2*PI*X_0-(0.785/0.125)*T)"
    deformation_velocity_function_0 = "(-0.785*((X_0+0.03125)/1.03125)*cos(2*PI*X_0-(0.785/0.125)*T))*N_0"
    deformation_velocity_function_1 = "(-0.785*((X_0+0.03125)/1.03125)*cos(2*PI*X_0-(0.785/0.125)*T))*N_1"
}
```

### Modifying Swimming Behavior

- **Phase offset**: Add phase to time term: `-(0.785/0.125)*T + phi`
- **Wave amplitude**: Change coefficient `0.125`
- **Wave speed**: Modify `0.785/0.125` (currently ≈ 6.28 rad/s)
- **Wavelength**: Change `2*PI` coefficient

## Hydrodynamic Force Analysis

The simulation computes for each fish:

- **Total force**: F = F_hydro + F_inertial
- **Torque**: About center of mass
- **Power**: P = F · v
- **Translational velocity**: COM velocity
- **Rotational velocity**: Angular velocity
- **Momentum**: Linear and angular

Data written to `Eel2dStr/Eel2d.curve`.

## Odor Transport Dynamics

This simulation includes **passive scalar transport** to model odor advection-diffusion coupled with fish swimming. The odor field C(x,y,t) satisfies:

```
∂C/∂t + u·∇C = κ∇²C + S(x,y,t)
```

### Key Features

- **Advection-Diffusion Solver**: IBAMR's `AdvDiffHierarchyIntegrator`
- **Coupling**: One-way (fluid affects odor, odor doesn't affect fluid)
- **Source Term**: Continuous Gaussian source at (-2.0, 0.0)
- **Schmidt Number**: Sc = ν/κ (configurable)
- **Visualization**: Odor contours, vortex-odor interaction

### Physical Parameters

```
KAPPA = 1.0e-3              # Molecular diffusivity
SCHMIDT = MU / (RHO * KAPPA)  # Schmidt number
```

### Odor Visualization

```bash
# Visualize odor concentration at iteration 200
python plot_odor_concentration.py 200

# Analyze odor plume statistics
python analyze_odor_plumes.py --stats

# Single iteration analysis
python analyze_odor_plumes.py 200
```

### Output Fields

- **C**: Odor concentration (normalized C* = (C-C_l)/(C_h-C_l))
- **∇C**: Concentration gradient (odor sharpness)
- **Mixing efficiency**: Vortex-enhanced spreading metric

**See [README_ODOR_DYNAMICS.md](README_ODOR_DYNAMICS.md) for complete documentation.**

## Visualization

### Using Python Scripts

```bash
# Plot fish positions only
python plot_eel_only.py

# Plot fish with fluid velocity/vorticity
python plot_combined_fluid_eel.py
```

### Using VisIt

```bash
visit -o viz_eel2d_Str/dumps.visit
```

Recommended pseudocolor plots:
- **Velocity magnitude**: Shows flow around fish
- **Vorticity**: Shows vortex structures
- **Pressure**: Shows pressure distribution

## Scientific Background

### References

#### Fish Swimming and Hydrodynamics

1. **Bhalla et al. (2013)**: "A unified mathematical framework and an adaptive numerical method for fluid-structure interaction with rigid, deforming, and elastic bodies." *Journal of Computational Physics*, 250:446-476.
   - Eel kinematics model used in this simulation

2. **Peskin (2002)**: "The immersed boundary method." *Acta Numerica*, 11:479-517.
   - Foundational immersed boundary method

3. **Griffith & Patankar (2020)**: "Immersed Methods for Fluid-Structure Interaction." *Annual Review of Fluid Mechanics*, 52:421-448.
   - Review of IB methods in FSI

#### Odor Transport and Vortex Dynamics

4. **"How does vortex dynamics help undulating bodies spread odor.pdf"** (included in repository)
   - Vortex-enhanced odor spreading mechanisms
   - Flapping kinematics effects on mixing
   - Quantification of mixing efficiency
   - *Directly applicable to this simulation's odor transport module*

5. **"Navigation in odor plumes How do the flapping kinematics modulate the odor landscape.pdf"** (included in repository)
   - Odor landscape structure in fish wakes
   - Effects of schooling formations on odor plumes
   - Sensory information content in turbulent odor fields
   - *Provides context for analyzing simulation results*

### Fish Schooling Hydrodynamics

This simulation enables study of:
- **Wake interactions**: How downstream fish interact with leader vortices
- **Energy savings**: Reduced swimming cost in formations
- **Lateral forces**: Side-to-side interactions between neighbors
- **Collective dynamics**: Emergent behaviors from local interactions

## Performance Considerations

### Memory Requirements

Approximate memory usage:
- **N=32**: ~2-4 GB (4 processes)
- **N=64**: ~16-32 GB (8 processes)
- **N=128**: ~128-256 GB (16+ processes)

### Computational Cost

- **N=64**: ~10-20 hours for END_TIME=30.0 (8 cores)
- Scales with: grid resolution, number of fish, refinement levels

### Optimization Tips

1. Use appropriate MPI process count (typically 4-16 for N=64)
2. Enable vorticity tagging for efficient AMR
3. Adjust `viz_dump_interval` (higher = less I/O)
4. Use restart files for long simulations

## Troubleshooting

### Common Issues

**"Out of memory" error**:
- Reduce `N` in input2d (N=32 or N=48)
- Reduce `MAX_LEVELS` to 2
- Increase MPI processes to distribute memory

**Simulation crashes**:
- Check `DT_MAX` is not too large (≤ 0.0001)
- Verify `CFL_MAX ≤ 0.3`
- Ensure vertex files match input2d configuration

**Slow performance**:
- Increase MPI process count
- Reduce output frequency
- Check load balancing in log files

## License

This code is distributed under the 3-clause BSD license, consistent with IBAMR.

Copyright (c) 2014-2024 by the IBAMR developers. All rights reserved.

## Contact and Support

For IBAMR-specific questions:
- IBAMR GitHub Issues: https://github.com/IBAMR/IBAMR/issues
- IBAMR Google Group: https://groups.google.com/g/ibamr-users

## Acknowledgments

This implementation follows the official IBAMR eel2d example and extends it to multi-fish schooling configurations. Thanks to the IBAMR development team for providing this powerful FSI framework.
