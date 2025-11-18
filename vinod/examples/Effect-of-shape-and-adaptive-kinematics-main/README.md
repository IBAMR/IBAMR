# Effects of Reynolds Number and Thickness on Undulatory Self-Propelled Foil

This repository contains an IBAMR-based simulation of an undulatory self-propelled foil (eel-like swimmer) with **adaptive kinematics** based on Reynolds number and foil thickness. The implementation studies the hydrodynamic performance across different flow regimes and geometric configurations.

## ⚠️ Important: Two Separate Implementations

**This repository contains the ADAPTIVE kinematics research code** for exploring Reynolds number and thickness effects.

For **exact reproduction of Gupta et al. (2022)** with fixed parameters (Re=5000, St=0.6), see the separate validation suite:

📁 **Validation Suite Location**: `/home/user/Validation-Gupta2022/`

| This Repository (Adaptive) | Validation Suite (Separate) |
|---------------------------|----------------------------|
| Variable Re (1000-10000+) | Fixed Re = 5000 |
| Adaptive kinematics | Exact Eqs. 3-6 from paper |
| Research & parameter studies | Benchmarking & verification |
| Multiple profiles evolving | 5 specific NACA profiles |

See `QUICKSTART.md` in the validation suite for benchmarking instructions.

## Overview

This code is **inspired by and extends** the research described in:

**Primary Reference**:
- **Zhang, D., Pan, G., Chao, L., & Zhang, Y. (2018)**. "Effects of Reynolds number and thickness on an undulatory self-propelled foil." *Physics of Fluids*, 30(7), 071902. DOI: [10.1063/1.5034439](https://doi.org/10.1063/1.5034439)

**Secondary Reference** (see separate validation suite):
- **Gupta, S., Sharma, A., Agrawal, A., Thompson, M. C., & Hourigan, K. (2022)**. "Anguilliform and carangiform fish-inspired hydrodynamic study for an undulating hydrofoil: Effect of shape and adaptive kinematics"

**Extensions in This Code** (beyond Zhang et al. 2018):
- Reynolds-dependent adaptive amplitude and frequency modulation
- Envelope power adaptation for anguilliform vs. carangiform modes
- Extended parameter study framework

### Key Features

1. **Reynolds Number Effects**: Simulates foil propulsion across different Reynolds numbers (Re = 1000 to 10000+)
2. **Thickness Variations**: Studies the effect of foil thickness ratio (h/L = 0.02 to 0.08+)
3. **Adaptive Kinematics**: Automatically adjusts swimming parameters (amplitude, frequency) based on Re and thickness
4. **Swimming Modes**: Supports anguilliform, carangiform, and intermediate swimming modes
5. **Performance Metrics**: Tracks thrust, power, swimming speed, and propulsive efficiency
6. **Shape Adaptation**: Models both slender (anguilliform) and thicker (carangiform) body shapes

## Physics and Implementation

### Adaptive Kinematics Model

The code implements Reynolds-dependent adaptive kinematics:

#### Amplitude Adaptation
```
A_adapted = A_base × f_Re(Re) × f_h(h/L)
```
- At **low Re** (Re < 5000): Amplitude increases due to viscous effects
- At **high Re** (Re > 5000): Amplitude decreases for better efficiency
- **Thicker foils** require larger amplitude to generate sufficient thrust

#### Frequency Adaptation
```
f_adapted = f_base × g_Re(Re) × g_h(h/L)
```
- At **low Re**: Frequency decreases (viscous damping)
- At **high Re**: Frequency can increase (inertial regime)
- **Thicker foils** undulate at lower frequencies

#### Envelope Power
The amplitude distribution along the body is controlled by:
```
A(s) = A_max × (s/L)^p
```
- **Anguilliform** (p ≈ 1): Linear amplitude increase
- **Carangiform** (p ≈ 2-3): Amplitude concentrated at tail

### Swimming Modes

| Mode | swimming_mode | Characteristics |
|------|---------------|-----------------|
| **Anguilliform** | 0.0 - 0.3 | Whole-body undulation, thin foil |
| **Mixed** | 0.3 - 0.7 | Intermediate behavior |
| **Carangiform** | 0.7 - 1.0 | Posterior body undulation, thicker foil |

## File Structure

```
├── IBEELKinematics.h              # Header file with new adaptive features
├── IBEELKinematics.cpp            # Implementation of adaptive kinematics
├── example.cpp                    # Main simulation driver
├── input2d                        # Original baseline configuration
├── input2d_Re1000_h004           # Low Re, standard thickness (anguilliform)
├── input2d_Re5609_h006           # Baseline Re, intermediate thickness (mixed)
├── input2d_Re10000_h008          # High Re, thick foil (carangiform)
├── eel2d.vertex                   # Lagrangian mesh definition
├── eel2d_straightswimmer.m       # MATLAB mesh generator
├── pycodeforvetexshift.py        # Vertex position adjustment tool
├── analyze_performance.py        # Python analysis script for results
└── CMakeLists.txt                # Build configuration
```

## Compilation

### Prerequisites
- IBAMR library (version 0.10.0 or later)
- PETSc, SAMRAI, and dependencies
- CMake 3.15 or later
- MPI (OpenMPI or MPICH)
- Python 3.x (for analysis)
- NumPy, Matplotlib (for analysis scripts)

### Build Instructions

```bash
# Create build directory
mkdir build
cd build

# Configure with CMake
cmake .. --debug-output

# Compile
make -j4

# Return to main directory
cd ..
```

## Running Simulations

### Single Simulation

Run a simulation with a specific configuration:

```bash
mpirun -np 6 ./build/main2d input2d_Re1000_h004
```

### Parameter Study

Run multiple simulations to study Reynolds and thickness effects:

```bash
# Low Reynolds number (viscous regime)
mpirun -np 6 ./build/main2d input2d_Re1000_h004

# Baseline Reynolds number
mpirun -np 6 ./build/main2d input2d

# Intermediate case
mpirun -np 6 ./build/main2d input2d_Re5609_h006

# High Reynolds number (inertial regime)
mpirun -np 6 ./build/main2d input2d_Re10000_h008
```

## Configuration Parameters

Key parameters in the input files:

### Reynolds Number and Thickness
```
reynolds_number    = 5609.0      # Reynolds number (UL/ν)
thickness_ratio    = 0.04        # Foil thickness ratio (h/L)
```

### Swimming Mode
```
swimming_mode      = 0.0         # 0 = anguilliform, 1 = carangiform
```

### Base Kinematics
```
base_amplitude     = 0.125       # Base tail amplitude
base_frequency     = 0.785       # Base undulation frequency
```

### Adaptive Features
```
enable_shape_adaptation = TRUE   # Enable Re-dependent adaptation
envelope_power         = 1.0     # Amplitude envelope exponent
tail_width_ratio       = 0.02    # Tail thickness (for carangiform)
```

### Performance Tracking
```
track_performance      = TRUE                    # Enable metrics logging
performance_log_file   = "performance_Re5609.dat"  # Output file
```

## Output and Analysis

### Simulation Outputs

Each simulation produces:

1. **Visualization Files** (in `viz_*/` directories)
   - Velocity fields, pressure, vorticity
   - Can be viewed with VisIt or ParaView

2. **Performance Metrics** (`performance_*.dat`)
   - Time-resolved thrust, power, speed, efficiency
   - Adapted amplitude and frequency
   - Strouhal number

3. **Structure Diagnostics** (in `Results_*/` directories)
   - Center of mass position and velocity
   - Drag, torque, power
   - Momentum components

### Analysis Tools

Analyze results using the provided Python script:

```bash
# Analyze all performance files in directory
python analyze_performance.py

# Analyze specific files
python analyze_performance.py performance_Re1000_h004.dat performance_Re10000_h008.dat
```

The script generates:
- **Time series plots**: Speed, thrust, power, efficiency vs. time
- **Kinematics plots**: Adapted amplitude and frequency evolution
- **Comparative plots**: Performance across different Re and h/L
- **Summary table**: Average performance metrics

### Key Performance Metrics

- **Swimming Speed**: Time-averaged forward velocity
- **Propulsive Efficiency**: η = (Thrust × Speed) / Power
- **Strouhal Number**: St = fA/U (optimal ≈ 0.2-0.4 for swimming)
- **Cost of Transport**: Power / (Speed × Weight)

## Expected Results

### Reynolds Number Effects

| Reynolds | Regime | Characteristics |
|----------|--------|-----------------|
| Re < 1000 | Viscous | Low speed, high power, low efficiency |
| Re ≈ 5000 | Transitional | Moderate performance, optimal for many fish |
| Re > 10000 | Inertial | Higher speed possible, better efficiency |

### Thickness Effects

| h/L | Swimming Mode | Characteristics |
|-----|---------------|-----------------|
| 0.02-0.04 | Anguilliform | High frequency, whole-body undulation |
| 0.04-0.06 | Mixed | Intermediate behavior |
| 0.06-0.08+ | Carangiform | Lower frequency, tail-focused propulsion |

### Adaptive Behavior

The code automatically adjusts kinematics based on Re and h/L:
- **Low Re + Thin foil**: Increases amplitude, decreases frequency
- **High Re + Thick foil**: Optimizes for efficiency, concentrates motion at tail
- **Intermediate cases**: Balances thrust generation and efficiency

## Physical Validation

The implementation includes several biologically-inspired features:

1. **Strouhal Number**: Adaptive kinematics target St ≈ 0.2-0.4, consistent with efficient swimmers
2. **Amplitude Envelope**: Matches observed fish kinematics (linear for anguilliform, polynomial for carangiform)
3. **Reynolds Scaling**: Follows empirical scaling laws from fish swimming studies
4. **Efficiency Trends**: Reproduces experimental observations of Re and thickness effects

## References

### Primary References

1. **Zhang, D., Pan, G., Chao, L., & Zhang, Y. (2018).** "Effects of Reynolds number and thickness on an undulatory self-propelled foil." *Physics of Fluids*, 30(7), 071902. https://doi.org/10.1063/1.5034439
   - **PDF Location**: `Effects of Reynolds number and thickness on an undulatory self propelled foil.pdf`
   - **GitHub**: https://github.com/vinodthale/Effect-of-shape-and-adaptive-kinematics

2. **Gupta, S., Sharma, A., Agrawal, A., Thompson, M. C., & Hourigan, K. (2022).** "Anguilliform and carangiform fish-inspired hydrodynamic study for an undulating hydrofoil: Effect of shape and adaptive kinematics"
   - **PDF Location**: `Anguilliform and carangiform fish-inspired hydrodynamic study for an undulating hydrofoil_ Effect of shape and adaptive kinematics.pdf`
   - **Validation Suite**: `/home/user/Validation-Gupta2022/`

### Supporting References

3. IBAMR: An adaptive and distributed-memory parallel implementation of the immersed boundary method
4. Lauder, G. V., & Tytell, E. D. (2006). Hydrodynamics of undulatory propulsion

### Parameter Verification

For detailed comparison of implementation parameters vs. Zhang et al. (2018), see:
- `/home/user/PARAMETER_VERIFICATION_Zhang2018.md`

## Troubleshooting

### Common Issues

1. **Simulation crashes or becomes unstable**
   - Reduce `DT_MAX` in input file
   - Increase grid resolution (`N`)
   - Check that Re and kinematics parameters are reasonable

2. **Performance metrics file is empty**
   - Ensure `track_performance = TRUE`
   - Check file permissions in output directory
   - Verify simulation ran long enough (t > 1.0)

3. **Compilation errors**
   - Verify IBAMR installation and environment variables
   - Check that all dependencies are properly linked
   - Ensure C++11 or later compiler

4. **Analysis script fails**
   - Install required Python packages: `pip install numpy matplotlib`
   - Check that performance files exist and are valid
   - Verify file format matches expected columns

## Contributing

To add new swimming modes or adaptive strategies:

1. Modify `calculateAdaptiveKinematics()` in `IBEELKinematics.cpp`
2. Add new parameters to input files
3. Update `IBEELKinematics.h` with new member variables
4. Test with parameter study

## License

This code is based on IBAMR and follows the 3-clause BSD license.
See COPYRIGHT file in IBAMR distribution for details.

## Contact

For questions about this implementation:
- Review the included research papers
- Check IBAMR documentation: https://ibamr.github.io
- Post issues on the repository

---

**Last Updated**: 2025
**IBAMR Version**: 0.10.0+
**Status**: Research Implementation
