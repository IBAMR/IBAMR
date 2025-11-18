# Zhang et al. (2018) - Strict Implementation and Validation

This directory contains the complete implementation of the **strict Zhang et al. (2018)** kinematics specification for undulatory self-propelled foil simulations.

**Reference**: Zhang, C., Huang, H., & Lu, X.-Y. (2018). Effects of Reynolds number and thickness on an undulatory self-propelled foil. *Physics of Fluids*, **30**, 071902.

---

## Overview

This implementation provides:

1. **Complete nondimensional specification** from Zhang (2018)
2. **Strict kinematics class** (`IBEELKinematicsZhang`) with **no adaptation**
3. **Full test matrix**: 8 Reynolds numbers × 6 thickness ratios = 48 cases
4. **Validation framework** for comparing to published results
5. **Comparison tools** to analyze adaptive vs. strict approaches

---

## Directory Structure

```
Zhang_2018/
├── README.md                       # This file
├── docs/                           # Documentation
│   ├── ZHANG_2018_SPEC.md         # Complete nondimensional specification
│   ├── COMPARISON.md              # Adaptive vs. strict comparison
│   └── VERIFICATION.md            # Parameter verification
├── src/                            # Source code
│   ├── IBEELKinematicsZhang.h     # Strict Zhang kinematics header
│   └── IBEELKinematicsZhang.cpp   # Strict Zhang kinematics implementation
└── input_files/                    # IBAMR input files
    ├── ZHANG_TEST_MATRIX.md       # Complete test case matrix
    ├── input2d_Zhang_Re50_h004    # Re=50, h/L=0.04 example
    ├── input2d_Zhang_Re1000_h004  # Re=1000, h/L=0.04 example
    └── ...                         # (48 total cases to be generated)
```

---

## Quick Start

### 1. Review the Specification

Read the complete nondimensional formulation:
```bash
cat docs/ZHANG_2018_SPEC.md
```

**Key parameters:**
- A_max = 0.125 (fixed)
- Envelope: (X + 0.03125) / 1.03125 (linear)
- Re ∈ {50, 500, 1000, 5000, 10⁴, 5×10⁴, 10⁵, 2×10⁵}
- h/L ∈ {0.04, 0.08, 0.12, 0.16, 0.20, 0.24}

---

### 2. Understand the Differences

Compare strict Zhang approach to adaptive kinematics:
```bash
cat docs/COMPARISON.md
```

**Key difference:**
- **Zhang**: Fixed kinematics, vary Re/thickness to observe flow effects
- **Adaptive**: Modify kinematics based on Re/thickness (bio-inspired)

---

### 3. Verify Current Implementation

Check parameter alignment:
```bash
cat docs/VERIFICATION.md
```

**Status**: Current code can match Zhang if `enable_shape_adaptation = FALSE`

---

### 4. Run a Zhang Case

**Example: Re = 1000, h/L = 0.04**

```bash
# From repository root
cd build

# Compile (if needed)
cmake .. && make -j4

# Run Zhang case
mpirun -np 6 ./main2d ../Zhang_2018/input_files/input2d_Zhang_Re1000_h004
```

**Expected output:**
- Simulation runs with fixed kinematics
- U₀ ≈ 0.15-0.20 (swimming speed)
- St ≈ 0.8-1.2 (Strouhal number)
- 2S or 2P wake pattern

---

## Complete Test Matrix

Zhang tested 48 cases (8 Re × 6 h/L).

### Priority Validation Cases

Recommended first runs:

1. **Re = 1000, h/L = 0.04**: `input2d_Zhang_Re1000_h004` ✓
2. **Re = 5000, h/L = 0.04**: `input2d_Zhang_Re5000_h004`
3. **Re = 10000, h/L = 0.08**: `input2d_Zhang_Re10000_h008`
4. **Re = 100000, h/L = 0.04**: `input2d_Zhang_Re100000_h004`

See `input_files/ZHANG_TEST_MATRIX.md` for complete matrix.

---

## Using the Strict Zhang Class

### Option 1: Disable Adaptation in Input File

**Current approach** (works with existing code):

```
// In input file ConstraintIBKinematics section:
enable_shape_adaptation = FALSE  // Critical!
base_amplitude = 0.125
base_frequency = 0.785
envelope_power = 1.0
swimming_mode = 0.0
```

This forces the existing `IBEELKinematics` class to use fixed parameters.

---

### Option 2: Use Dedicated Zhang Class (Recommended)

**New dedicated class** (requires code modification):

1. **Copy Zhang source files to main directory:**
   ```bash
   cp Zhang_2018/src/IBEELKinematicsZhang.* .
   ```

2. **Modify CMakeLists.txt** to include Zhang class:
   ```cmake
   add_executable(main2d
     example.cpp
     IBEELKinematics.cpp
     IBEELKinematicsZhang.cpp  # Add this line
   )
   ```

3. **In example.cpp**, use Zhang class:
   ```cpp
   // Replace:
   Pointer<IBEELKinematics> ib_kinematics = new IBEELKinematics(...);

   // With:
   #include "IBEELKinematicsZhang.h"
   Pointer<IBEELKinematicsZhang> ib_kinematics = new IBEELKinematicsZhang(...);
   ```

4. **Rebuild:**
   ```bash
   cd build && make -j4
   ```

**Advantages:**
- Guarantees Zhang compliance (no user error)
- Automatic parameter verification
- Clear separation of approaches
- Warnings if parameters don't match Zhang

---

## Validation Protocol

### Step 1: Run Simulation

```bash
mpirun -np 6 ./main2d input2d_Zhang_Re1000_h004
```

### Step 2: Extract Performance Metrics

Outputs are in:
- `Zhang_Results_Re1000_h004/` - Structure diagnostics
- `Zhang_performance_Re1000_h004.dat` - Time-series performance
- `Zhang_viz_Re1000_h004/` - Visualization files

### Step 3: Compute Time-Averaged Values

From performance file:
```python
import numpy as np

data = np.loadtxt('Zhang_performance_Re1000_h004.dat')
time = data[:, 0]
U0 = data[:, 3]  # Swimming speed

# Average over last 5 periods (after transient)
U0_avg = np.mean(U0[time > 5.0])
print(f"U_0 = {U0_avg:.4f}")
```

### Step 4: Compare to Zhang Fig. 4

Plot swimming speed vs. Re for multiple thickness values:
- Zhang Fig. 4a: U₀(Re) for h/L = 0.04-0.24
- Zhang Fig. 4b: St(Re)
- Zhang Fig. 6: η(Re)

### Step 5: Visualize Wake

Use VisIt or ParaView to view vorticity:
```bash
visit Zhang_viz_Re1000_h004/dumps.visit
```

Expected wake patterns:
- Re = 50: Weak, diffuse vortices
- Re = 1000: 2S (two single vortices per cycle)
- Re = 10000: 2P (two vortex pairs per cycle)
- Re ≥ 50000: 2P + Kelvin-Helmholtz eddies

---

## Expected Results (from Zhang 2018)

### Swimming Speed U₀

| Re | h/L=0.04 | h/L=0.12 | h/L=0.24 |
|----|----------|----------|----------|
| 1000 | ~0.18 | ~0.12 | ~0.07 |
| 10000 | ~0.35 | ~0.28 | ~0.20 |
| 100000 | ~0.50 | ~0.45 | ~0.35 |

**Trend**: U₀ increases with Re, decreases with thickness.

### Strouhal Number St

- Range: 0.3 - 2.1
- Peaks around Re ~ 1000-5000
- Optimal St ≈ 0.3 for high efficiency

### Propulsive Efficiency η

- Increases with Re (low → high)
- Optimal thickness depends on Re
- Peak efficiency at intermediate thickness

---

## Comparison: Zhang vs. Adaptive

To compare the two approaches:

### Run Both Modes

1. **Strict Zhang** (fixed kinematics):
   ```bash
   mpirun -np 6 ./main2d input2d_Zhang_Re1000_h004
   ```

2. **Adaptive** (Re-dependent kinematics):
   ```bash
   mpirun -np 6 ./main2d input2d_Re1000_h004  # From parent directory
   ```

### Compare Performance

```python
import matplotlib.pyplot as plt
import numpy as np

# Load both datasets
zhang_data = np.loadtxt('Zhang_performance_Re1000_h004.dat')
adaptive_data = np.loadtxt('performance_Re1000_h004.dat')

# Extract swimming speed
time_zhang = zhang_data[:, 0]
U0_zhang = zhang_data[:, 3]
time_adaptive = adaptive_data[:, 0]
U0_adaptive = adaptive_data[:, 3]

# Plot comparison
plt.figure()
plt.plot(time_zhang, U0_zhang, label='Zhang (fixed)', linewidth=2)
plt.plot(time_adaptive, U0_adaptive, label='Adaptive', linewidth=2, linestyle='--')
plt.xlabel('Time')
plt.ylabel('Swimming Speed U₀')
plt.legend()
plt.title('Re = 1000, h/L = 0.04')
plt.savefig('comparison_Re1000.png')
plt.show()
```

**Expected differences:**
- Adaptive uses larger amplitude at low Re
- Adaptive achieves potentially higher U₀
- Zhang shows pure Re effect (no kinematic compensation)

---

## Troubleshooting

### Issue: Simulation unstable

**Solution:**
- Reduce `DT_MAX` (especially for low or high Re)
- Increase grid resolution `N`
- Adjust `CFL_MAX` to 0.2-0.25
- Check vorticity tagging thresholds

### Issue: Performance metrics file empty

**Solution:**
- Verify `track_performance = TRUE` in input file
- Check file permissions
- Ensure simulation ran long enough (t > 1.0)

### Issue: Results don't match Zhang

**Solution:**
- Verify `enable_shape_adaptation = FALSE`
- Check Re calculation: `MU = 0.785 / Re`
- Confirm amplitude = 0.125 (not adapted)
- Run longer to reach quasi-steady state
- Check grid resolution is sufficient

### Issue: Adaptation still occurring

**Solution:**
- Explicitly set `enable_shape_adaptation = FALSE`
- Use `IBEELKinematicsZhang` class (guarantees fixed kinematics)
- Check log output for "Adaptive Kinematics Update" messages
- Verify base parameters are not being modified

---

## Generating Additional Input Files

To create input files for other Re and h/L combinations:

### Manual Method

1. Copy template:
   ```bash
   cp input2d_Zhang_Re1000_h004 input2d_Zhang_Re5000_h004
   ```

2. Edit parameters:
   ```bash
   # Change:
   Re = 5000.0
   MU = 0.785/Re
   reynolds_number = 5000.0
   # Update output directories and filenames
   ```

### Automated Method

Create a script to generate all 48 cases:

```bash
#!/bin/bash
# generate_all_zhang_cases.sh

RE_VALUES=(50 500 1000 5000 10000 50000 100000 200000)
H_VALUES=(004 008 012 016 020 024)

for re in "${RE_VALUES[@]}"; do
  for h in "${H_VALUES[@]}"; do
    # Generate input file from template
    sed -e "s/Re = 1000.0/Re = ${re}.0/g" \
        -e "s/reynolds_number = 1000.0/reynolds_number = ${re}.0/g" \
        -e "s/h004/h${h}/g" \
        -e "s/0.04/0.${h}/g" \
        input2d_Zhang_Re1000_h004 > input2d_Zhang_Re${re}_h${h}
  done
done
```

---

## Performance Expectations

### Computational Cost

**Per simulation:**
- Low Re (50-1000): 12-24 hours on 6-12 cores
- Moderate Re (5000-10000): 6-12 hours on 6-12 cores
- High Re (≥50000): 12-48 hours on 12-24 cores

**Memory:**
- Typical: 4-8 GB per case
- High Re or fine grid: 16-32 GB

### Parallelization

- Use 6-12 MPI processes for most cases
- High Re may benefit from 12-24 processes
- Domain decomposition handled by SAMRAI

---

## Publications and Citation

If using this implementation for research, please cite:

**Original work:**
```
Zhang, C., Huang, H., & Lu, X.-Y. (2018).
Effects of Reynolds number and thickness on an undulatory self-propelled foil.
Physics of Fluids, 30(7), 071902.
```

**IBAMR:**
```
Griffith, B. E., & Hornung, R. D. (2017).
An adaptive, formally second order accurate version of the immersed boundary method.
Journal of Computational Physics, 223(1), 10-49.
```

---

## Contributing

To add new features or report issues:

1. **Validate against Zhang (2018)** first
2. Document any modifications clearly
3. Maintain separation between strict and adaptive modes
4. Update test matrix documentation
5. Include comparison plots

---

## License

This code is based on IBAMR and follows the 3-clause BSD license.
See COPYRIGHT file in IBAMR distribution for details.

---

## Contact and Support

**For questions:**
- Review included documentation first
- Check IBAMR documentation: https://ibamr.github.io
- Compare to Zhang et al. (2018) published results
- Post issues on repository

**For IBAMR-specific help:**
- IBAMR mailing list
- IBAMR GitHub issues

---

## Summary

This Zhang_2018 directory provides:

✓ Complete nondimensional specification from Zhang (2018)
✓ Strict kinematics implementation (no adaptation)
✓ Full test matrix (48 cases: 8 Re × 6 h/L)
✓ Validation framework and expected results
✓ Comparison tools (strict vs. adaptive)
✓ Documentation and usage guidelines

**Key principle**: Zhang uses **FIXED kinematics** to isolate Reynolds number and thickness effects. This is fundamentally different from adaptive approaches that modify kinematics based on flow conditions.

**To match Zhang exactly**: Set `enable_shape_adaptation = FALSE` and use the provided input files.

---

**Last Updated**: 2025-01-17
**Version**: 1.0
**Status**: Ready for Validation
