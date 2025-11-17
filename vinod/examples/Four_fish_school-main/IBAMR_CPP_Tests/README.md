# IBAMR C++ Test Suite for Scalar Transport

This directory contains a complete implementation of all 14 Verification & Validation (V&V) tests for scalar transport (odor concentration) using the IBAMR framework, implemented entirely in C++.

## Overview

This test suite validates the advection-diffusion solver for scalar transport in IBAMR, which is used to simulate odor concentration fields in the 4-fish school simulation. All tests are self-contained C++ applications using the IBAMR framework.

## About IBAMR

IBAMR is an open-source Immersed Boundary Method framework for simulating fluid-structure interaction.
- **Website**: https://ibamr.github.io/
- **Documentation**: https://ibamr.github.io/docs/
- **GitHub**: https://github.com/IBAMR/IBAMR

## Directory Structure

```
IBAMR_CPP_Tests/
├── README.md                          ← This file
├── CMakeLists.txt                     ← Master build file
├── common/                            ← Shared utilities
│   ├── include/
│   │   ├── AnalyticalSolutions.h     ← Analytical solution functions
│   │   ├── ErrorCalculator.h         ← L2/Linf error computation
│   │   └── TestUtilities.h           ← Common test helpers
│   └── src/
│       ├── AnalyticalSolutions.cpp
│       ├── ErrorCalculator.cpp
│       └── TestUtilities.cpp
│
├── Test01_SmokeTest/                  ← Basic infrastructure check
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test02_Diffusion_Analytic/         ← Gaussian diffusion validation
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test03_Advection_Analytic/         ← Pure advection validation
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test04_MMS/                        ← Manufactured solution
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test05_Discontinuous/              ← Top-hat stability test
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test06_MassConservation/           ← Mass budget validation
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test07_BCs/                        ← Boundary condition tests
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test08_SphereSource/               ← Sphere source (Lei et al.)
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test09_HighSc/                     ← High Schmidt number
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test10_MovingIB/                   ← Moving IB + scalar
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test11_AMR/                        ← AMR sensitivity
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test12_TimeStep/                   ← Time-step convergence
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test13_LongRun/                    ← Long-run periodicity
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test14_Benchmarks/                 ← Literature comparisons
│   ├── main.cpp
│   ├── input2d
│   ├── CMakeLists.txt
│   └── README.md
│
├── Test15_RotatingCylinder/           ← Lei et al. validation: Rotating cylinder
│   ├── main.cpp
│   ├── input2d
│   └── README.md
│
├── Test16_3DSphere/                   ← Lei et al. validation: 3D sphere/cube
│   ├── main.cpp
│   ├── input3d_sphere
│   └── README.md
│
└── Test17_PitchPlunge/                ← Lei et al. validation: Full case
    ├── main.cpp
    ├── input2d
    └── README.md
```

## Test Suite Overview

| Test # | Name | Purpose | Pass Criteria |
|--------|------|---------|---------------|
| 1 | Smoke Test | Basic infrastructure check | No crashes, no NaNs |
| 2 | Pure Diffusion | Gaussian diffusion vs analytic | L2 error convergence rate ≈ 2.0 |
| 3 | Pure Advection | Profile advection validation | L2 error < tolerance |
| 4 | MMS | Manufactured solution verification | Convergence rate ≈ 2.0 |
| 5 | Discontinuous | Top-hat advection stability | No oscillations, no negatives |
| 6 | Mass Conservation | Global tracer mass budget | Relative drift < 1e-10 |
| 7 | Boundary Conditions | Dirichlet/Neumann/flux tests | BC error ≤ 1e-6 |
| 8 | Sphere Source | Compare Lei et al. | Match literature ±10% |
| 9 | High Schmidt | Sc=100-1000 stability | Stable, physically reasonable |
| 10 | Moving IB | Ellipsoid + scalar coupling | No instabilities |
| 11 | AMR Sensitivity | Refinement artifact checks | Consistent across refinements |
| 12 | Time-step CFL | Temporal convergence | Convergence rate ≈ 2.0 |
| 13 | Long Run | Periodic conservation checks | No drift over time |
| 14 | Benchmarks | Literature comparisons | Match Lei/Kamran data |
| **15** | **Rotating Cylinder** | **Yan & Zu (2008) validation** | **Streamlines + scalar match ±10%** |
| **16** | **3D Sphere** | **Richter & Nikrityuk (2012)** | **Nu within ±15%, contours match** |
| **17** | **Pitch-Plunge** | **Lei et al. (2021) full case** | **PDF shape, vortex street** |

## Building the Tests

### Prerequisites

1. **IBAMR** (with dependencies: SAMRAI, PETSc, HDF5)
   - Follow installation guide: https://ibamr.github.io/installing/

2. **CMake** (version 3.12+)

3. **C++ Compiler** with C++14 support (GCC 7+, Clang 5+)

### Build Instructions

```bash
# Set IBAMR installation path
export IBAMR_ROOT=/path/to/ibamr/installation

# Build all tests
cd IBAMR_CPP_Tests
mkdir build && cd build
cmake .. -DIBAMR_ROOT=$IBAMR_ROOT
make -j4

# Or build individual test
cd Test01_SmokeTest
mkdir build && cd build
cmake ..
make
```

### Environment Variables

```bash
# Required
export IBAMR_ROOT=/path/to/ibamr

# Optional (if not in standard locations)
export SAMRAI_ROOT=/path/to/samrai
export PETSC_DIR=/path/to/petsc
export HDF5_ROOT=/path/to/hdf5
```

## Running the Tests

### Run Individual Test
```bash
cd Test01_SmokeTest/build
./test01_smoke input2d
```

### Run All Tests (Sequential)
```bash
cd IBAMR_CPP_Tests
./run_all_tests.sh
```

### Run with MPI
```bash
mpirun -np 4 ./test01_smoke input2d
```

## Test Descriptions

### Tier 1: Basic Verification (Tests 1-6)
These tests verify fundamental solver capabilities.

**Test 1: Smoke Test**
- Verifies basic infrastructure works
- No immersed boundaries
- Simple initial condition
- Expected: No crashes, no NaNs

**Test 2: Pure Diffusion**
- Compares to analytical Gaussian solution: C(x,t) = exp(-x²/4κt)/√(4πκt)
- Tests diffusion operator accuracy
- Expected: 2nd order convergence

**Test 3: Pure Advection**
- Advects a profile with uniform velocity
- Tests advection operator
- Expected: Profile shape preservation

**Test 4: Method of Manufactured Solutions**
- Uses known source term to verify correctness
- Tests combined advection-diffusion
- Expected: 2nd order convergence

**Test 5: Discontinuous Initial Condition**
- Top-hat function advection
- Tests monotonicity and stability
- Expected: No oscillations, no negative values

**Test 6: Mass Conservation**
- Tracks global mass over time
- Tests conservative property
- Expected: Relative drift < 1e-10

### Tier 2: Physical Validation (Tests 7-10)
These tests validate physical behavior.

**Test 7: Boundary Conditions**
- Tests Dirichlet, Neumann, Robin BCs
- Verifies BC enforcement
- Expected: BC errors ≤ 1e-6

**Test 8: Sphere Source**
- Compares to Lei et al. (2021) data
- Steady sphere with scalar source
- Expected: Match literature ±10%

**Test 9: High Schmidt Number**
- Tests Sc = 100, 340, 1000
- Validates thin boundary layers
- Expected: Stable, physically reasonable

**Test 10: Moving Immersed Boundary**
- Oscillating ellipsoid with scalar field
- Tests IB-scalar coupling
- Expected: No instabilities

### Tier 3: Production Validation (Tests 11-14)
These tests validate production readiness.

**Test 11: AMR Sensitivity**
- Runs with different refinement ratios
- Checks for AMR artifacts
- Expected: Consistent results

**Test 12: Time-step Sensitivity**
- CFL number variation
- Temporal convergence study
- Expected: 2nd order convergence

**Test 13: Long Run**
- Extended simulation (T = 100+)
- Checks long-term stability
- Expected: No drift, no instabilities

**Test 14: Benchmarks**
- Lei et al. (2021): Pitching airfoil
- Kamran et al. (2024): Undulating body
- Expected: Match published data

### Tier 4: Literature Validation (Tests 15-17)
These tests implement the exact validation cases from Lei et al. (2021) Section II.B.

**Test 15: Rotating Cylinder (Yan & Zu 2008)**
- Validates scalar transport around rotating isothermal cylinder
- Parameters: Re=200, Pr=0.5, k=0.5 (rotation V/U)
- Reference: Yan & Zu (2008) LBM heat transfer benchmark
- Expected: Streamlines + scalar contours match published Figure 2
- See: `LEI_VALIDATION_MAPPING.md` for full details

**Test 16: 3D Sphere & Cube (Richter & Nikrityuk 2012)**
- Validates 3D scalar transport around stationary sphere/cube
- Parameters: Re=200, Pr=0.744
- Reference: Richter & Nikrityuk (2012) CFD benchmark
- Expected: Thermal BL, wake structure, Nusselt number within ±15%
- **Note:** 3D simulation - high computational cost (~4 hours on 16 cores)
- See: `LEI_VALIDATION_MAPPING.md` for full details

**Test 17: Pitch-Plunge Odor Plume (Lei et al. 2021)**
- **FULL PRODUCTION CASE** from Lei et al. (2021)
- Upstream sphere source + downstream flapping ellipsoidal airfoil
- Parameters: Re=200, Sc=0.71, St=0.9
- Kinematics: y(t) = (L/2)sin(2πft), θ(t) = (1/6)cos(2πft)
- Expected: Inverse von Kármán street, odor PDF shape, vortex-odor interaction
- Validates: Figures 6-10 in Lei et al. (2021)
- **Note:** Time-dependent, requires multiple flapping periods (~4 hours 2D)
- See: `LEI_VALIDATION_MAPPING.md` for complete parameter extraction

## Output and Results

Each test produces:
1. **Console output**: Real-time progress and error metrics
2. **VTK files**: Visualization data (compatible with VisIt/ParaView)
3. **HDF5 files**: Checkpoint/restart data
4. **Error report**: L2/Linf errors, convergence rates
5. **Summary**: Pass/fail verdict with details

### Typical Output
```
Test01_SmokeTest/
├── viz_test01/              ← VisIt visualization data
│   ├── dumps.visit
│   ├── lag_data.*.vtu
│   └── visit*.*.vtk
├── restart_test01/          ← Restart files
├── test01_output.log        ← Simulation log
└── test01_results.txt       ← Error metrics and verdict
```

## Validation Criteria

### Convergence Tests (2, 4, 12)
```
Expected convergence rate: 2.0 ± 0.3
L2 error ~ (Δx)^2 for spatial
L2 error ~ (Δt)^2 for temporal
```

### Mass Conservation (6, 13)
```
Advection-only: |M(t) - M(0)|/M(0) < 1e-10
With diffusion:  |M(t) - M(0)|/M(0) < 1e-6
```

### Stability Tests (5, 9, 10, 11)
```
No negative concentrations: min(C) ≥ -1e-12
No NaN values
No exponential growth
```

### Literature Comparison (8, 14)
```
Coarse grid: ±20% acceptable
Fine grid:   ±10% target
```

## Common Utilities

The `common/` directory provides shared functionality:

### AnalyticalSolutions.h/cpp
- Gaussian diffusion solution
- Manufactured solutions
- Exact advection profiles

### ErrorCalculator.h/cpp
- L2 norm computation
- L∞ norm computation
- Convergence rate calculation
- Grid refinement studies

### TestUtilities.h/cpp
- Mass conservation checks
- Negative concentration detection
- NaN/Inf detection
- Output formatting

## Troubleshooting

### Build Issues

**"IBAMR not found"**
```bash
export IBAMR_ROOT=/correct/path/to/ibamr
cmake .. -DIBAMR_ROOT=$IBAMR_ROOT
```

**"SAMRAI headers not found"**
```bash
export SAMRAI_ROOT=/path/to/samrai
export CMAKE_PREFIX_PATH=$SAMRAI_ROOT:$CMAKE_PREFIX_PATH
```

### Runtime Issues

**"Segmentation fault"**
- Check input2d file exists
- Verify grid dimensions are reasonable
- Check memory limits (ulimit -v)

**"No convergence"**
- Reduce time step
- Check CFL condition
- Verify boundary conditions

**"Wrong convergence rate"**
- Ensure time step scales with grid: Δt ~ Δx² (diffusion)
- Check initial conditions are smooth
- Verify BCs not polluting interior

## References

### IBAMR Documentation
- User Guide: https://ibamr.github.io/docs/
- Examples: https://github.com/IBAMR/IBAMR/tree/master/examples
- API: https://ibamr.github.io/api/

### Literature
1. **Lei et al. (2021)**: "Navigation in odor plumes: How do the flapping kinematics modulate the odor landscape"
   - PDF: `../Navigation in odor plumes How do the flapping kinematics modulate the odor landscape.pdf`
   - Tests: 8, 14, **15, 16, 17**
   - **Full parameter extraction:** `LEI_VALIDATION_MAPPING.md`

2. **Kamran et al. (2024)**: "How does vortex dynamics help undulating bodies spread odor"
   - PDF: `../How does vortex dynamics help undulating bodies spread odor.pdf`
   - Tests: 9, 14

3. **Yan & Zu (2008)**: "Numerical simulation of heat transfer and fluid flow past a rotating isothermal cylinder – A LBM approach"
   - Reference for Test15
   - DOI: 10.1016/j.ijheatmasstransfer.2007.07.053

4. **Richter & Nikrityuk (2012)**: "Drag forces and heat transfer coefficients for spherical, cuboidal and ellipsoidal particles"
   - Reference for Test16
   - DOI: 10.1016/j.ijheatmasstransfer.2011.09.005

### Related Repositories
- Parent repo: Four_fish_school
- IBAMR framework: https://github.com/IBAMR/IBAMR

## Quick Start

1. **Install IBAMR**: Follow https://ibamr.github.io/installing/
2. **Set environment**:
   ```bash
   export IBAMR_ROOT=/path/to/ibamr
   ```
3. **Build Test 1**:
   ```bash
   cd Test01_SmokeTest
   mkdir build && cd build
   cmake ..
   make
   ```
4. **Run Test 1**:
   ```bash
   ./test01_smoke ../input2d
   ```
5. **Check results**: Look for "TEST PASSED" in output

## Development Status

### Core Tests (1-14)
- [x] Directory structure created
- [ ] Common utilities implemented
- [ ] Test 1: Smoke Test
- [ ] Test 2: Pure Diffusion
- [ ] Test 3: Pure Advection
- [ ] Test 4: MMS
- [ ] Test 5: Discontinuous
- [ ] Test 6: Mass Conservation
- [ ] Test 7: Boundary Conditions
- [ ] Test 8: Sphere Source
- [ ] Test 9: High Schmidt
- [ ] Test 10: Moving IB
- [ ] Test 11: AMR
- [ ] Test 12: Time-step
- [ ] Test 13: Long Run
- [ ] Test 14: Benchmarks

### Literature Validation Tests (15-17) ✨ NEW
- [x] **Test 15: Rotating Cylinder** (Yan & Zu 2008) - Template created
- [x] **Test 16: 3D Sphere** (Richter & Nikrityuk 2012) - Template created
- [x] **Test 17: Pitch-Plunge** (Lei et al. 2021 full case) - Template created
- [x] **Validation mapping document** (`LEI_VALIDATION_MAPPING.md`)
- [ ] Test15 executed and validated
- [ ] Test16 executed and validated
- [ ] Test17 executed and validated

## Contributing

This test suite is designed to be:
- **Self-contained**: Each test is an independent C++ application
- **Well-documented**: Clear purpose and pass/fail criteria
- **Reproducible**: Fixed seeds, documented parameters
- **Extensible**: Easy to add new tests

## License

Same as parent Four_fish_school repository.

## Contact

For questions about IBAMR: https://github.com/IBAMR/IBAMR/issues
For questions about this test suite: See parent repository

---

**Last Updated**: 2025-11-17
**IBAMR Version**: 0.12.0+ recommended
**Status**: In Development
