# IBAMR C++ Test Suite - Quick Start Guide

## What You Have

A complete **C++ implementation framework** for 14 Verification & Validation (V&V) tests for scalar transport using the IBAMR framework.

## Status

- **Common Utilities**: ✅ Fully implemented
- **Test01 (Smoke Test)**: ✅ Fully implemented
- **Tests 02-14**: ✅ Template structures created (ready for implementation)
- **Build System**: ✅ Complete CMake setup
- **Scripts**: ✅ Build and run automation

## Quick Start in 3 Steps

### 1. Set Environment

```bash
export IBAMR_ROOT=/path/to/your/ibamr/installation
```

### 2. Build Tests

```bash
cd IBAMR_CPP_Tests
./build_all_tests.sh
```

### 3. Run Test 1

```bash
cd Test01_SmokeTest
../build/test01_smoke input2d
```

## What's Implemented

### ✅ Fully Implemented

1. **Common Utilities** (`common/`)
   - `AnalyticalSolutions.h/cpp` - Analytical solutions for validation
   - `ErrorCalculator.h/cpp` - L2/Linf error computation, convergence rates
   - `TestUtilities.h/cpp` - Logging, timing, formatting utilities

2. **Test 01: Smoke Test** (`Test01_SmokeTest/`)
   - Complete C++ implementation
   - Input file configured
   - Comprehensive README
   - Validates: Variable registration, BCs, I/O, stability

### ✅ Template Structures (Ready for Implementation)

Tests 02-14 have:
- Template `main.cpp` files with proper includes
- Basic `input2d` configuration files
- README templates with descriptions
- Ready to be filled in following Test01 pattern

## Directory Structure

```
IBAMR_CPP_Tests/
├── README.md              ← Comprehensive documentation
├── QUICK_START.md         ← This file
├── CMakeLists.txt         ← Master build file
├── build_all_tests.sh     ← Build automation
├── run_all_tests.sh       ← Run automation
│
├── common/                ← Shared utilities (COMPLETE)
│   ├── include/
│   │   ├── AnalyticalSolutions.h
│   │   ├── ErrorCalculator.h
│   │   └── TestUtilities.h
│   └── src/
│       ├── AnalyticalSolutions.cpp
│       ├── ErrorCalculator.cpp
│       └── TestUtilities.cpp
│
├── Test01_SmokeTest/      ← FULLY IMPLEMENTED
│   ├── main.cpp           ← Complete implementation
│   ├── input2d            ← Configured
│   └── README.md          ← Detailed guide
│
├── Test02_Diffusion_Analytic/  ← TEMPLATE
├── Test03_Advection_Analytic/  ← TEMPLATE
├── Test04_MMS/                 ← TEMPLATE
├── Test05_Discontinuous/       ← TEMPLATE
├── Test06_MassConservation/    ← TEMPLATE
├── Test07_BCs/                 ← TEMPLATE
├── Test08_SphereSource/        ← TEMPLATE
├── Test09_HighSc/              ← TEMPLATE
├── Test10_MovingIB/            ← TEMPLATE
├── Test11_AMR/                 ← TEMPLATE
├── Test12_TimeStep/            ← TEMPLATE
├── Test13_LongRun/             ← TEMPLATE
└── Test14_Benchmarks/          ← TEMPLATE
```

## Next Steps

### Immediate (Today)

1. **Install IBAMR** (if not already installed)
   - Follow: https://ibamr.github.io/installing/

2. **Set environment variable**:
   ```bash
   export IBAMR_ROOT=/path/to/ibamr
   ```

3. **Build and test**:
   ```bash
   cd IBAMR_CPP_Tests
   ./build_all_tests.sh
   cd Test01_SmokeTest
   ../build/test01_smoke input2d
   ```

4. **Check results**:
   - Look for "TEST PASSED" in console output
   - Check `test01_results.txt` for detailed results
   - Visualize with: `visit viz_test01/dumps.visit`

### This Week

1. **Implement Test02**: Pure Diffusion
   - Use Test01 as template
   - Add analytical solution comparison
   - Implement convergence study

2. **Implement Test03**: Pure Advection
   - Add uniform velocity field
   - Test advection operator

3. **Implement Test04**: MMS
   - Add manufactured source term
   - Verify combined advection-diffusion

### This Month

Complete Tier 1 tests (Tests 1-6):
- ✅ Test01: Smoke Test
- ⏩ Test02: Pure Diffusion
- ⏩ Test03: Pure Advection
- ⏩ Test04: MMS
- ⏩ Test05: Discontinuous
- ⏩ Test06: Mass Conservation

## Implementation Pattern

Each test follows this structure (see `Test01_SmokeTest/main.cpp`):

```cpp
int main(int argc, char* argv[])
{
    // 1. Initialize IBAMR/SAMRAI
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // 2. Create integrators
    // - INSStaggeredHierarchyIntegrator (for velocity)
    // - AdvDiffHierarchyIntegrator (for scalar)

    // 3. Set up grid and hierarchy

    // 4. Set initial/boundary conditions

    // 5. Time integration loop

    // 6. Compute errors / validate results

    // 7. Report pass/fail

    return test_passed ? 0 : 1;
}
```

## Common Utilities Usage

### Analytical Solutions
```cpp
#include "AnalyticalSolutions.h"

// Gaussian diffusion
double C_exact = AnalyticalSolutions::gaussianDiffusion2D(
    x, y, t, kappa, x0, y0, C0);

// Manufactured solution
double C_mms = AnalyticalSolutions::manufacturedSolution2D(x, y, t);
```

### Error Calculation
```cpp
#include "ErrorCalculator.h"

// Compute L2 error
double L2_error = ErrorCalculator::computeL2Error(
    patch_hierarchy, computed_idx, exact_idx);

// Compute convergence rate
double rate = ErrorCalculator::computeConvergenceRate(
    grid_spacings, errors);

// Check for problems
bool has_nan = ErrorCalculator::checkForNaNInf(patch_hierarchy, var_idx);
int neg_count = ErrorCalculator::checkForNegatives(patch_hierarchy, var_idx);
```

### Test Utilities
```cpp
#include "TestUtilities.h"

// Logging
TestUtils::printTestHeader("My Test", 1);
TestUtils::printProgress("Computing errors...");
ResultLogger logger("results.txt");
logger.logError("L2 error", 1.23e-6);

// Timing
Timer timer;
timer.start();
// ... do work ...
timer.stop();
std::cout << "Elapsed: " << timer.getFormattedTime() << std::endl;
```

## Building Individual Tests

```bash
# Option 1: Use master build
cd IBAMR_CPP_Tests
./build_all_tests.sh

# Option 2: Build specific test
cd Test01_SmokeTest
mkdir build && cd build
cmake ..
make
```

## Running Tests

```bash
# Single test
cd Test01_SmokeTest
../build/test01_smoke input2d

# With MPI
mpirun -np 4 ../build/test01_smoke input2d

# All tests
cd IBAMR_CPP_Tests
./run_all_tests.sh
```

## Troubleshooting

### "IBAMR not found"
```bash
export IBAMR_ROOT=/path/to/ibamr
export CMAKE_PREFIX_PATH=$IBAMR_ROOT:$CMAKE_PREFIX_PATH
```

### "undefined reference to..."
Ensure IBAMR was built with same compiler and dependencies.

### "No such file: input2d"
Run from test directory: `cd Test01_SmokeTest`

## Customizing Tests

### Modify Grid Resolution
In `input2d`:
```
domain_boxes = [ (0,0) , (127,127) ]  // Change from (63,63)
```

### Change Time Step
```
DT = 0.005  // Smaller time step
END_TIME = 2.0  // Longer simulation
```

### Modify Initial Conditions
```
OdorInitialConditions {
   function = "sin(pi*X_0)*sin(pi*X_1)"  // Different IC
}
```

## Expected Test Progression

### Phase 1: Infrastructure (Week 1)
- ✅ Test01: Smoke Test → Validates infrastructure works

### Phase 2: Operators (Week 2)
- Test02: Pure Diffusion → Validates diffusion operator
- Test03: Pure Advection → Validates advection operator

### Phase 3: Verification (Week 3)
- Test04: MMS → Verifies combined operators
- Test05: Discontinuous → Tests stability
- Test06: Mass Conservation → Validates conservation

### Phase 4: Physical Validation (Week 4)
- Test07: BCs → Boundary conditions
- Test08: Sphere Source → Literature comparison
- Test09: High Sc → Extreme parameters
- Test10: Moving IB → IB coupling

### Phase 5: Production (Week 5+)
- Test11: AMR → Adaptivity
- Test12: Time-step → CFL sensitivity
- Test13: Long Run → Long-term stability
- Test14: Benchmarks → Full validation

## Key Files

- **README.md** - Comprehensive documentation
- **CMakeLists.txt** - Build configuration
- **common/include/*.h** - Utility headers
- **Test01_SmokeTest/main.cpp** - Reference implementation
- **build_all_tests.sh** - Build automation
- **run_all_tests.sh** - Test automation

## Getting Help

1. **IBAMR Documentation**: https://ibamr.github.io/docs/
2. **IBAMR Examples**: https://github.com/IBAMR/IBAMR/tree/master/examples
3. **Test01 README**: `Test01_SmokeTest/README.md` for detailed guide
4. **Main README**: `README.md` for comprehensive documentation

## Summary

✅ **What's Done**:
- Complete test framework structure
- Common utilities fully implemented
- Test01 fully implemented
- All 14 tests have template structures
- Build system configured
- Automation scripts ready

⏩ **Next**:
1. Install/configure IBAMR
2. Build test suite
3. Run Test01 to validate setup
4. Implement remaining tests (02-14) following Test01 pattern

---

**Repository**: Four_fish_school/IBAMR_CPP_Tests
**Status**: Framework Complete, Ready for Implementation
**Last Updated**: 2025-11-17
