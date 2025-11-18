# Test 01: Smoke Test

## Purpose

This test validates the basic scalar transport infrastructure in IBAMR. It is the simplest test and should be run first to ensure that the fundamental components are working correctly.

## Test Description

- **Test Type**: Infrastructure validation
- **Physics**: Pure diffusion (no advection, no immersed boundaries)
- **Initial Condition**: Gaussian blob centered at origin
- **Boundary Conditions**: Homogeneous Dirichlet (C = 0 at all boundaries)
- **Domain**: 2D square [-1, 1] × [-1, 1]
- **Grid**: 64 × 64 uniform grid
- **Time**: t = 0 to 1.0 with dt = 0.01

## What This Test Checks

1. **Variable Registration**: Scalar field variable is properly registered
2. **Initialization**: Initial conditions are applied correctly
3. **Boundary Conditions**: BCs are enforced without errors
4. **Time Integration**: Solver advances through all time steps
5. **Stability**: No NaN/Inf values during simulation
6. **Negativity**: Minimal or no negative concentration values
7. **I/O**: Visualization data is written correctly

## Expected Behavior

The Gaussian blob should:
- Diffuse radially outward
- Decrease in amplitude
- Spread out over time
- Not produce negative values
- Not crash or produce NaN

## Pass/Fail Criteria

The test PASSES if:
- [x] Simulation completes all time steps without crashing
- [x] No NaN or Inf values detected
- [x] Negative values < 1% of total cells
- [x] Visualization files are created
- [x] Mass drift < 1e-6 (for no-flux BCs, mass should be conserved)

The test FAILS if:
- [ ] Simulation crashes
- [ ] NaN or Inf values are detected
- [ ] Excessive negative values
- [ ] Cannot read/write data

## Building and Running

### Build

```bash
cd Test01_SmokeTest
mkdir build && cd build
cmake ..
make
```

### Run

```bash
./test01_smoke ../input2d
```

Or with MPI:
```bash
mpirun -np 4 ./test01_smoke ../input2d
```

## Expected Output

The test will produce:

1. **Console Output**:
   ```
   ================================================================================
     TEST 1: Smoke Test - Scalar Transport Infrastructure
   ================================================================================
     Start time: 2025-11-17 XX:XX:XX
   ================================================================================

   [INFO] Initializing application...
   [INFO] Creating integrators...
   [INFO] Creating grid geometry...
   ...
   [INFO] Initial statistics:
     C_min = 0.000000
     C_max = 1.000000
     Total mass = 0.314159

   ----------------------------------------
     Time Integration
   ----------------------------------------
   Iteration 0 / 100, t = 0.000000
   Iteration 10 / 100, t = 0.100000
   ...

   ----------------------------------------
     Final Statistics
   ----------------------------------------
     C_min = 0.000000
     C_max = 0.876543
     Total mass = 0.314157
     Mass drift = 6.37e-06

   ----------------------------------------
     Test Verdict
   ----------------------------------------
     [PASS] No crashes
     [PASS] No NaN/Inf
     [PASS] Negatives < 1% of cells
     [PASS] Completed all time steps

   Elapsed time: 12.34s

   ================================================================================
     TEST RESULT: PASSED
     All checks passed
     End time: 2025-11-17 XX:XX:XX
   ================================================================================
   ```

2. **Files Created**:
   - `test01.log` - IBAMR log file
   - `test01_results.txt` - Test results summary
   - `viz_test01/` - Directory with VisIt visualization files
     - `dumps.visit` - VisIt session file
     - `visit*.vtk` - VTK files for each dump

3. **Result File** (`test01_results.txt`):
   ```
   Test Results Log
   Generated: 2025-11-17 XX:XX:XX
   ============================================================

   End time: 1.000000
   Time step: 0.010000
   Number of steps: 100
   Diffusion coefficient: 0.001000
   Domain X: -1.000000 to 1.000000
   Domain Y: -1.000000 to 1.000000
   Initial C_min: 0.000000
   Initial C_max: 1.000000
   Initial mass: 0.314159
   Final C_min: 0.000000
   Final C_max: 0.876543
   Final mass: 0.314157
   Mass drift: 6.37e-06

   Test: Test 01: Smoke Test
   Result: PASSED
   Time: 2025-11-17 XX:XX:XX
   ------------------------------------------------------------

   Total elapsed time: 12.34s
   ```

## Visualization

To visualize results with VisIt:

```bash
visit -o viz_test01/dumps.visit
```

Or with ParaView:
```bash
paraview viz_test01/visit*.vtk
```

You should see the Gaussian blob diffusing outward over time.

## Common Issues

### Issue: "IBAMR not found"
**Solution**: Set `IBAMR_ROOT` environment variable:
```bash
export IBAMR_ROOT=/path/to/ibamr
cmake ..
```

### Issue: Compilation errors about missing headers
**Solution**: Ensure IBAMR, SAMRAI, and PETSc are properly installed and paths are set.

### Issue: Negative concentrations
**Diagnosis**: This can happen if:
- Time step is too large (reduce DT)
- Diffusion coefficient is too large
- Grid is too coarse
- Initial conditions are too sharp

**Solution**: Reduce DT or increase grid resolution.

### Issue: Mass not conserved
**Diagnosis**: Check boundary conditions. Dirichlet BCs allow mass to leave domain.

**Solution**: This is expected with Dirichlet BCs. For true conservation, use Neumann BCs or periodic BCs.

## Modifications

To test different scenarios, modify `input2d`:

1. **Periodic BCs**: Change to `periodic_dimension = 1, 1`
2. **Finer grid**: Change to `domain_boxes = [ (0,0) , (127,127) ]`
3. **Longer time**: Change `END_TIME = 2.0`
4. **Different IC**: Modify `OdorInitialConditions` function

## Next Steps

After Test 01 passes:
1. Run Test 02 (Pure Diffusion) for quantitative validation
2. Run Test 03 (Pure Advection) to test advection operator
3. Run Test 04 (MMS) for combined advection-diffusion validation

## References

- IBAMR Documentation: https://ibamr.github.io/docs/
- IBAMR Examples: https://github.com/IBAMR/IBAMR/tree/master/examples

---

**Test Status**: Implemented
**Last Updated**: 2025-11-17
