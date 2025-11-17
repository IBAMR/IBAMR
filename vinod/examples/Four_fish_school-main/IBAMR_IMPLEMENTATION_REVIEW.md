# IBAMR Implementation Review for Odor Dynamics
## Comprehensive Analysis Against Test Plan Requirements

**Date:** 2025-11-17
**Branch:** `claude/review-ibamr-implementation-01FQ34vzZ7obX15xZ3zNcEjj`
**IBAMR Documentation:** https://ibamr.github.io/about

---

## Executive Summary

### Overall Status: 🟡 Framework Complete, Implementation Pending

The IBAMR-based scalar transport implementation has:

✅ **COMPLETE:**
- Comprehensive test framework architecture (14 tests)
- Fully implemented common utilities (896 lines)
- Reference implementation (Test01 - 363 lines)
- Proper IBAMR syntax and API usage
- Analytical solutions for validation

⏳ **IN PROGRESS:**
- Tests 02-14 are template skeletons (67 lines each)
- Need full implementation following Test01 pattern

---

## 1. IBAMR Syntax & API Review

### ✅ CORRECT Implementation in Test01

#### Initialization (Lines 38-45, Test01_SmokeTest/main.cpp:38)
```cpp
IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
```
**Status:** ✅ Correct
**Analysis:** Proper IBAMR/SAMRAI/PETSc initialization with MPI communicator.

#### Integrator Setup (Lines 84-95, Test01_SmokeTest/main.cpp:84)
```cpp
Pointer<INSHierarchyIntegrator> navier_stokes_integrator =
    new INSStaggeredHierarchyIntegrator(...);

Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator =
    new AdvDiffHierarchyIntegrator(...);

adv_diff_integrator->setAdvectionVelocity(
    navier_stokes_integrator->getAdvectionVelocityVariable());
navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
```
**Status:** ✅ Correct
**Analysis:**
- Proper two-way coupling between velocity and scalar fields
- Follows IBAMR best practices for coupled solvers
- Registration ensures synchronized time stepping

#### Time Integration Scheme (input2d:97)
```
diffusion_time_stepping_type = "CRANK_NICOLSON"
convective_time_stepping_type = "ADAMS_BASHFORTH"
```
**Status:** ✅ EXCELLENT
**Analysis:**
- Crank-Nicolson is 2nd-order implicit (perfect for high Schmidt numbers)
- Matches test plan requirement for implicit diffusion (Section 5.2)
- Will handle Sc = 340 (water) and Sc = 0.7 (air) per arXiv 2408.16136v1

#### Boundary Conditions (Lines 150-157, Test01_SmokeTest/main.cpp:150)
```cpp
RobinBcCoefStrategy<NDIM>* C_bc_coef = nullptr;
if (!(periodic_shift.min() > 0) && input_db->keyExists("OdorBcCoefs"))
{
    C_bc_coef = new muParserRobinBcCoefs(...);
    adv_diff_integrator->setPhysicalBcCoef(C_bc_coef);
}
```
**Status:** ✅ Correct
**Analysis:**
- Robin BC formulation: `a*C + b*dC/dn = g`
- Supports Dirichlet (a=1, b=0), Neumann (a=0, b=1), Robin (mixed)
- Matches test plan requirements (Section 2)

#### Grid Hierarchy (Lines 102-131, Test01_SmokeTest/main.cpp:102)
```cpp
Pointer<CartesianGridGeometry<NDIM>> grid_geometry = ...;
Pointer<PatchHierarchy<NDIM>> patch_hierarchy = ...;
Pointer<BergerRigoutsos<NDIM>> box_generator = ...;
```
**Status:** ✅ Correct
**Analysis:**
- Proper SAMRAI grid hierarchy construction
- BergerRigoutsos clustering algorithm for AMR
- Supports test plan Section 7 (AMR tests)

---

## 2. Common Utilities Analysis

### Analytical Solutions (AnalyticalSolutions.cpp - 176 lines)

#### ✅ Gaussian Diffusion 1D/2D/3D (Lines 7-65)
```cpp
double gaussianDiffusion2D(double x, double y, double t, double kappa,
                          double x0, double y0, double C0)
{
    double r_squared = (x-x0)*(x-x0) + (y-y0)*(y-y0);
    double denom = 4.0 * PI * kappa * t;
    double exponent = -r_squared / (4.0 * kappa * t);
    return (C0 / denom) * std::exp(exponent);
}
```
**Validation:** ✅ Matches analytical solution
**Formula:** `C(x,y,t) = (C₀/4πκt) * exp(-r²/4κt)`
**Test Plan:** Section 1.1 (Pure Diffusion Test)

#### ✅ Method of Manufactured Solutions (Lines 93-121)
```cpp
double manufacturedSolution2D(double x, double y, double t) {
    return std::exp(-t) * std::sin(PI * x) * std::sin(PI * y);
}

double manufacturedSource2D(double x, double y, double t,
                           double kappa, double u, double v) {
    // S = dC/dt - kappa*Laplacian(C) + u*dC/dx + v*dC/dy
    double laplacian_C = -2.0 * PI * PI * C;
    ...
}
```
**Validation:** ✅ Correct partial derivatives
**Test Plan:** Section 1.4 & Section 13 (MMS test)
**Example Match:** User's suggested MMS: `C(x,y,t) = sin(πx)sin(πy)e^(-2π²αt)`

#### ✅ Sphere Source 3D (Lines 139-155)
```cpp
double sphereSource3D(...) {
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    if (r < R + EPSILON) {
        return Q / (4.0 * PI * kappa * R);
    }
    return Q / (4.0 * PI * kappa * r);  // 1/r decay
}
```
**Validation:** ✅ Matches steady-state Laplace solution
**Formula:** `C(r) = Q/(4πκr)` for r > R
**Test Plan:** Section 3.1 (Steady point/spherical source)

### Error Calculator (ErrorCalculator.cpp - 319 lines)

#### ✅ L2/Linf Error Computation
**Functions:** `computeL2Error()`, `computeLinfError()`, `computeL1Error()`
**Implementation:** Uses SAMRAI `HierarchyCellDataOpsReal`
**Test Plan:** All tests require error norms

#### ✅ Convergence Rate Calculation
```cpp
double computeConvergenceRate(const std::vector<double>& grid_spacings,
                             const std::vector<double>& errors)
{
    // Log-log linear regression: log(error) = p*log(h) + C
    linearRegression(log_h, log_err, slope, intercept);
    return slope;  // Returns p (convergence order)
}
```
**Validation:** ✅ Correct methodology
**Test Plan:** Section 8.1 (Grid Convergence Index)

#### ✅ Mass Conservation Checks
```cpp
double computeTotalMass(patch_hierarchy, var_idx);
int checkForNegatives(patch_hierarchy, var_idx, tolerance);
bool checkForNaNInf(patch_hierarchy, var_idx);
```
**Test Plan:** Section 9.1-9.2 (Conservation & diagnostics)

### Test Utilities (TestUtilities.cpp - 401 lines)

#### ✅ Dimensionless Numbers
```cpp
double computeSchmidtNumber(double viscosity, double diffusivity);
double computePecletNumber(double velocity, double length, double diffusivity);
double computeCFL_advection(double velocity, double dt, double dx);
double computeCFL_diffusion(double diffusivity, double dt, double dx);
```
**Test Plan:** Section 5 (Temporal stability & timestep constraints)

---

## 3. Test Coverage vs. Comprehensive Test Plan

### Section Mapping (User's 16-Section Test Plan)

| User Section | Test Plan Item | IBAMR Test | Status | Files |
|--------------|----------------|------------|--------|-------|
| **0** | Practical testing notes | All | ✅ | README.md, utilities |
| **1.1** | Pure diffusion 1D & 2D | Test02 | ⏳ Template | Test02_Diffusion_Analytic/ |
| **1.2** | Pure advection 1D | Test03 | ⏳ Template | Test03_Advection_Analytic/ |
| **1.3** | Combined adv-diff 1D | Test03 | ⏳ Template | (combine with Test03) |
| **1.4** | MMS | Test04 | ⏳ Template | Test04_MMS/ |
| **2** | BC tests (all types) | Test07 | ⏳ Template | Test07_BCs/ |
| **3.1** | Sphere source | Test08 | ⏳ Template | Test08_SphereSource/ |
| **3.2** | Convective plume | Test08 | ⏳ Template | (extend Test08) |
| **4.1** | Prescribed velocity | Test03 | ⏳ Template | Test03_Advection_Analytic/ |
| **4.2** | INS + scalar | Test08 | ⏳ Template | Test08_SphereSource/ |
| **4.3** | IB + scalar BC | Test10 | ⏳ Template | Test10_MovingIB/ |
| **5.1** | CFL advection | Test12 | ⏳ Template | Test12_TimeStep/ |
| **5.2** | Diffusion timestep | Test12 | ⏳ Template | Test12_TimeStep/ |
| **6.1** | High Schmidt Sc 0.7→340 | Test09 | ⏳ Template | Test09_HighSc/ |
| **7.1** | AMR conservation | Test11 | ⏳ Template | Test11_AMR/ |
| **7.2** | Coarse-fine error | Test11 | ⏳ Template | Test11_AMR/ |
| **8.1** | Grid refinement GCI | Test02,04 | ⏳ Template | ErrorCalculator has GCI |
| **8.2** | Time-step convergence | Test12 | ⏳ Template | Test12_TimeStep/ |
| **9.1** | Global mass conservation | Test06 | ⏳ Template | Test06_MassConservation/ |
| **9.2** | Local flux checks | Test06 | ⏳ Template | Test06_MassConservation/ |
| **9.3** | PDF / intermittency | Test14 | ⏳ Template | Test14_Benchmarks/ |
| **10.1** | Rotating cylinder | Test14 | ⏳ Template | Test14_Benchmarks/ |
| **10.2** | 3D sphere/cube | Test14 | ⏳ Template | Test14_Benchmarks/ |
| **10.3** | Undulating body | Test14 | ⏳ Template | Test14_Benchmarks/ |
| **11** | Sensitivity tests | Test11 | ⏳ Template | Test11_AMR/ |
| **12** | Production checklist | All | ⏳ Template | (aggregate) |

### Summary by Tier

#### ✅ **Tier 1: Infrastructure (Test01 only)**
- Test01_SmokeTest: **FULLY IMPLEMENTED** ✅
  - Variable registration ✅
  - Basic I/O ✅
  - Boundary conditions ✅
  - No crashes, no NaNs ✅

#### ⏳ **Tier 2: Unit Operators (Tests 02-06)**
All have templates and analytical solutions ready:

- **Test02 (Diffusion):**
  - ✅ `gaussianDiffusion2D()` implemented
  - ⏳ Main loop needs implementation
  - ⏳ Convergence study loop needed

- **Test03 (Advection):**
  - ✅ `advectedGaussian2D()` implemented
  - ⏳ Uniform velocity field setup needed

- **Test04 (MMS):**
  - ✅ `manufacturedSolution2D()` ✅
  - ✅ `manufacturedSource2D()` ✅
  - ⏳ Source term injection needed

- **Test05 (Discontinuous):**
  - ✅ `topHat2D()` implemented
  - ⏳ Monotonicity checks needed

- **Test06 (Mass Conservation):**
  - ✅ `computeTotalMass()` ✅
  - ⏳ Time-series tracking needed

#### ⏳ **Tier 3: Physical Validation (Tests 07-10)**

- **Test07 (BCs):**
  - ✅ Robin BC framework ready
  - ⏳ Need 3 sub-tests (Dirichlet, Neumann, flux)

- **Test08 (Sphere Source):**
  - ✅ `sphereSource3D()` analytical solution ✅
  - ⏳ INS velocity solver coupling needed
  - 📚 Compare with Lei et al. (10308831)

- **Test09 (High Sc):**
  - ✅ Crank-Nicolson already configured ✅
  - ⏳ Test loop for Sc = [0.7, 100, 340, 1000]

- **Test10 (Moving IB):**
  - ✅ IB framework exists (`IBEELKinematics.cpp`)
  - ⏳ Coupling with scalar field needed

#### ⏳ **Tier 4: Production Readiness (Tests 11-14)**

- **Test11 (AMR):**
  - ✅ SAMRAI AMR infrastructure ready
  - ⏳ Coarse-fine flux matching checks needed

- **Test12 (Time-step):**
  - ✅ CFL computation functions ready
  - ⏳ Temporal convergence study loop

- **Test13 (Long Run):**
  - ✅ Infrastructure ready
  - ⏳ Extended time loop + drift monitoring

- **Test14 (Benchmarks):**
  - ⏳ Rotating cylinder (Yan & Zu comparison)
  - ⏳ 3D sphere (Richter & Nikrityuk)
  - ⏳ Undulating body (arXiv 2408.16136v1)

---

## 4. Critical Implementation Gaps

### 🔴 HIGH PRIORITY

1. **Test02-06 Implementation (Unit Operators)**
   - **Impact:** Cannot verify solver correctness
   - **Effort:** ~2-3 days per test (following Test01 pattern)
   - **Blockers:** None; all utilities ready

2. **Test04 MMS Source Term Injection**
   - **Impact:** Gold standard for verification missing
   - **Effort:** 1 day
   - **Note:** Section 13 of test plan explicitly requires this

3. **Test09 High-Sc Validation**
   - **Impact:** Production use requires Sc=340 (water)
   - **Effort:** 1 day (solver already implicit)
   - **Reference:** arXiv 2408.16136v1 used Sc=340 for water

### 🟡 MEDIUM PRIORITY

4. **Test08 Sphere Source + INS Coupling**
   - **Impact:** Literature validation (Lei et al.)
   - **Effort:** 2-3 days
   - **Reference:** 10308831 (Lei, Crimaldi & Li)

5. **Test10 Moving IB Coupling**
   - **Impact:** Core production feature
   - **Effort:** 2-3 days
   - **Note:** `IBEELKinematics.cpp` exists but needs scalar integration

6. **Test11 AMR Coarse-Fine Conservation**
   - **Impact:** AMR production use
   - **Effort:** 2 days
   - **Note:** Section 7.1 requires <0.1% mass loss during regrid

### 🟢 LOW PRIORITY (Refinement)

7. **Test14 Full Benchmark Suite**
   - Rotating cylinder (Yan & Zu)
   - 3D sphere (Richter & Nikrityuk)
   - Undulating body (arXiv 2408.16136v1)
   - **Effort:** 3-5 days
   - **Note:** Can use coarse-grid validation initially

---

## 5. IBAMR-Specific Recommendations

### ✅ Already Following Best Practices

1. **Implicit Diffusion:** Crank-Nicolson configured ✅
2. **Proper Coupling:** Two-way NS + AdvDiff registration ✅
3. **Conservative Operators:** SAMRAI conservative interpolation ✅
4. **PETSc Solvers:** FGMRES for Helmholtz equations ✅

### 🔧 Suggested Enhancements

1. **Add Slope Limiters for Test05**
   ```
   convective_limiter = "MINMOD"  // or "SUPERBEE"
   ```
   Prevents oscillations in discontinuous transport (Section 1.2 of plan)

2. **AMR Tagging for Scalar Gradients**
   ```cpp
   // In input2d, add:
   tagging_method = "GRADIENT_DETECTOR"
   Gradient_tol = 0.1
   ```
   Ensures AMR refines around scalar fronts (Test11)

3. **Checkpoint/Restart for Long Runs**
   Already configured in Test01 input file ✅

4. **HDF5 Output for Post-Processing**
   ```
   uses_visit = TRUE  // Already set ✅
   ```

---

## 6. Input File Validation (input2d)

### ✅ Correct Settings (Test01_SmokeTest/input2d)

| Parameter | Value | Validation | Test Plan Ref |
|-----------|-------|------------|---------------|
| `diffusion_coefficient` | 0.001 | ✅ Reasonable | Sc = ν/κ |
| `DT` | 0.01 | ✅ Safe CFL | Section 5 |
| `END_TIME` | 1.0 | ✅ Reasonable | - |
| `diffusion_time_stepping_type` | CRANK_NICOLSON | ✅ 2nd-order implicit | Section 5.2 |
| `convective_difference_type` | CENTERED | ✅ 2nd-order | Section 1.1 |
| `helmholtz_solver_db/ksp_type` | fgmres | ✅ Recommended | IBAMR docs |
| `rel_residual_tol` | 1.0e-6 | ✅ Good | Section 14 |

### ⚠️ Considerations for Production

1. **Grid Resolution:**
   - Current: 64×64
   - For high-Sc: Need to resolve `δ ~ √(κT)` (Section 6.1)
   - Recommendation: Start with 128×128 for Test09

2. **Time Step:**
   - Current CFL_adv = 0.01 (very safe)
   - Can increase to ~0.5 for efficiency (Section 5.1)
   - Implicit diffusion removes diffusive CFL limit ✅

3. **AMR Levels:**
   - Current: `max_levels = 1` (no AMR)
   - For Test11: Need `max_levels = 3`

---

## 7. Gap Analysis: Test Plan Section-by-Section

### Section 1: Unit/Operator Tests

| Test Plan Item | Implementation | Status | Gap |
|----------------|----------------|--------|-----|
| 1.1 Pure diffusion 1D & 2D | `gaussianDiffusion2D()` exists | ⏳ | Need convergence loop |
| 1.2 Pure advection 1D | `advectedGaussian2D()` exists | ⏳ | Need uniform velocity setup |
| 1.3 Combined adv-diff 1D | Analytical solution exists | ⏳ | Need test implementation |
| 1.4 MMS | Source term implemented | ⏳ | Need injection mechanism |

**Estimated Effort:** 4-5 days

### Section 2: Boundary Conditions

| BC Type | Framework | Status | Gap |
|---------|-----------|--------|-----|
| Dirichlet | Robin with a=1, b=0 | ✅ | Need sub-test |
| Neumann | Robin with a=0, b=1 | ✅ | Need sub-test |
| No-flux wall | Same as Neumann | ✅ | Need sub-test |
| IB Dirichlet | IB framework exists | ⏳ | Need IB+scalar coupling |

**Estimated Effort:** 2-3 days

### Section 3: Source/Emitter Tests

| Test | Analytical Solution | Status | Gap |
|------|---------------------|--------|-----|
| Steady sphere source | `sphereSource3D()` | ✅ | Need quiescent fluid test |
| Convective plume | Gaussian plume theory | ⏳ | Need INS coupling |

**Estimated Effort:** 2-3 days

### Section 4: Coupling Tests

| Test | Components | Status | Gap |
|------|------------|--------|-----|
| Prescribed velocity | Analytical u(x,y) | ⏳ | Need velocity field injection |
| INS + scalar | Both solvers exist | ⏳ | Need coupling test |
| IB + scalar | IB code exists | ⏳ | Need integration |

**Estimated Effort:** 3-4 days

### Section 5: Temporal Stability

| Test | Utility | Status | Gap |
|------|---------|--------|-----|
| CFL advection | `computeCFL_advection()` | ✅ | Need stability sweep |
| Diffusion timestep | Crank-Nicolson (unconditionally stable) | ✅ | Verify only |

**Estimated Effort:** 1 day

### Section 6: High-Schmidt Regime

| Sc Value | Air/Water | Status | Gap |
|----------|-----------|--------|-----|
| 0.7 | Air | ⏳ | Need test run |
| 340 | Water | ⏳ | Need test run |
| 1000 | Extreme | ⏳ | Need test run |

**Estimated Effort:** 1-2 days (solver already configured)

### Section 7: AMR Tests

| Test | Framework | Status | Gap |
|------|-----------|--------|-----|
| Coarse-fine conservation | SAMRAI AMR | ✅ | Need mass tracking across regrid |
| Error at boundaries | SAMRAI interpolation | ✅ | Need local error computation |

**Estimated Effort:** 2 days

### Section 8: Convergence Studies

| Test | Utility | Status | Gap |
|------|---------|--------|-----|
| Grid refinement (GCI) | `computeConvergenceRate()` | ✅ | Need multi-grid loop |
| Time-step convergence | Same utility | ✅ | Need multi-dt loop |

**Estimated Effort:** 2 days

### Section 9: Conservation & Diagnostics

| Test | Utility | Status | Gap |
|------|---------|--------|-----|
| Global mass | `computeTotalMass()` | ✅ | Need time-series |
| Local flux | Manual patch loop | ⏳ | Need implementation |
| PDF/intermittency | Post-processing | ⏳ | Need statistics module |

**Estimated Effort:** 2-3 days

### Section 10: Benchmark Reproduction

| Benchmark | Reference | Status | Effort |
|-----------|-----------|--------|--------|
| Rotating cylinder | Yan & Zu (via Lei 10308831) | ⏳ | 2 days |
| 3D sphere | Richter (via Lei 10308831) | ⏳ | 2 days |
| Undulating body | arXiv 2408.16136v1 | ⏳ | 3 days |

**Estimated Effort:** 5-7 days

### Section 11-12: Sensitivity & Production Checklist

These aggregate previous tests. Estimated effort depends on completeness of Tests 1-10.

---

## 8. Recommended Implementation Roadmap

### Phase 1: Core Verification (Week 1-2)
**Goal:** Prove solver is mathematically correct

1. **Test02: Pure Diffusion** (2 days)
   - Implement 3-grid convergence study
   - Verify 2nd-order convergence
   - **Pass Criteria:** Convergence rate ≈ 2.0 ± 0.3

2. **Test03: Pure Advection** (1 day)
   - Uniform velocity field
   - Check translation accuracy
   - **Pass Criteria:** Shape preservation

3. **Test04: MMS** (2 days)
   - Inject source term in AdvDiffHierarchyIntegrator
   - Convergence study
   - **Pass Criteria:** Convergence rate ≈ 2.0 ± 0.3

4. **Test06: Mass Conservation** (1 day)
   - Time-series mass tracking
   - **Pass Criteria:** Drift < 1e-10

### Phase 2: Physical Validation (Week 3)
**Goal:** Match literature and handle production scenarios

5. **Test07: Boundary Conditions** (2 days)
   - 3 sub-tests: Dirichlet, Neumann, flux
   - **Pass Criteria:** BC error ≤ 1e-6

6. **Test09: High Schmidt** (2 days)
   - Sc = [0.7, 100, 340, 1000]
   - **Pass Criteria:** Stable, no oscillations

7. **Test08: Sphere Source** (3 days)
   - INS + scalar coupling
   - Compare with Lei et al. (10308831)
   - **Pass Criteria:** Match ±10%

### Phase 3: Production Readiness (Week 4)
**Goal:** Ensure robustness for long runs and AMR

8. **Test11: AMR** (2 days)
   - Multi-level tests
   - Conservation across regridding
   - **Pass Criteria:** Mass loss < 0.1%

9. **Test12: Time-step Convergence** (1 day)
   - Temporal order verification
   - **Pass Criteria:** Convergence rate ≈ 2.0

10. **Test13: Long Run** (1 day)
    - Extended integration (T = 50-100)
    - **Pass Criteria:** No secular drift

### Phase 4: Full Validation (Week 5+)
**Goal:** Complete benchmark comparisons

11. **Test10: Moving IB** (3 days)
    - Couple with `IBEELKinematics.cpp`
    - **Pass Criteria:** No instabilities

12. **Test14: Benchmarks** (5-7 days)
    - Rotating cylinder
    - 3D sphere
    - Undulating body (arXiv 2408.16136v1)
    - **Pass Criteria:** Match published data ±10-20%

---

## 9. Key Findings

### ✅ Strengths

1. **Excellent Framework:**
   - Well-architected test suite
   - Comprehensive common utilities
   - Reference implementation (Test01) is high-quality

2. **Correct IBAMR Usage:**
   - Proper API calls
   - Conservative schemes
   - Implicit diffusion for high-Sc

3. **Complete Analytical Solutions:**
   - All major test cases have exact solutions
   - MMS framework ready
   - Error calculators implemented

### ⚠️ Weaknesses

1. **Implementation Incomplete:**
   - Only 1 of 14 tests fully implemented (7%)
   - Core verification tests (02-06) are critical gaps

2. **No Benchmark Validation Yet:**
   - Cannot verify against Lei et al. (10308831)
   - Cannot verify against arXiv 2408.16136v1

3. **Production Scenarios Untested:**
   - High-Sc (Sc=340) not validated
   - AMR not tested
   - Moving IB + scalar uncoupled

---

## 10. Production Readiness Assessment

### Current Status: 🔴 NOT PRODUCTION-READY

**Reasons:**
1. Core verification tests incomplete (Tests 02-06)
2. No literature validation (Test08, Test14)
3. High-Sc regime untested (Test09)
4. AMR untested (Test11)

### Path to Production:

**Minimum Viable Tests (MVP):**
- ✅ Test01 (Smoke) — DONE
- ⏳ Test02 (Diffusion) — CRITICAL
- ⏳ Test04 (MMS) — CRITICAL
- ⏳ Test06 (Mass Conservation) — CRITICAL
- ⏳ Test09 (High-Sc) — CRITICAL for water simulations

**Estimated Time to MVP:** 2 weeks (assuming full-time development)

**Full Production Readiness:** 4-6 weeks

---

## 11. Checklist for User's Comprehensive Test Plan

Mapping to **Section 12 — Production Checklist** from test plan:

- [ ] Unit operator tests: diffusion, advection, advection–diffusion (analytic) — **Test02, Test03**
- [ ] MMS test: manufactured solution global order verified — **Test04**
- [ ] BC tests: Dirichlet, Neumann, no-flux, IB Dirichlet/Neumann — **Test07**
- [ ] Source tests: steady sphere source, plume-like experiments — **Test08**
- [ ] Coupling: scalar advected by analytic velocity and by computed INS velocity — **Test03, Test08**
- [ ] IB coupling: moving IB with scalar BCs — **Test10**
- [ ] High-Sc tests run: Sc up to production target (e.g., 340) — **Test09**
- [ ] AMR tests: coarse-fine conservation & error — **Test11**
- [ ] Grid/time convergence & GCI — **Test02, Test04, Test12**
- [ ] Mass conservation & flux checks — **Test06**
- [ ] Benchmark reproductions: rotating cylinder + sphere/cube + undulating-foil trends — **Test14**

**Current Progress:** 1/11 complete (9%)

---

## 12. Immediate Next Steps

### For Next Commit:

1. **Implement Test02 (Pure Diffusion)**
   - Copy Test01 structure
   - Replace initial condition with `gaussianDiffusion2D()`
   - Add convergence study loop (3 grids)
   - Compute L2 error and convergence rate
   - **Expected:** 2-3 days

2. **Implement Test04 (MMS)**
   - Add `manufacturedSource2D()` to RHS
   - Convergence study
   - **Expected:** 2 days

3. **Implement Test06 (Mass Conservation)**
   - Time-series mass tracking
   - **Expected:** 1 day

### For Production Use:

Complete Minimum Viable Tests (MVP) list above.

---

## 13. References Validation

### ✅ Papers Cited in Test Plan

1. **Lei et al. (2021)** — 10308831
   - Rotating cylinder validation (Test14)
   - Sphere/cube validation (Test14)
   - PDF/whiff statistics (Test14)

2. **Kamran et al. (2024)** — arXiv 2408.16136v1
   - Undulating body (Test14)
   - High-Sc analysis (Test09)
   - Sc = 0.7 (air) and Sc = 340 (water)

3. **Yan & Zu** — (via Lei et al.)
   - Rotating cylinder benchmark

4. **Richter & Nikrityuk** — (via Lei et al.)
   - 3D sphere heat transfer analog

### Test Plan Alignment:

User's test plan explicitly references these papers and suggests:
- "reproduce their validation figures"
- "same parameters as in arXiv 2408.16136"
- "compare scalar contour/line plots to Yan & Zu"

**Current Status:** Framework ready, implementations pending.

---

## 14. Conclusion

The IBAMR implementation has:

✅ **Solid foundation:**
- Correct IBAMR syntax and API usage
- Complete analytical solution library
- Robust error calculation tools
- Proper Crank-Nicolson implicit diffusion

⏳ **Critical gaps:**
- Only 7% of tests implemented (1/14)
- Core verification missing (Tests 02-06)
- No benchmark validation
- Production scenarios untested

**Recommendation:**
Follow the **4-week roadmap** in Section 8 to achieve production readiness. Prioritize Tests 02, 04, 06, 09 for MVP (Minimum Viable Product).

---

**Review completed by:** Claude Code
**Date:** 2025-11-17
**Status:** Framework ✅ | Implementation ⏳ | Production ❌
