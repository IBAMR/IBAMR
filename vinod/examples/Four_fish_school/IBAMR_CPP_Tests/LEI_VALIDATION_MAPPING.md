# Lei et al. (2021) Validation Test Mapping

## Document Purpose

This document provides a **precise mapping** between the validation tests extracted from **Lei, Crimaldi & Li (2021) "Navigation in odor plumes: How do the flapping kinematics modulate the odor landscape"** (AIAA 2021-2817) and their IBAMR C++ implementation.

**Paper Location:** `Navigation in odor plumes How do the flapping kinematics modulate the odor landscape.pdf`

---

## Complete Validation Suite Extracted from Lei et al. (2021)

The paper contains **THREE distinct validation cases** in Section II.B "Solver validation":

1. **TEST 08** (proposed) → 2-D Rotating Cylinder (Yan & Zu 2008)
2. **TEST 09** (proposed) → 3-D Sphere & Cube (Richter & Nikrityuk 2012)
3. **TEST 10** (proposed) → Full Pitch-Plunge Odor Problem (Lei et al. 2021)

---

## IBAMR Implementation Mapping

### ✅ TEST 15: Rotating Cylinder Validation

**Maps to:** Lei et al. (2021) Section II.B, Figure 2

**Reference:** Yan & Zu (2008) - Lattice Boltzmann Method
- *"Numerical simulation of heat transfer and fluid flow past a rotating isothermal cylinder – A LBM approach"*
- International Journal of Heat and Mass Transfer
- DOI: 10.1016/j.ijheatmasstransfer.2007.07.053

**Parameters:**
| Parameter | Value | Source |
|-----------|-------|--------|
| Re (Reynolds) | 200 | Lei et al. Table/Text |
| Pr (Prandtl) | 0.5 | Lei et al. Figure 2 caption |
| k (rotation) | 0.5 | V/U ratio |
| Geometry | 2D circular cylinder, D=1.0 | Paper Section II.B |
| Grid | 672 × 416 (from main simulation) | Lei et al. Figure 4 |
| BC (scalar) | Dirichlet on cylinder (isothermal) | Yan & Zu methodology |
| BC (velocity) | Rotating no-slip, u_θ = V | Yan & Zu methodology |

**Expected Output:**
- Streamlines matching Yan & Zu Figure 2
- Temperature (scalar) contours matching published shape
- Recirculation bubble ~1.5-2.0 D
- Thermal boundary layer ~0.1 D

**Implementation:**
- Directory: `IBAMR_CPP_Tests/Test15_RotatingCylinder/`
- Executable: `test15_rotating_cylinder`
- Input file: `input2d`
- Status: ✅ Template created, ready for execution

---

### ✅ TEST 16: 3D Sphere & Cube Validation

**Maps to:** Lei et al. (2021) Section II.B, Figure 3

**Reference:** Richter & Nikrityuk (2012) - CFD benchmark
- *"Drag forces and heat transfer coefficients for spherical, cuboidal and ellipsoidal particles in cross flow at sub-critical Reynolds numbers"*
- International Journal of Heat and Mass Transfer, 55(4), 1343-1354
- DOI: 10.1016/j.ijheatmasstransfer.2011.09.005

**Parameters:**
| Parameter | Value | Source |
|-----------|-------|--------|
| Re (Reynolds) | 200 | Lei et al. Section II.B |
| Pr (Prandtl) | 0.744 | Lei et al. Figure 3 caption |
| Sc (Schmidt) | 0.744 | IBAMR notation (same as Pr) |
| Geometry (Case 1) | 3D sphere, D=1.0 | Richter & Nikrityuk |
| Geometry (Case 2) | 3D cube, L=1.0 | Richter & Nikrityuk |
| Grid | 352 × 224 × 160 | Lei et al. Figure 4 |
| BC (scalar) | Dirichlet on surface (C = C_h) | Richter & Nikrityuk |
| BC (velocity) | No-slip on surface | Standard |

**Expected Output:**
- Scalar contours on symmetry plane matching Figure 3
- Thermal boundary layer thickness matching reference
- Wake scalar spreading matching reference
- Nusselt number Nu ≈ 15-20 (sphere), 18-23 (cube)

**Implementation:**
- Directory: `IBAMR_CPP_Tests/Test16_3DSphere/`
- Executable: `test16_sphere` (sphere case), `test16_cube` (cube case - future)
- Input file: `input3d_sphere`, `input3d_cube`
- Status: ✅ Template created, sphere case ready
- Note: **3D simulation** - high computational cost

---

### ✅ TEST 17: Pitch-Plunge Full Case

**Maps to:** Lei et al. (2021) Section II.C, Figures 4-10

**This is the MAIN PHYSICS CASE** from the paper, not just solver validation.

**Reference:** Lei et al. (2021) - Original work
- Full odor plume navigation problem
- Upstream sphere source + downstream flapping airfoil
- Vortex-odor interaction

**Parameters:**
| Parameter | Value | Source |
|-----------|-------|--------|
| Re (Reynolds) | 200 | Lei et al. Table 1 |
| Sc (Schmidt) | 0.71 | Lei et al. Table 1 |
| St (Strouhal) | 0.9 | Lei et al. Table 1 |
| **Sphere Source** | | |
| - Diameter | D = L | Section II.C |
| - Location | (-3L, 0) | Figure 4 |
| - BC | C = C_h (Dirichlet) | Section II.C |
| **Airfoil (2D)** | | |
| - Semi-major axis | a_x = L | Section II.C |
| - Semi-minor axis | a_y = 0.24L | Section II.C |
| - Plunge | y(t) = (L/2)sin(2πft) | Eq. in paper |
| - Pitch | θ(t) = (1/6)cos(2πft) | Eq. in paper (≈9.5°) |
| - Scalar BC | ∂C/∂n = 0 (Neumann) | Section II.C |
| **Airfoil (3D)** | | |
| - Spanwise | a_z = 2L (from L/0.5) | Section II.C |
| **Grid** | | |
| - 2D | 672 × 416 | Figure 4 |
| - 3D | 352 × 224 × 160 | Figure 4 |
| - Refinement | Two nested regions (sphere + airfoil) | Figure 4 |

**Expected Output (Figures 6-10):**

**Figure 6 - 2D Vorticity & Odor:**
- Inverse von Kármán vortex street
- Vortex spacing ~0.5-1.0 L
- Odor trapped in vortex cores
- Asymmetric wake structure

**Figures 7-8 - Odor PDF:**
- Nondimensional C* = (C - C_l)/(C_h - C_l)
- PDF peak at C* ≈ 0 (background)
- PDF trough at C* ≈ 0.2 (wake dilution)
- Tail to C* = 1 (near source)

**Figures 9-10 - 3D Structure:**
- Two branches of horseshoe vortices
- Bifurcation angle ≈ 37°
- Odor iso-surface (C* = 0.05) follows vortex rings
- Faster dissipation than 2D

**Implementation:**
- Directory: `IBAMR_CPP_Tests/Test17_PitchPlunge/`
- Executable: `test17_pitch_plunge_2d` (2D), `test17_pitch_plunge_3d` (3D)
- Input files: `input2d`, `input3d`
- IB meshes: Sphere (.vertex) + ellipse/ellipsoid (.vertex)
- Status: ✅ Template created, kinematics defined

---

## Parameter Summary Table

Comparison of all three validation tests:

| Test | Reference | Re | Pr/Sc | Geometry | Dimension | Grid | Difficulty |
|------|-----------|----|----|----------|-----------|------|-----------|
| **Test15** | Yan & Zu 2008 | 200 | Pr=0.5 | Rotating cylinder D=1 | 2D | 672×416 | Advanced |
| **Test16** | Richter & Nikrityuk 2012 | 200 | Pr=0.744 | Sphere/Cube D=1 | **3D** | 352×224×160 | Advanced |
| **Test17** | Lei et al. 2021 | 200 | **Sc=0.71** | Sphere + ellipsoid | 2D/3D | 672×416 (2D) | Expert |

**Key Observations:**
- All use **Re = 200** (consistent laminar flow regime)
- Test15-16 use **Prandtl number** (heat transfer, air-like)
- Test17 uses **Schmidt number** (mass transfer, odor)
- Test17 is **time-dependent** (pitch-plunge motion)
- Test16-17 have **3D versions** (computationally expensive)

---

## Validation Workflow

### Phase 1: Solver Validation (Test15-16)
1. **Test15**: Rotating cylinder → validates scalar transport around curved, rotating boundary
2. **Test16**: 3D sphere/cube → validates 3D extension, complex geometries

### Phase 2: Production Case (Test17)
3. **Test17**: Full pitch-plunge → validates complete physics (source + IB + kinematics)

**Recommended Order:**
1. Run Test15 (2D, fast, ~1 hour)
2. Run Test16 sphere (3D, medium, ~4 hours with 16 cores)
3. Run Test17 2D (2D, medium, ~4 hours)
4. (Optional) Run Test16 cube if sphere passes
5. (Optional) Run Test17 3D if 2D passes (~48 hours)

---

## Acceptance Criteria

### Test15 Acceptance
- ✅ Steady state reached (|∂C/∂t| < 1e-6)
- ✅ No NaN/Inf values
- ✅ Scalar bounded 0 ≤ C ≤ 1
- ✅ Streamlines visually match Yan & Zu Figure 2
- ✅ Recirculation bubble length within ±20%
- ✅ Thermal BL thickness within ±20%

### Test16 Acceptance
- ✅ 3D simulation completes without divergence
- ✅ Scalar field stable (no oscillations)
- ✅ Symmetry plane contours match Richter & Nikrityuk Figure 3
- ✅ Nusselt number within ±15% of reference
- ✅ Wake structure reasonable (length, width)

### Test17 Acceptance
- ✅ Time-periodic steady state reached (multiple flapping cycles)
- ✅ Inverse von Kármán vortex street visible
- ✅ Odor plume from sphere source
- ✅ Odor modulation by wake vortices
- ✅ **PDF shape:** Peak at C*≈0, trough at C*≈0.2 (±0.05)
- ✅ (3D only) Horseshoe vortex bifurcation angle ≈ 37° (±10°)

---

## Differences from Original Paper

### Grid Resolution
**Paper:** Adaptive refinement strategy not fully specified
**IBAMR:** Uses AMR with 2-3 levels, refinement boxes around IB and wake

### Rotation Implementation (Test15)
**Paper:** Yan & Zu uses LBM with rotating boundary
**IBAMR:** Uses immersed boundary with prescribed velocity

### 3D Geometry (Test16)
**Paper:** Richter & Nikrityuk uses body-fitted mesh
**IBAMR:** Uses immersed boundary method (Cartesian grid)

### IB Scalar BC (Test17)
**Paper:** Airfoil has zero-gradient ∂C/∂n = 0
**IBAMR:** Implemented via Neumann BC or feedback forcing

---

## Reference Data Availability

| Test | Reference | Data Availability |
|------|-----------|------------------|
| **Test15** | Yan & Zu 2008 | Figures only (no digital data) |
| **Test16** | Richter & Nikrityuk 2012 | Figures + Nusselt correlations |
| **Test17** | Lei et al. 2021 | Figures 6-10 (PDF provided) |

**Validation Strategy:**
- **Qualitative:** Visual comparison with published figures
- **Quantitative:** Global metrics (Nu, drag, wake length, PDF shape)
- **Acceptance:** ±10-20% error acceptable for coarse grids

---

## Connection to Existing Tests

These **three new tests (15-17)** complement the existing suite:

| Existing Test | Focus | Relation to Test15-17 |
|---------------|-------|----------------------|
| Test08 | Sphere source (steady) | Simpler version of Test16 (2D vs 3D) |
| Test09 | High Schmidt number | Validates Sc parameter used in Test17 |
| Test10 | Moving IB + scalar | Framework for Test17 pitch-plunge |
| Test14 | Benchmarks | Collection test, Test15-17 are detailed |

**Test15-17 provide:**
- Direct literature comparison (Yan & Zu, Richter & Nikrityuk, Lei et al.)
- 3D validation (Test16)
- Full production case (Test17)
- Explicit parameter extraction from Lei et al. (2021)

---

## Building and Running

### Build Literature Validation Tests Only
```bash
cd IBAMR_CPP_Tests
mkdir build && cd build
cmake .. -DBUILD_LITERATURE_ONLY=ON
make -j4
```

### Build All Tests (Including 15-17)
```bash
cmake .. -DBUILD_ALL_TESTS=ON
make -j8
```

### Run Tests
```bash
# Test15 (2D, ~1 hour, 4 cores)
cd Test15_RotatingCylinder/build
mpirun -np 4 ./test15_rotating_cylinder ../input2d

# Test16 (3D, ~4 hours, 16 cores)
cd Test16_3DSphere/build
mpirun -np 16 ./test16_sphere ../input3d_sphere

# Test17 (2D, ~4 hours, 8 cores)
cd Test17_PitchPlunge/build
mpirun -np 8 ./test17_pitch_plunge_2d ../input2d
```

---

## Verification Checklist

Before claiming "Lei et al. (2021) validated":

### Test15 (Rotating Cylinder)
- [ ] Executable compiles
- [ ] Simulation runs to completion
- [ ] Streamlines match Yan & Zu Figure 2
- [ ] Scalar contours match reference
- [ ] Documented in test results

### Test16 (3D Sphere)
- [ ] 3D executable compiles
- [ ] Simulation runs (may take hours)
- [ ] Symmetry plane matches Richter & Nikrityuk Figure 3
- [ ] Nusselt number computed and within ±15%
- [ ] Documented in test results

### Test17 (Pitch-Plunge)
- [ ] 2D executable compiles
- [ ] Kinematics correctly implemented
- [ ] Simulation reaches periodic steady state
- [ ] Vortex street matches Lei et al. Figure 6
- [ ] Odor PDF matches Figures 7-8 (peak/trough)
- [ ] (Optional) 3D version matches Figures 9-10
- [ ] Documented in test results

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-11-17 | AI Assistant | Initial extraction from Lei et al. (2021) |
| | | | Tests 15-17 created and mapped |

---

## References

### Primary Paper
**Lei, H., Crimaldi, J. P., & Li, C. (2021)**
*Navigation in odor plumes: How do the flapping kinematics modulate the odor landscape*
AIAA Aviation 2021 Forum, Paper 2021-2817
DOI: 10.2514/6.2021-2817

### Test15 Reference
**Yan, Y. Y., & Zu, Y. Q. (2008)**
*Numerical simulation of heat transfer and fluid flow past a rotating isothermal cylinder – A LBM approach*
International Journal of Heat and Mass Transfer
DOI: 10.1016/j.ijheatmasstransfer.2007.07.053

### Test16 Reference
**Richter, A., & Nikrityuk, P. A. (2012)**
*Drag forces and heat transfer coefficients for spherical, cuboidal and ellipsoidal particles in cross flow at sub-critical Reynolds numbers*
International Journal of Heat and Mass Transfer, 55(4), 1343-1354
DOI: 10.1016/j.ijheatmasstransfer.2011.09.005

---

**Status:** ✅ Complete parameter extraction and test mapping
**Next Steps:** Execute tests, compare results, document validation
