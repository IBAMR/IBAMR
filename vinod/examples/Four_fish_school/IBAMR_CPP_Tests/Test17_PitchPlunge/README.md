# Test 17: Pitch-Plunge Odor Plume Navigation (Lei et al. 2021)

## Purpose

Validates the **full production case** from Lei et al. (2021) - the canonical odor plume navigation problem with upstream sphere source and downstream flapping (pitch-plunge) ellipsoidal airfoil.

This test corresponds to **Lei et al. (2021) Section II.C, Figures 4-10** - the main physics case demonstrating how flapping kinematics modulate the odor landscape.

## Physical Problem

### Overview

This is a **complete physics simulation** combining:
1. **Odor source**: Upstream sphere emitting odor at concentration C_h
2. **Flapping body**: Downstream ellipsoidal airfoil performing pitch-plunge motion
3. **Wake-odor interaction**: Vortices transport and mix odor
4. **Chemotactic navigation**: (future) Fish responds to odor gradients

### Geometry

**Upstream Sphere (Odor Source):**
- Diameter: **D = L** (same as airfoil characteristic length)
- Center: (-3L, 0) or (-3L, 0, 0) for 3D
- Odor concentration: **C = C_h** (Dirichlet boundary condition)

**Downstream Ellipsoidal Airfoil:**

**2D Case (Ellipse):**
```
(x/a_x)² + (y/a_y)² = 1
```
- Semi-major axis: **a_x = L**
- Semi-minor axis: **a_y = 0.24L**
- Center at origin (when at rest)

**3D Case (Ellipsoid):**
```
(x/a_x)² + (y/a_y)² + (z/a_z)² = 1
```
- **a_x = L** (streamwise)
- **a_y = 0.24L** (lateral)
- **a_z = 2L** (spanwise, from L/0.5 in paper)

**Spacing:**
- Sphere-to-airfoil center: **3L** (upstream source)

### Motion Kinematics

The airfoil performs **pitch-plunge** motion:

**Plunge (heave):**
```
y(t) = (L/2) sin(2πft)
```

**Pitch (angle of attack):**
```
θ(t) = (1/6) cos(2πft)    [radians ≈ 9.5°]
```

**Nondimensional time:**
```
t* = tU/L
```

where:
- **f** = flapping frequency
- **U** = free-stream velocity
- **L** = characteristic length (airfoil semi-major axis)

**Phase relationship:** Pitch leads plunge by 90° (cos vs sin)

### Governing Equations

**Momentum (Navier-Stokes):**
```
∂u/∂t + u·∇u = -∇p + (1/Re)∇²u
∇·u = 0
```

**Scalar Transport (Odor):**
```
∂C/∂t + u·∇C = (1/Sc)∇²C
```

where:
- **Re** = Reynolds number = UL/ν
- **Sc** = Schmidt number = ν/κ (odor diffusivity)

### Dimensionless Parameters

From **Lei et al. (2021) Table 1, Section II.C**:

| Parameter | Symbol | Value | Definition |
|-----------|--------|-------|------------|
| **Reynolds number** | Re | 200 | Re = UL/ν |
| **Schmidt number** | Sc | 0.71 | Sc = ν/κ |
| **Strouhal number** | St | 0.9 | St = fL/U |
| **Plunge amplitude** | h_0 | L/2 | Max displacement |
| **Pitch amplitude** | θ_0 | 1/6 rad | ≈ 9.5° |

**Derived quantities:**
- ν = UL/Re = UL/200
- κ = ν/Sc = (UL/200)/0.71 ≈ 0.007 UL
- f = StU/L = 0.9U/L

**Typical values (if U=1, L=1):**
- ν = 0.005
- κ = 0.007
- f = 0.9 Hz

## Boundary Conditions

### Velocity Field

| Boundary | Condition | Description |
|----------|-----------|-------------|
| **Inlet** | u = (U, 0, 0) | Uniform flow |
| **Sphere surface** | No-slip | u = 0 |
| **Airfoil surface** | No-slip + prescribed motion | u = u_IB(t) |
| **Outlet** | ∂u/∂n = 0 | Neumann (convective) |
| **Lateral** | ∂u/∂n = 0 | Neumann |

### Scalar Field (Odor)

| Boundary | Condition | Description |
|----------|-----------|-------------|
| **Inlet** | C = C_l | Background odor = 0 |
| **Sphere** | C = C_h | Odor source = 1 |
| **Airfoil** | ∂C/∂n = 0 | **Zero-gradient (Neumann)** |
| **Outlet** | ∂C/∂n = 0 | Neumann |
| **Lateral** | ∂C/∂n = 0 | Neumann |

**Critical:** Airfoil has **zero-gradient BC** for odor (not Dirichlet). This allows odor to flow past the fish without being absorbed or emitted.

## Grid Resolution

From **Lei et al. (2021)**:

**2D Grid:**
- **672 × 416** nonuniform grid
- Two high-resolution regions:
  1. Dense grid around sphere source
  2. Dense grid around airfoil and wake

**3D Grid:**
- **352 × 224 × 160**
- Similar nested refinement strategy

**IBAMR AMR Strategy:**
- Level 0 (coarse): Δx ≈ L/32
- Level 1 (medium): Δx ≈ L/128 around sphere and airfoil
- Level 2 (fine): Δx ≈ L/256 in wake and near boundaries
- Adaptive refinement tracks vortices and odor gradients

## Expected Output

This test validates against **Figures 6-10** in Lei et al. (2021).

### Figure 6: Inverse von Kármán Vortex Street (2D)

✅ **Vorticity field ω = ∇×u:**
- Alternating vortex street in wake
- Vortex spacing ~ 0.5-1.0 L
- Asymmetric due to pitch-plunge motion
- Clear dipole structure

✅ **Odor field C(x,y,t):**
- Plume from sphere source
- Modulation by wake vortices
- Odor trapped in vortex cores
- Downstream mixing and dispersion

### Figure 7-8: Odor Probability Distribution Functions

**Nondimensional odor:**
```
C* = (C - C_l) / (C_h - C_l)
```

✅ **PDF shape:**
- Peak at C* ≈ 0 (background)
- Trough at C* ≈ 0.2 (wake dilution)
- Tail extending to C* = 1 (near source)

✅ **Spatial variation:**
- Upstream of airfoil: Gaussian-like
- In wake: Bimodal (vortex cores vs ambient)
- Far downstream: Diffusive tail

### Figure 9-10: 3D Horseshoe Vortices

✅ **Vortex structure:**
- Two branches of horseshoe vortices
- Bifurcation angle: **≈ 37°** (from paper)
- Vortex rings shed periodically

✅ **Odor iso-surfaces:**
- C* = 0.05 iso-surface follows vortex rings
- 3D structure more complex than 2D
- Faster odor dissipation due to 3D mixing

## Pass/Fail Criteria

The test **PASSES** if:

### ✅ **Stability & Convergence:**
- [ ] No NaN/Inf values
- [ ] Simulation completes full period (t* = 0 to t* = 5+)
- [ ] No solver divergence
- [ ] Time-periodic steady state reached

### ✅ **Flow Field (Vorticity):**
- [ ] Inverse von Kármán street visible
- [ ] Vortex spacing ~ 0.5-1.0 L
- [ ] Alternating sign vortices
- [ ] Wake extends 5-10 L downstream

### ✅ **Odor Field:**
- [ ] Plume from sphere source
- [ ] Odor modulation by wake
- [ ] Bounded: C_l ≤ C ≤ C_h
- [ ] No negative concentrations

### ✅ **Quantitative Metrics (2D):**
- [ ] Vortex shedding frequency ≈ f (St = 0.9)
- [ ] PDF peak at C* ≈ 0
- [ ] PDF trough at C* ≈ 0.15-0.25
- [ ] Mean odor <C*> in wake ≈ 0.1-0.3

### ✅ **Quantitative Metrics (3D):**
- [ ] Horseshoe vortex bifurcation angle ≈ 37° (±10°)
- [ ] Two distinct vortex branches
- [ ] Odor dissipates faster than 2D

## Building and Running

```bash
cd Test17_PitchPlunge
mkdir build && cd build
cmake ..
make

# 2D case
./test17_pitch_plunge_2d ../input2d

# 3D case (very expensive!)
mpirun -np 32 ./test17_pitch_plunge_3d ../input3d
```

## Expected Runtime

**2D Case:**
- Coarse grid: ~1 hour (4 cores)
- Medium grid: ~4 hours (8 cores)
- Fine grid: ~12 hours (16 cores)

**3D Case:**
- Coarse grid: ~8 hours (32 cores)
- Fine grid: ~48 hours (64 cores)

**Recommendation:** Start with 2D, then 3D if needed.

## Output Files

```
Test17_PitchPlunge/
├── viz_test17_2d/           ← 2D visualization
│   ├── dumps.visit
│   ├── lag_data.*.vtu       ← IB meshes (sphere + airfoil)
│   └── visit*.*.vtk         ← Vorticity + odor fields
├── viz_test17_3d/           ← 3D visualization
├── test17_odor_pdf.dat      ← Odor PDF histogram
├── test17_vorticity_wake.dat ← Wake vorticity profile
├── test17_centerline.dat    ← Odor along centerline
└── test17_results.txt       ← Validation verdict
```

## Visualization

### With VisIt/ParaView (2D):

1. **Vorticity field ω:**
   - Pseudocolor plot
   - Range: [-20, 20] (adjust as needed)
   - Overlay IB boundaries

2. **Odor field C:**
   - Pseudocolor with transparency
   - Range: [0, 1]
   - Show iso-contours at C = 0.1, 0.5, 0.9

3. **Combined plot:**
   - Vorticity as background
   - Odor contours overlaid
   - Streamlines showing flow

4. **Animation:**
   - Play through time steps
   - Observe vortex shedding
   - Track odor transport

### With VisIt/ParaView (3D):

1. **Volume rendering:**
   - Odor field with transparency
   - Color by C*

2. **Iso-surfaces:**
   - C* = 0.05 (vortex-odor structure)
   - Color by vorticity magnitude

3. **Slice plots:**
   - y=0 plane (compare with 2D)
   - z=0 plane (spanwise symmetry)

## Implementation Notes

### IB Meshes Required

**Sphere (.vertex file):**
- 2D: Circle with ~128 points
- 3D: Icosphere with ~2000-5000 triangles

**Ellipsoidal Airfoil (.vertex file):**
- 2D: Ellipse with ~256 points
- 3D: Ellipsoid with ~5000-10000 triangles

### Kinematics Prescription

```cpp
// In IBAMR, prescribe IB motion as:
void getKinematics(double t, double& y_center, double& theta) {
    double f = St * U / L;          // Frequency
    y_center = 0.5 * L * sin(2.0 * M_PI * f * t);
    theta = (1.0/6.0) * cos(2.0 * M_PI * f * t);
}
```

### Zero-Gradient BC on Airfoil

Critical for fish-odor interaction:
- Fish does not emit odor
- Fish does not absorb odor
- Odor flows around fish boundary layer
- Implemented via Neumann BC: ∂C/∂n = 0

### Post-Processing

**Odor PDF:**
```python
# Sample odor field in wake region
C_samples = sample_field_in_box(x=[2L, 5L], y=[-L, L])
C_star = (C_samples - C_l) / (C_h - C_l)
hist, bins = np.histogram(C_star, bins=50, range=[0, 1])
pdf = hist / hist.sum()
```

**Vortex bifurcation angle (3D):**
```python
# Identify vortex cores using Q-criterion
Q = 0.5 * (||Omega||^2 - ||S||^2)  # Omega=vorticity, S=strain
vortex_cores = Q > threshold

# Fit lines to vortex branches
angle = arctan(Δy / Δx)
```

## Physical Interpretation

### Why This Test Matters

1. **Complete physics case**: Combines all elements (source, flow, transport)
2. **Wake-odor interaction**: Validates coupling between vortices and scalar
3. **Navigation context**: Foundation for chemotactic fish behavior
4. **Benchmark**: Reproduces published results

### Key Physics Insights

**2D vs 3D:**
- **2D**: Stronger vortices, slower odor dissipation
- **3D**: Weaker vortices, faster mixing, more realistic

**Strouhal Number (St = 0.9):**
- **St < 0.2**: Poor propulsion, weak wake
- **St ≈ 0.3**: Optimal propulsion (fish, birds)
- **St = 0.9**: Strong wake, enhanced mixing (used here for odor transport study)

**Schmidt Number (Sc = 0.71):**
- Air-like (Sc ~ 0.7)
- Water would be Sc ~ 340 (Test09 validates high Sc)

## References

### Primary Reference

**Lei, H., Crimaldi, J. P., & Li, C. (2021)**
*"Navigation in odor plumes: How do the flapping kinematics modulate the odor landscape"*
AIAA Aviation 2021 Forum, Paper 2021-2817
DOI: 10.2514/6.2021-2817

**Key sections:**
- **Section II.C**: Problem setup (Table 1, parameters)
- **Figure 4**: Computational domain
- **Figure 5**: Kinematics definition
- **Figures 6-10**: Results to reproduce

### Related Papers

**Kamran et al. (2024):**
*"How does vortex dynamics help undulating bodies spread odor"*
arXiv:2408.16136v1
- Extension to undulating body kinematics
- Higher Sc (water-like)

### IBAMR Examples

- `IBAMR/examples/IB/explicit/ex2` - Oscillating cylinder
- `IBAMR/examples/IB/explicit/ex3` - Flapping airfoil
- `IBAMR/examples/adv_diff/ex1` - Passive scalar with IB

## Related Tests

- **Test15_RotatingCylinder**: 2D cylinder validation
- **Test16_3DSphere**: 3D validation
- **Test08_SphereSource**: Sphere source only (no flapping)
- **Test10_MovingIB**: IB-scalar coupling (simpler case)

## Extensions & Future Work

1. **Closed-loop chemotaxis:**
   - Fish kinematics = f(∇C, ∂C/∂t)
   - Gradient-following behavior

2. **Multiple fish:**
   - School formation effects
   - Odor trail following

3. **Unsteady source:**
   - Pulsed odor emission
   - Moving prey/predator

4. **Environmental realism:**
   - Background turbulence
   - Stratification
   - Obstacles

## Status

- [x] Test specification created
- [ ] Directory structure created
- [ ] main.cpp implemented (2D)
- [ ] main.cpp implemented (3D)
- [ ] input2d configured
- [ ] input3d configured
- [ ] Sphere vertex file
- [ ] Ellipse/ellipsoid vertex file
- [ ] Kinematics functions
- [ ] Post-processing scripts
- [ ] Test executed and validated
- [ ] Results documented

---

**Test ID:** Test17
**Type:** Full Production Case (Lei et al. 2021)
**Category:** Odor Plume Navigation - Pitch-Plunge Kinematics
**Difficulty:** Expert (moving IB + source + scalar + validation)
**Last Updated:** 2025-11-17
