# Test 16: 3D Sphere & Cube Validation (Richter & Nikrityuk 2012)

## Purpose

Validates 3-D advection-diffusion solver against published CFD benchmark from **Richter & Nikrityuk (2012)** for scalar transport around stationary sphere and cube.

This test corresponds to **Lei et al. (2021) Section II.B, Figure 3** - "Solver validation" using 3D immersed boundaries with scalar transport.

## Physical Problem

### Geometry

**Case 1: Sphere**
- Diameter: **D = 1.0**
- Center: (0, 0, 0)
- Stationary

**Case 2: Cube**
- Edge length: **L = 1.0**
- Center: (0, 0, 0)
- Stationary

### Governing Equations

**Momentum (Navier-Stokes):**
```
∂u/∂t + u·∇u = -∇p + (1/Re)∇²u
∇·u = 0
```

**Scalar Transport (Heat/Odor):**
```
∂C/∂t + u·∇C = α∇²C
```

where:
- α = ν/Pr (thermal diffusivity in Richter & Nikrityuk)
- α = ν/Sc (scalar diffusivity for odor in IBAMR context)

### Dimensionless Parameters

From **Richter & Nikrityuk (2012)** and **Lei et al. (2021)**:

| Parameter | Symbol | Value | Definition |
|-----------|--------|-------|------------|
| **Reynolds number** | Re | 200 | Re = UD/ν |
| **Prandtl number** | Pr | 0.744 | Pr = ν/α |
| **Schmidt number** | Sc | 0.744 | Sc = ν/α (IBAMR notation) |

**Computed:**
- ν = UD/Re = U(1.0)/200 = U/200
- α = ν/Pr = ν/0.744 ≈ 1.344ν

**Note:** Pr = 0.744 is for **air** at ~300K. This validates the solver for air-like fluids.

## Boundary Conditions

### Velocity Field

| Boundary | Condition | Description |
|----------|-----------|-------------|
| **Inlet** | u = (U, 0, 0) | Uniform flow in x-direction |
| **Sphere/Cube surface** | No-slip | u = 0 |
| **Outlet** | ∂u/∂n = 0 | Neumann (zero-gradient) |
| **Sides (y, z)** | ∂u/∂n = 0 | Neumann (zero-gradient) |

### Scalar Field

| Boundary | Condition | Description |
|----------|-----------|-------------|
| **Inlet** | C = C_l | Dirichlet (background = 0) |
| **Sphere/Cube surface** | C = C_h | Dirichlet (isothermal = 1) |
| **Outlet** | ∂C/∂n = 0 | Neumann (zero-gradient) |
| **Sides** | ∂C/∂n = 0 | Neumann (zero-gradient) |

## Grid Resolution

From **Lei et al. (2021)** main simulation:
- 3D nonuniform grid: **352 × 224 × 160**

**IBAMR Recommended:**
- Coarse level: Δx ≈ D/32
- Fine level (AMR): Δx ≈ D/128 near surface
- Finest level: Δx ≈ D/256 in boundary layer

**Computational Cost:**
- 3D is ~100× more expensive than 2D
- Fine grid may require 16-32 cores
- Runtime: 4-12 hours for production run

## Expected Output

### Qualitative Validation (from Figure 3)

✅ **Scalar Contours on Symmetry Plane (y=0):**
- Thermal boundary layer visible
- Wake spreading downstream
- Asymmetric contours due to advection

✅ **Wake Structure:**
- Recirculation bubble (Re=200)
- Scalar transport in wake
- Vortex shedding may occur

### Quantitative Metrics

| Metric | Sphere | Cube | Tolerance |
|--------|--------|------|-----------|
| **Thermal BL thickness** | ~0.1 D | ~0.08 L | ±20% |
| **Wake length** | ~2.0 D | ~1.5 L | ±15% |
| **C = 0.5 isotherm location** | Match ref | Match ref | ±3% |
| **Nusselt number** | Nu ≈ 15-20 | Nu ≈ 18-23 | ±10% |

### Comparison Strategy

1. **Visual comparison** with Richter & Nikrityuk Figure 3
2. **Line plots** along centerline (y=0, z=0)
3. **Scalar isosurfaces** at C = 0.1, 0.5, 0.9
4. **Nusselt number** (dimensionless heat transfer):
   ```
   Nu = (convective heat transfer) / (conductive heat transfer)
   Nu = ∫∫ (∂C/∂n)|_surface dS
   ```

## Pass/Fail Criteria

The test **PASSES** if:

✅ **Stability:**
- [ ] No NaN/Inf values
- [ ] Simulation reaches steady state
- [ ] No solver divergence

✅ **Physical Correctness:**
- [ ] Scalar bounded: 0 ≤ C ≤ 1
- [ ] No negative concentrations
- [ ] Wake structure forms
- [ ] Thermal boundary layer visible

✅ **Quantitative Match (Sphere):**
- [ ] Nusselt number within ±15% of reference
- [ ] C = 0.5 isotherm position within ±5%
- [ ] Wake length within ±20%
- [ ] Boundary layer thickness reasonable

✅ **Quantitative Match (Cube):**
- [ ] Similar metrics as sphere
- [ ] Corner effects visible
- [ ] Edge vortices present

## Building and Running

```bash
cd Test16_3DSphere
mkdir build && cd build
cmake ..
make

# Sphere case
./test16_sphere ../input3d_sphere

# Cube case
./test16_cube ../input3d_cube
```

### With MPI (Recommended for 3D)
```bash
mpirun -np 16 ./test16_sphere ../input3d_sphere
mpirun -np 16 ./test16_cube ../input3d_cube
```

## Expected Runtime

**Sphere (coarse grid 64³):**
- 1 core: ~8 hours
- 4 cores: ~2 hours
- 16 cores: ~30 minutes

**Sphere (fine grid 128³):**
- 16 cores: ~4 hours
- 32 cores: ~2 hours

**Cube:**
- Similar to sphere (slightly more expensive due to corners)

## Output Files

```
Test16_3DSphere/
├── viz_test16_sphere/       ← Sphere visualization
│   ├── dumps.visit
│   ├── lag_data.*.vtu       ← IB mesh
│   └── visit*.*.vtk         ← 3D fields
├── viz_test16_cube/         ← Cube visualization
├── test16_centerline.dat    ← Centerline profiles
├── test16_nusselt.dat       ← Nusselt number vs time
└── test16_results.txt       ← Validation verdict
```

## Visualization with VisIt/ParaView

**Recommended plots:**

1. **Slice plot (y=0 plane):**
   - Scalar contours
   - Velocity vectors
   - Compare with Richter & Nikrityuk Figure 3

2. **Isosurfaces:**
   - C = 0.1, 0.5, 0.9
   - Color by velocity magnitude

3. **Volume rendering:**
   - Scalar field transparency
   - Shows 3D wake structure

4. **Streamlines:**
   - Seeded upstream
   - Color by scalar concentration

## Implementation Notes

### 3D Immersed Boundary

**Sphere vertex file:**
- Triangulated surface mesh
- ~5000-10000 triangles for smooth representation
- Can use MATLAB/Python to generate icosphere

**Cube vertex file:**
- 6 faces × N² points per face
- Sharp edges require dense point spacing
- Corner resolution critical for accuracy

### Nusselt Number Calculation

The Nusselt number quantifies heat/mass transfer:

```cpp
// Integrate normal gradient over IB surface
double total_flux = 0.0;
for (each IB point) {
    Vector3d normal = get_normal(point);
    double dC_dn = compute_gradient_normal(C_field, point, normal);
    double dS = get_surface_area(point);
    total_flux += dC_dn * dS;
}

double Nu = total_flux / (reference_flux);
```

## Physical Interpretation

### Why This Test Matters

1. **3D validation:** Most production cases are 3D
2. **Literature benchmark:** Direct comparison to published CFD
3. **IB method test:** Complex 3D geometries (especially cube)
4. **Wake dynamics:** 3D wakes differ from 2D
5. **Boundary layer resolution:** Tests AMR effectiveness

### Sphere vs Cube Differences

**Sphere:**
- Smooth surface → smooth BL
- Axisymmetric (for aligned flow)
- Well-studied reference case
- Lower drag than cube

**Cube:**
- Sharp edges → flow separation
- Corner vortices
- More challenging for IB method
- Higher drag due to form drag

### Reynolds Number Effects (Re = 200)

At Re = 200:
- Flow is **laminar** but may have weak vortex shedding
- Boundary layer is **well-developed**
- Wake is **2-3 diameters long**
- Recirculation bubble forms behind object

## References

### Primary Reference
**Richter, A., & Nikrityuk, P. A. (2012)**
*"Drag forces and heat transfer coefficients for spherical, cuboidal and ellipsoidal particles in cross flow at sub-critical Reynolds numbers"*
International Journal of Heat and Mass Transfer, 55(4), 1343-1354.
DOI: 10.1016/j.ijheatmasstransfer.2011.09.005

### Validation Context
**Lei, H., Crimaldi, J. P., & Li, C. (2021)**
*"Navigation in odor plumes: How do the flapping kinematics modulate the odor landscape"*
AIAA Aviation 2021 Forum, Paper 2021-2817
DOI: 10.2514/6.2021-2817
**Section II.B, Figure 3** - Uses Richter & Nikrityuk for 3D solver validation

### IBAMR 3D Examples
- `IBAMR/examples/IB/explicit/ex0` - 3D sphere in flow
- `IBAMR/examples/adv_diff/ex0` - 3D scalar transport

## Related Tests

- **Test15_RotatingCylinder**: 2D validation (Yan & Zu 2008)
- **Test08_SphereSource**: 2D sphere (different focus)
- **Test10_MovingIB**: IB-scalar coupling (2D)
- **Test17_PitchPlunge**: Full production case

## Known Limitations

1. **Computational cost:** 3D is expensive, use coarse grid for quick validation
2. **No exact reference data:** Visual comparison with published figures
3. **Cube corners:** Challenging for IB method, may show artifacts
4. **Vortex shedding:** May occur at Re=200 (depends on geometry details)

## Grid Convergence Study (Optional)

To verify convergence, run with multiple grid levels:

| Grid | Δx | Points | Estimated Time (16 cores) |
|------|-----|--------|---------------------------|
| Coarse | D/32 | 64³ ≈ 260k | 30 min |
| Medium | D/64 | 128³ ≈ 2M | 4 hours |
| Fine | D/128 | 256³ ≈ 17M | 30 hours |

**Expected:** Nusselt number converges to within 5% between medium and fine grids.

## Status

- [x] Test specification created
- [ ] Directory structure created
- [ ] main.cpp implemented (sphere case)
- [ ] main.cpp implemented (cube case)
- [ ] input3d_sphere configured
- [ ] input3d_cube configured
- [ ] Sphere vertex file generated
- [ ] Cube vertex file generated
- [ ] Tests executed and validated
- [ ] Results documented

---

**Test ID:** Test16
**Type:** Literature Validation (Richter & Nikrityuk 2012)
**Category:** 3D Solver Validation - Sphere & Cube
**Difficulty:** Advanced (3D + immersed boundaries)
**Dimension:** 3D
**Last Updated:** 2025-11-17
