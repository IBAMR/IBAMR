# Test 15: Rotating Cylinder Validation (Yan & Zu 2008)

## Purpose

Validates 2-D advection-diffusion solver against published heat-transfer benchmark from **Yan & Zu (2008)** Lattice Boltzmann Method results.

This test corresponds to **Lei et al. (2021) Section II.B, Figure 2** - "Solver validation" using rotating isothermal cylinder.

## Physical Problem

### Geometry
- 2-D circular cylinder
- Diameter: **D = 1.0** (reference length)
- Rotation: Counter-clockwise
- Domain: Sufficiently large to capture wake

### Governing Equations

**Momentum (Navier-Stokes):**
```
∂u/∂t + u·∇u = -∇p + (1/Re)∇²u
∇·u = 0
```

**Scalar Transport (Temperature ≡ Odor):**
```
∂C/∂t + u·∇C = α∇²C
```

where:
- α = ν/Pr (thermal/scalar diffusivity)
- The temperature equation is **mathematically identical** to odor equation

### Dimensionless Parameters

From **Yan & Zu (2008)** and **Lei et al. (2021)**:

| Parameter | Symbol | Value | Definition |
|-----------|--------|-------|------------|
| **Reynolds number** | Re | 200 | Re = UD/ν |
| **Prandtl number** | Pr | 0.5 | Pr = ν/α |
| **Rotation parameter** | k | 0.5 | k = V/U |

where:
- **U** = free-stream velocity (inlet)
- **V** = tangential velocity of rotating cylinder
- **D** = cylinder diameter = 1.0
- **ν** = kinematic viscosity
- **α** = scalar diffusivity = ν/Pr

**Thus:**
- ν = UD/Re = U(1.0)/200 = U/200
- α = ν/Pr = (U/200)/0.5 = U/100
- V = kU = 0.5U

## Boundary Conditions

### Velocity Field

| Boundary | Condition | Description |
|----------|-----------|-------------|
| **Inlet** | u = (U, 0) | Uniform flow |
| **Cylinder surface** | No-slip + rotation | u_θ = V = 0.5U |
| **Outlet** | ∂u/∂n = 0 | Neumann (zero-gradient) |
| **Lateral (top/bottom)** | ∂u/∂n = 0 | Neumann (zero-gradient) |

### Scalar Field (Temperature/Odor)

| Boundary | Condition | Description |
|----------|-----------|-------------|
| **Inlet** | C = C_l | Dirichlet (background = 0) |
| **Cylinder surface** | C = C_h | Dirichlet (isothermal = 1) |
| **Outlet** | ∂C/∂n = 0 | Neumann (zero-gradient) |
| **Lateral** | ∂C/∂n = 0 | Neumann (zero-gradient) |

## Grid Resolution

From **Lei et al. (2021)** main simulation:
- Base grid: **672 × 416** (2D nonuniform)

**IBAMR Recommended:**
- Coarse level: Δx ≈ D/64
- Fine level (AMR): Δx ≈ D/256
- Refinement around cylinder and wake

## Expected Output

### Qualitative Validation (from Figure 2)

✅ **Streamlines:**
- Recirculation bubble behind cylinder
- Asymmetric wake due to rotation
- Match Yan & Zu (2008) Figure 2

✅ **Scalar (Temperature) Contours:**
- Thermal boundary layer on cylinder
- Wake transport of heated fluid
- Asymmetry consistent with rotation direction
- Steeper gradients on advancing side (rotation + flow)

### Quantitative Metrics

| Metric | Expected | Tolerance |
|--------|----------|-----------|
| **Recirculation length** | ~1.5-2.0 D | ±10% |
| **Thermal boundary layer thickness** | ~0.1 D | ±20% |
| **Stagnation point location** | Depends on k | Visual match |
| **Wake asymmetry** | Consistent with k=0.5 | Qualitative |
| **Scalar field L2 error** | N/A (visual comparison) | ±10% deviation |

### Comparison Strategy

Since exact numerical data from Yan & Zu (2008) may not be available:
1. **Visual comparison** of streamlines and scalar contours with Figure 2
2. **Line plots** along centerline (y=0) for velocity and scalar
3. **Radial profiles** at x = 1D, 2D, 5D downstream
4. **Check physical reasonableness:**
   - No negative concentrations
   - Scalar bounded: C_l ≤ C ≤ C_h
   - Smooth fields, no oscillations
   - Mass/heat flux balance

## Pass/Fail Criteria

The test **PASSES** if:

✅ **Stability:**
- [ ] No NaN/Inf values
- [ ] Simulation reaches steady state (|∂C/∂t| < 1e-6)
- [ ] No solver divergence

✅ **Physical Correctness:**
- [ ] Scalar field bounded: 0 ≤ C ≤ 1
- [ ] No negative concentrations (C ≥ -1e-12)
- [ ] Wake asymmetry consistent with rotation direction
- [ ] Recirculation bubble forms

✅ **Qualitative Match to Yan & Zu (2008):**
- [ ] Streamline pattern matches reference
- [ ] Scalar contour shape matches reference
- [ ] Thermal boundary layer thickness reasonable (±20%)
- [ ] Wake structure similar to published figures

✅ **Quantitative (if reference data available):**
- [ ] Line plot comparisons within ±10% error
- [ ] Drag coefficient within ±15%
- [ ] Nusselt number (heat transfer) within ±15%

## Building and Running

```bash
cd Test15_RotatingCylinder
mkdir build && cd build
cmake ..
make
./test15_rotating_cylinder ../input2d
```

### With MPI
```bash
mpirun -np 4 ./test15_rotating_cylinder ../input2d
```

## Expected Runtime

- Coarse grid (128×128): ~5 min
- Medium grid (256×256): ~20 min
- Fine grid (512×512): ~2 hours

## Output Files

```
Test15_RotatingCylinder/
├── viz_test15/              ← Visualization data
│   ├── dumps.visit
│   ├── lag_data.*.vtu       ← IB cylinder mesh
│   └── visit*.*.vtk         ← Scalar + velocity fields
├── restart_test15/          ← Restart files
├── test15_output.log        ← Simulation log
├── test15_centerline.dat    ← Centerline profiles
└── test15_results.txt       ← Validation verdict
```

## Visualization with VisIt/ParaView

**Recommended plots:**
1. **Streamlines** - colored by velocity magnitude
2. **Scalar contours** - with IB boundary overlay
3. **Vorticity field** - to identify recirculation
4. **Line plots** - centerline y=0 for C(x) and u(x)

## Implementation Notes

### Rotating Cylinder IB Specification

The cylinder is defined using IBAMR's immersed boundary method:

**Lagrangian mesh (.vertex file):**
- Circular boundary: x² + y² = R² where R = D/2 = 0.5
- N = 128 boundary points (sufficient resolution)

**Velocity prescription:**
```
u_θ = V = kU = 0.5U
u_r = 0
```

In Cartesian coordinates at point (x, y) on cylinder:
```
u_x = -u_θ * sin(θ) = -V * y/R
u_y =  u_θ * cos(θ) =  V * x/R
```

### Scalar Boundary Condition on IB

Two approaches:

**Option 1: Source term near IB**
- Add Gaussian source around cylinder
- Smoother for IBAMR scalar transport

**Option 2: Feedback forcing**
- Force C → C_h near IB
- Requires careful tuning of stiffness

This test uses **Option 1** for robustness.

## Physical Interpretation

### Why This Test Matters

1. **Coupled advection-diffusion:** Flow field affects scalar, validating coupling
2. **Curved boundary:** Tests IB method with non-Cartesian geometry
3. **Rotating boundary:** More complex than stationary IB
4. **Published benchmark:** Direct comparison to literature

### Rotation Effects

- **k = 0.5** (moderate rotation): Wake asymmetry, shifted stagnation point
- **Advancing side** (rotation + flow): Thinner boundary layer, higher gradients
- **Receding side** (rotation - flow): Thicker boundary layer, lower gradients

## References

### Primary Reference
**Yan, Y. Y., & Zu, Y. Q. (2008)**
*"Numerical simulation of heat transfer and fluid flow past a rotating isothermal cylinder – A LBM approach"*
International Journal of Heat and Mass Transfer
DOI: 10.1016/j.ijheatmasstransfer.2007.07.053

### Validation Context
**Lei, H., Crimaldi, J. P., & Li, C. (2021)**
*"Navigation in odor plumes: How do the flapping kinematics modulate the odor landscape"*
AIAA Aviation 2021 Forum, Paper 2021-2817
DOI: 10.2514/6.2021-2817
**Section II.B, Figure 2** - Uses Yan & Zu for solver validation

### IBAMR Documentation
- IB Method: https://ibamr.github.io/docs/ibmethod
- Advection-Diffusion: https://ibamr.github.io/docs/advdiff
- Examples: `IBAMR/examples/adv_diff/ex0`

## Related Tests

- **Test08_SphereSource**: Steady source validation (different geometry)
- **Test10_MovingIB**: IB-scalar coupling (different focus)
- **Test14_Benchmarks**: Other literature comparisons
- **Test16_3DSphere**: 3D extension of scalar-IB validation

## Known Limitations

1. **No exact reference data**: Visual comparison only (unless Yan & Zu data obtained)
2. **2D only**: Real flows are 3D (Test16 addresses this)
3. **Steady rotation**: Not time-varying like fish swimming
4. **Temperature ≠ Odor**: Different Prandtl/Schmidt numbers in practice
   - This test: Pr = 0.5 (air-like)
   - Production: Sc = 0.71-340 (water-like) → See Test09

## Status

- [x] Test specification created
- [ ] Directory structure created
- [ ] main.cpp implemented
- [ ] input2d configured
- [ ] CMakeLists.txt created
- [ ] Cylinder vertex file generated
- [ ] Test executed and validated
- [ ] Results documented

---

**Test ID:** Test15
**Type:** Literature Validation (Yan & Zu 2008)
**Category:** Solver Validation - Rotating Cylinder
**Difficulty:** Advanced (rotating IB + scalar)
**Last Updated:** 2025-11-17
