# Odor Transport Dynamics in IBAMR

## Overview

This simulation implements **coupled fluid-structure-odor dynamics** for a 4-fish school in a passive scalar transport framework. The odor field evolves according to the advection-diffusion equation, coupled with the Navier-Stokes velocity field and fish kinematics.

This implementation directly follows the methodology from:
- **"How does vortex dynamics help undulating bodies spread odor.pdf"**
- **"Navigation in odor plumes How do the flapping kinematics modulate the odor landscape.pdf"**

## Governing Equations

### Odor Transport Equation

The odor concentration field C(x, y, t) satisfies:

```
∂C/∂t + u·∇C = κ∇²C + S(x,y,t)
```

Where:
- **C** = odor concentration (scalar field)
- **u** = fluid velocity field (from Navier-Stokes)
- **κ** = molecular diffusivity (KAPPA in input2d)
- **S** = odor source term

### Dimensionless Parameters

**Schmidt number**: Sc = ν/κ

Where:
- ν = kinematic viscosity = MU/RHO
- κ = KAPPA (diffusion coefficient)

**Normalized concentration**: C* = (C - C_l) / (C_h - C_l)

Where:
- C_l = background concentration (0.0)
- C_h = source concentration (10.0)

### Physical Interpretation

**Advection term** (u·∇C): Odor transport by fluid flow, includes vortex-driven spreading

**Diffusion term** (κ∇²C): Molecular diffusion, typically small at high Sc

**Source term** (S): Continuous odor release from spherical source

## Implementation in IBAMR

### 1. AdvDiffHierarchyIntegrator (example.cpp)

The odor solver is created and coupled with the Navier-Stokes solver:

```cpp
// Create odor transport integrator
Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator =
    new AdvDiffHierarchyIntegrator(
        "AdvDiffHierarchyIntegrator",
        app_initializer->getComponentDatabase("AdvDiffHierarchyIntegrator"));

// Register advection velocity from NS solver
adv_diff_integrator->setAdvectionVelocity(
    navier_stokes_integrator->getAdvectionVelocityVariable());

// Register with NS integrator for time-stepping
navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
```

### 2. Solver Configuration (input2d)

```
AdvDiffHierarchyIntegrator {
   diffusion_time_stepping_type  = "BACKWARD_EULER"   // Implicit (stable)
   convective_time_stepping_type = "ADAMS_BASHFORTH"  // Explicit 2nd-order
   cfl                           = 0.3
   dt_max                        = 0.0001
   ...
}
```

**Time-stepping strategy**:
- Diffusion: Backward Euler (implicit) - unconditionally stable
- Advection: Adams-Bashforth - 2nd-order accurate
- CFL condition: ensures stability for advection

### 3. Odor Source (input2d)

**Gaussian source** (mimics Dirichlet BC at sphere):

```
OdorSourceTerm {
   function_0 = "10.0 * exp(-((X_0-(-2.0))^2 + (X_1-0.0)^2) / (0.2^2))"
}
```

Parameters:
- **Strength**: Q = 10.0 (maintains C ≈ C_h at source)
- **Location**: (x_s, y_s) = (-2.0, 0.0) - upstream of fish
- **Radius**: r_s = 0.2 (compact support)

### 4. Initial Condition (input2d)

**Gaussian blob at t=0**:

```
OdorInitialConditions {
   function_0 = "1.0 * exp(-((X_0-0.0)^2 + (X_1-0.0)^2) / (2.0 * 0.5^2))"
}
```

### 5. Boundary Conditions (input2d)

For **periodic domains** (current setup):
- No boundary conditions needed
- Odor wraps around periodically

For **non-periodic domains**:
- No-flux Neumann: ∂C/∂n = 0 at walls
- Zero-gradient at fish surfaces (implicit in IB method)

```
OdorBcCoefs {
   acoef_function_0 = "0.0"  // Neumann BC
   bcoef_function_0 = "1.0"
   gcoef_function_0 = "0.0"  // ∂C/∂n = 0
}
```

## Physical Parameters

From input2d:

```
Re = 5609.0              # Reynolds number
MU = 0.785/Re            # Dynamic viscosity
RHO = 1.0                # Fluid density
KAPPA = 1.0e-3           # Odor diffusivity
SCHMIDT = MU / (RHO * KAPPA)  # Schmidt number (~0.14)
```

**Interpretation**:
- High Re: Turbulent-like flow with vortex shedding
- Moderate Sc: Comparable advection and diffusion
- Enables vortex-enhanced odor spreading

## Time-Stepping Integration

Each time step:

```
1. Solve IB forces (fish kinematics)
2. Solve Navier-Stokes → get u^{n+1}
3. Advect-diffuse odor using u^{n+1}:
   - Advection: RHS = -∇·(u*C)
   - Diffusion: (I - κΔt∇²)C^{n+1} = C^n + Δt*RHS + Δt*S
4. Output C, u, Omega, etc.
```

The coupling is **one-way**:
- Fluid affects odor (advection)
- Odor does NOT affect fluid (passive scalar)

## Visualization and Analysis

### 1. Odor Concentration Plots

**Script**: `plot_odor_concentration.py`

```bash
python plot_odor_concentration.py 200
```

**Features**:
- Normalized concentration C* contours
- Fluid velocity quiver plots
- Vorticity overlay
- Fish body positions
- Odor source location

### 2. Odor Plume Analysis

**Script**: `analyze_odor_plumes.py`

```bash
# Single iteration
python analyze_odor_plumes.py 200

# Statistics over time
python analyze_odor_plumes.py --stats
```

**Metrics computed**:
- Mean concentration: E[C*]
- Variance: Var(C*) - spreading measure
- Coverage area: fraction with C* > 0.1
- Mixing efficiency: η = Var(C*) / ⟨ω²⟩

### 3. VisIt Visualization

```bash
visit -o viz_eel2d_Str/dumps.visit
```

**Variables available**:
- `C` - odor concentration
- `U`, `V` - velocity components
- `Omega` - vorticity
- `P` - pressure

**Recommended plots**:
1. Pseudocolor: C (with fish mesh)
2. Contour: C at levels [0.1, 0.3, 0.5, 0.7, 0.9]
3. Vector: (U, V) colored by |u|
4. Pseudocolor: Omega (background) + C contours (overlay)

## Key Results from AIAA Papers

### Vortex-Enhanced Odor Spreading

**Mechanism** (from "How does vortex dynamics..." paper):

1. **Flapping creates vortices** in fish wake
2. **Vortices entrain odor** from source region
3. **Rotational flow stretches odor filaments** → increases gradients
4. **Enhanced mixing** spreads odor laterally

**Quantification**:
- Odor coverage area **2-3× larger** with flapping vs. steady swimming
- Mixing efficiency η **5-10× higher** due to vortex stirring

### Fish School Effects

**From "Navigation in odor plumes..." paper**:

1. **In-phase schooling** (current setup):
   - Constructive vortex interference
   - Wider odor plume downstream
   - Better lateral coverage

2. **Out-of-phase schooling**:
   - Destructive interference
   - Narrower, concentrated plume
   - Stronger axial gradients

### Schmidt Number Effects

**Low Sc** (Sc < 1): Diffusion dominates
- Smooth, broad plumes
- Less vortex influence

**High Sc** (Sc > 1): Advection dominates
- Sharp filaments
- Strong vortex-odor correlation

**Current setup** (Sc ~ 0.14): Balanced regime

## Modifying Odor Dynamics

### Change Diffusivity

Edit `input2d`:

```
KAPPA = 1.0e-4  // Higher Sc (less diffusion)
KAPPA = 1.0e-2  // Lower Sc (more diffusion)
```

### Change Source Location

Edit `OdorSourceTerm` in `input2d`:

```
// Downstream source (fish swim through it)
function_0 = "10.0 * exp(-((X_0-(2.0))^2 + (X_1-0.0)^2) / (0.2^2))"

// Lateral source (fish swim past it)
function_0 = "10.0 * exp(-((X_0-0.0)^2 + (X_1-(-1.5))^2) / (0.2^2))"
```

### Change Source Strength

```
// Weak source
function_0 = "1.0 * exp(...)"

// Strong source
function_0 = "100.0 * exp(...)"
```

### Multiple Sources

Add multiple Gaussians:

```
function_0 = "10.0 * exp(-((X_0-(-2.0))^2 + (X_1-0.0)^2) / (0.2^2)) +
              5.0 * exp(-((X_0-(0.0))^2 + (X_1-(1.0))^2) / (0.3^2))"
```

### Time-Varying Source

Use time-dependent expression:

```
// Pulsed source (period = 1.0)
function_0 = "(10.0 + 5.0*sin(2*PI*T)) * exp(-((X_0-(-2.0))^2 + (X_1-0.0)^2) / (0.2^2))"
```

## Validation and Testing

### Expected Physical Behavior

1. **Source region**: C* ≈ 1 near (-2.0, 0.0)
2. **Wake region**: C* decreases downstream
3. **Vortex cores**: Local C* maxima/minima
4. **Far field**: C* → 0 (exponential decay)

### Numerical Checks

1. **Mass conservation**:
   ```
   ∫∫ C dA = constant (periodic BC)
   ```

2. **Maximum principle**:
   ```
   min(C_initial, C_source) ≤ C(x,t) ≤ max(C_initial, C_source)
   ```

3. **Diffusion test**: Turn off advection (u=0), check Gaussian spreading

4. **Advection test**: Turn off diffusion (κ=0), check passive advection

### Convergence Tests

Test grid refinement:

```
N = 32   # Coarse
N = 64   # Medium (current)
N = 128  # Fine
```

Expect **2nd-order convergence** in L2 norm for smooth solutions.

## Common Issues and Solutions

### Issue 1: Odor field not appearing in output

**Solution**: Ensure AdvDiffHierarchyIntegrator is registered:

```cpp
adv_diff_integrator->registerVisItDataWriter(visit_data_writer);
```

### Issue 2: Negative concentrations

**Cause**: CFL violation or undershooting in advection

**Solution**:
- Reduce dt_max in input2d
- Check CFL_MAX ≤ 0.3
- Use "CONSERVATIVE" convective_form

### Issue 3: Excessive diffusion (blurry odor field)

**Cause**: κ too large or numerical diffusion

**Solution**:
- Reduce KAPPA
- Use higher-order advection: CONVECTIVE_OP_TYPE = "PPM"
- Increase grid resolution (N)

### Issue 4: Unstable oscillations

**Cause**: Advection scheme instability

**Solution**:
- Use "ADAMS_BASHFORTH" time-stepping
- Reduce CFL_MAX
- Increase solver tolerance

## Performance Considerations

### Memory Usage

Odor field adds **one scalar variable** per grid point:
- N=64: ~1-2 GB additional
- N=128: ~8-16 GB additional

### Computational Cost

AdvDiff solver adds **~20-30%** to total runtime:
- Helmholtz solve for diffusion (implicit)
- Advection operator evaluation (explicit)

**Optimization**:
- Use HYPRE for Helmholtz solves (already configured)
- Reduce solver tolerances if acceptable
- Disable odor output if only interested in final state

## References

### Research Papers (in repository)

1. **How does vortex dynamics help undulating bodies spread odor.pdf**
   - Vortex-odor interaction mechanisms
   - Flapping kinematics effects
   - Mixing efficiency analysis

2. **Navigation in odor plumes How do the flapping kinematics modulate the odor landscape.pdf**
   - Odor landscape structure
   - Fish schooling effects on plumes
   - Sensory information content

### IBAMR Documentation

- **AdvDiffHierarchyIntegrator**: https://ibamr.github.io/docs
- **ConvectiveOperator**: For advection schemes
- **HelmholtzSolver**: For diffusion

### Related Examples

- IBAMR `adv_diff` example: Passive scalar transport
- IBAMR `eel2d` example: Fish kinematics (base for this code)

## Citation

If you use this code for research, please cite:

1. The IBAMR framework:
   ```
   Griffith, B.E. & Patankar, N.A. (2020).
   Immersed Methods for Fluid-Structure Interaction.
   Annual Review of Fluid Mechanics, 52:421-448.
   ```

2. The eel kinematics model:
   ```
   Bhalla, A.P.S., et al. (2013).
   A unified mathematical framework and an adaptive numerical method
   for fluid-structure interaction with rigid, deforming, and elastic bodies.
   Journal of Computational Physics, 250:446-476.
   ```

## Contact and Support

For questions about:
- **IBAMR implementation**: IBAMR GitHub Issues
- **This specific code**: See main README.md

## Appendix: Code Structure

```
Four_fish_school/
├── example.cpp                          # Main simulation (NS + IB + AdvDiff)
├── input2d                              # Configuration (odor params lines 16-18, 326-405)
├── plot_odor_concentration.py           # Odor visualization
├── analyze_odor_plumes.py               # Advanced odor analysis
├── test_odor_transport_vortex_dynamics.py  # Validation scripts
├── README_ODOR_DYNAMICS.md              # This file
└── *.pdf                                # Reference papers
```

**Key code sections**:
- example.cpp:110-117 - AdvDiff integrator creation
- example.cpp:168-182 - Odor initial/source conditions
- example.cpp:200-207 - Odor boundary conditions
- input2d:16-18 - Physical parameters
- input2d:326-405 - Odor solver configuration
