# Odor Transport Test: Vortex Dynamics in Fish Schools

## Overview

This test demonstrates **how vortex dynamics help undulating bodies spread odor**, based on the research paper:

**"Collective Chemotactic Behavior in Fish Schools"**
*Authors: Maham Kamran, Amirhossein Fardi, Chengyu Li, Muhammad Saif Ullah Khalid*
*arXiv:2408.16136*

## Scientific Background

### The Odor Transport Equation

After computing the velocity field **u** from the IBAMR fluid-structure interaction simulation, the odor concentration field **c(x,y,t)** is governed by the **convection-diffusion equation**:

```
∂c/∂t + u·∇c = D∇²c
```

where:
- **c(x,y,t)** = odor concentration field (scalar)
- **u = (u_x, u_y)** = velocity field from fish swimming (computed by IBAMR)
- **D** = molecular diffusion coefficient
- **u·∇c** = convection term (advection by vortices created by fish)
- **D∇²c** = diffusion term (molecular diffusion)

### Physical Interpretation

1. **Convection term (u·∇c)**: Represents how the vortices created by undulating fish bodies **actively transport** odor molecules through the fluid.

2. **Diffusion term (D∇²c)**: Represents passive **molecular diffusion** of odor (random thermal motion).

3. **Key Question**: How much do vortices enhance odor spreading compared to pure diffusion?

## Test Design

### Comparison Study

The test compares **two scenarios**:

1. **WITH Vortex Dynamics** (Full equation):
   ```
   ∂c/∂t + u·∇c = D∇²c
   ```
   - Uses velocity field from IBAMR simulation
   - Odor is advected by vortices created by fish

2. **WITHOUT Vortices** (Pure diffusion):
   ```
   ∂c/∂t = D∇²c
   ```
   - No velocity field (u = 0)
   - Only molecular diffusion

### Initial Condition

- **Gaussian point source** at position (x₀, y₀):
  ```
  c(x,y,0) = A exp(-r²/(2σ²))
  ```
  where r² = (x-x₀)² + (y-y₀)²

### Numerical Method

The solver uses **finite differences** with:
- **Upwind scheme** for convection (stable for advection-dominated flows)
- **Central differences** for diffusion (2nd-order accurate)
- **Operator splitting**: convection step → diffusion step
- **CFL condition** for stability: CFL = max(|u|Δt/Δx) < 0.5

### Metrics

1. **Spreading width σ(t)**: Standard deviation of concentration distribution
2. **Enhancement factor**: Ratio of spreading widths
   ```
   Enhancement = σ_vortex(t) / σ_diffusion(t)
   ```
3. **Total mass conservation**: ∫∫ c(x,y,t) dx dy = constant

## Usage

### Prerequisites

Install required Python packages:
```bash
pip install numpy matplotlib pyvista scipy
```

For LaTeX rendering (optional):
```bash
# Ubuntu/Debian
sudo apt-get install texlive-latex-extra dvipng

# Or disable LaTeX by setting USE_LATEX = False in the script
```

### Running the Test

```bash
python3 test_odor_transport_vortex_dynamics.py
```

### Expected Input Files

The test reads IBAMR output files:

1. **Eulerian (fluid) data**:
   ```
   ExportEULERIANData/visit_eulerian_db__XXXX/*.vtk
   ```
   - Contains velocity field (U_x, U_y)
   - Contains vorticity field (Omega)

2. **Lagrangian (fish) data**:
   ```
   ExportLagrangianData/visit_lagrangian_db__NN__XXXX.vtk
   ```
   - Contains fish body positions

### Configuration

Edit parameters in the script:

```python
# Odor parameters
D_ODOR = 0.001                 # Diffusion coefficient
SOURCE_X = -2.0                # Initial odor location
SOURCE_Y = 0.0

# Grid resolution
NX = 200                       # Grid points in x
NY = 150                       # Grid points in y

# Frames to process
FRAME_START = 0
FRAME_END = 100
FRAME_SKIP = 10
```

## Output

### Generated Files

All output is saved to `odor_transport_test/` directory:

1. **Comparison frames**: `comparison_XXXX.png`
   - Left panel: Odor WITH vortex dynamics
   - Middle panel: Odor WITHOUT vortices (pure diffusion)
   - Right panel: Enhancement ratio map

2. **Summary plot**: `summary.png`
   - Top: Spreading width σ(t) vs time
   - Bottom: Enhancement factor vs time

### Interpretation

- **Enhancement > 1**: Vortices **enhance** odor spreading
- **Enhancement < 1**: Vortices **suppress** spreading (rare)
- **Enhancement ≈ 1**: Diffusion dominates

Typical results for fish schools: **Enhancement = 2-5x**

## Physical Insights

### Why do vortices enhance odor spreading?

1. **Vortex stretching**: Vortices stretch odor filaments, increasing surface area for diffusion

2. **Chaotic advection**: Complex vortex interactions create chaotic mixing

3. **Increased effective diffusivity**: Convection acts as an "eddy diffusivity"
   ```
   D_effective ≈ D + D_eddy
   where D_eddy ~ u·L (velocity × length scale)
   ```

4. **Persistent vortex structures**: Counter-rotating vortices from fish undulation create coherent transport

### Applications

1. **Collective chemotaxis**: Fish can use vortices to communicate chemically
2. **Predator evasion**: Schools can disperse alarm pheromones rapidly
3. **Foraging**: Enhanced odor spreading helps locate food sources
4. **Environmental sensing**: Improved sampling of chemical gradients

## Technical Details

### Finite Difference Stencils

**Convection (Upwind)**:
```
if u > 0:  ∂c/∂x ≈ (c_i - c_{i-1}) / Δx
if u < 0:  ∂c/∂x ≈ (c_{i+1} - c_i) / Δx
```

**Diffusion (Central)**:
```
∂²c/∂x² ≈ (c_{i+1} - 2c_i + c_{i-1}) / Δx²
```

### Boundary Conditions

- **Neumann (zero flux)**: ∂c/∂n = 0 at all boundaries
- Physical interpretation: No odor escapes the domain

### Stability Conditions

1. **Convection CFL**: Δt < Δx / max(|u|)
2. **Diffusion**: Δt < Δx² / (4D)
3. **Combined**: Δt = min(Δt_conv, Δt_diff)

## Reference

Please cite the original paper if you use this test:

```bibtex
@article{kamran2024collective,
  title={Collective Chemotactic Behavior in Fish Schools},
  author={Kamran, Maham and Fardi, Amirhossein and Li, Chengyu and Khalid, Muhammad Saif Ullah},
  journal={arXiv preprint arXiv:2408.16136},
  year={2024}
}
```

## Related Files

- `test_odor_transport_vortex_dynamics.py` - Main test script
- `plot_combined_fluid_eel.py` - Visualization of vorticity + fish
- `IBEELKinematics.cpp` - IBAMR fish kinematics

## Author

Generated for the Four Fish School IBAMR simulation project.

## License

This test is provided for research and educational purposes.
