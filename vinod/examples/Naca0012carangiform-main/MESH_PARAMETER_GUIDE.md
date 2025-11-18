# NACA0012 Carangiform Mesh Generation Guide

## Overview
This guide explains the relationships between domain size, grid spacing, and Lagrangian mesh generation for the NACA0012 carangiform swimming simulation.

---

## 1. Background Eulerian Mesh Parameters

### Domain Size and Grid Resolution

The computational domain is defined by:

```matlab
% Background Mesh Properties
Ny = 16*4*4*4 = 1024;     % Grid points in y-direction
Ly = 4;                    % Domain height (4L, where L=1)
Nx = 32*4*4*4 = 2048;     % Grid points in x-direction
Lx = 8;                    % Domain length (8L)
```

### Grid Spacing Calculation

The uniform grid spacing is calculated as:

```matlab
dy = Ly/Ny = 4/1024 = 0.00390625 ≈ 3.906e-3
dx = Lx/Nx = 8/2048 = 0.00390625 ≈ 3.906e-3
```

**Key Insight:** The grid is **isotropic** (dx = dy), which is important for:
- Uniform IB kernel spreading/interpolation
- Consistent resolution in all directions
- Accurate force distribution

### Domain Layout

```
Domain: 8L × 4L
         ┌─────────────────────────────────────┐
         │                                     │
    -1.52│              FLOW →                 │ 2.48
         │        ──────────                   │
         │       /NACA0012  \                  │
         │      └────────────┘                 │
         │       centered at origin            │
         └─────────────────────────────────────┘
      -6.52                                  1.48

Fish position: x ∈ [-0.5, 0.5], y ≈ [-0.06, 0.06]
Clearance: ~6.5L upstream, ~1.5L downstream, ±1.5L lateral
```

---

## 2. NACA0012 Geometry Parameters

### Fish Body Length
```matlab
L = 1;  % Non-dimensional fish body length
```

### Thickness Parameters
```matlab
thickness_scale = 0.2489;                    % Scale factor
t_max = 0.12 * thickness_scale;              % Effective max thickness
t_max = 0.02987 ≈ 0.03                       % ~3% thickness (scaled)
```

**Why scale down from 12% to 3%?**
- Standard NACA0012 is 12% thick (t/c = 0.12)
- Real fish are thinner: typically 4-8% thickness
- Scaling by 0.2489 gives effective thickness ≈ 3%
- This matches eel2d point density and produces realistic fish proportions

### NACA0012 Thickness Distribution
For a point at normalized position x/L ∈ [0, 1]:

```matlab
y_t(x/L) = 5*t_max*(0.2969*√(x/L) - 0.1260*(x/L)
                    - 0.3516*(x/L)² + 0.2843*(x/L)³
                    - 0.1015*(x/L)⁴)
```

This gives the **half-thickness** at each location.
**Full thickness** = 2 × y_t × L

---

## 3. Lagrangian Mesh Generation

### Number of Cross-Sections
```matlab
BodyNx = ceil(L/dx) = ceil(1/0.00390625) = 256 cross-sections
```

**Meaning:** The fish body is discretized into 256 vertical slices along its length.

### Cross-Section Point Distribution

At each cross-section i (where i = 1 to 256):

1. **Streamwise position:**
   ```matlab
   x = (i-1)*dx
   x_norm = x/L        % Normalized position: 0 → 1
   ```

2. **Centerline displacement (Carangiform):**
   ```matlab
   A(x/L) = 0.02 - 0.0825*(x/L) + 0.1625*(x/L)²
   y_centerline = A(x/L) * cos(2π*x/L)  % at t=0
   ```

3. **Thickness at this location:**
   ```matlab
   y_t = 5*t_max*(0.2969*√(x/L) - 0.1260*(x/L) - ...)
   height = 2*y_t*L
   ```

4. **Number of points across thickness:**
   ```matlab
   NumPointsInHeight = max(1, ceil(height/dy))
   ```

   This ensures the Lagrangian mesh has approximately **dy spacing** normal to the centerline.

5. **Generate points:**
   - Points are distributed vertically (normal to flow)
   - Spacing: approximately dy between points
   - Creates "up" and "down" sides of the cross-section

### Total Lagrangian Points

For your vertex file:
```
Total IB points = 2938
Cross-sections  = 256
Points/section  ≈ 11.5 (average)
```

**Why varying points per section?**
- Near head/tail: thinner → fewer points
- At maximum thickness: thicker → more points
- Ensures consistent resolution: ~dy spacing everywhere

---

## 4. Key Relationships

### 4.1 Grid Spacing → Lagrangian Resolution

```
dx = dy = Ly/Ny = Lx/Nx
       ↓
BodyNx = ceil(L/dx)     [streamwise discretization]
       ↓
NumPointsInHeight = ceil(height/dy)  [normal discretization]
```

### 4.2 Resolution Requirements

For accurate IB simulations:
- **Lagrangian spacing ≈ Eulerian spacing:** Ensures proper IB kernel support
- **IB_4 kernel** requires ~4 grid points per IB wavelength
- Your setup: dx = dy = ds_lag ≈ 0.0039 ✓ Well-resolved

### 4.3 Domain Size → Fish Placement

```
Domain: Lx × Ly = 8L × 4L
Fish:   1L × 0.06L (approximately)
Ratio:  Fish occupies ~1/8 × 1/67 of domain
```

This provides sufficient room for:
- Wake development downstream
- Undisturbed flow upstream
- Lateral clearance for transverse motion

---

## 5. AMR (Adaptive Mesh Refinement) Details

From `input2d`:
```c
N = 64                    // Coarsest grid: 64 cells in x
MAX_LEVELS = 3            // 4 total levels (0,1,2,3)
REF_RATIO = 4             // Each level 4× finer

Level 0: 64 × 32 cells    (coarsest)
Level 1: 256 × 128 cells  (4× refinement)
Level 2: 1024 × 512 cells (4× refinement)
Level 3: 4096 × 2048 cells (finest)  ← Matches your Nx × Ny if fully refined
```

**Effective resolution at finest level:**
```matlab
dx_finest = Lx/(2*N*4^3) = 8/(2*64*64) = 8/8192 = 0.000977
```

Wait, this doesn't match your dx = 0.00390625!

### Resolution Clarification

Your MATLAB parameters suggest:
```matlab
Nx = 2048, Lx = 8  →  dx = 8/2048 = 0.00390625
Ny = 1024, Ly = 4  →  dy = 4/1024 = 0.00390625
```

The `input2d` file has:
```c
N = 64, MAX_LEVELS = 3, REF_RATIO = 4
→ Finest grid: 64*4^3 × 32*4^3 = 4096 × 2048
```

**Domain from input2d:**
```c
x_lo = -6.52, -1.52
x_up =  1.48,  2.48
→ Lx = 1.48-(-6.52) = 8.0
→ Ly = 2.48-(-1.52) = 4.0
```

So finest spacing should be:
```
dx = 8.0/4096 = 0.001953125
dy = 4.0/2048 = 0.001953125
```

This is **2× finer** than your MATLAB parameters!

---

## 6. Recommended Resolution Matching

### Option A: Match input2d (Finer Resolution)
```matlab
% For MAX_LEVELS=3, REF_RATIO=4, N=64
Nx = 2*64*4^3 = 8192;    % 8192 grid points in x
Ny = 64*4^3 = 4096;      % 4096 grid points in y
dx = Lx/Nx = 0.000977;   % Finer spacing
dy = Ly/Ny = 0.000977;
```

This will generate ~5876 IB points (2× more resolution).

### Option B: Match MATLAB (Coarser, Current)
Adjust input2d:
```c
N = 128                  // Double coarse grid
MAX_LEVELS = 2           // Reduce to 3 levels total
```

This gives: 128×4^2 = 2048 cells in x at finest level ✓

---

## 7. Generating the .vertex File

### Run the Generator
```bash
matlab -batch "naca0012_swimmer_generator"
```

### Output Files
1. **naca0012carangiform.vertex** - Lagrangian mesh coordinates
2. **amplitude_envelope_carangiform.png** - A(x/L) visualization
3. **backbone_motion_carangiform.png** - Swimming kinematics
4. **mesh_visualization_carangiform.png** - IB point distribution

### Vertex File Format
```
Line 1: NumMatPoints (e.g., 2938)
Line 2+: x_coord    y_coord
         -0.500000  0.020000
         -0.500000  0.016094
         ...
```

---

## 8. Parameter Sensitivity

### Effect of dx/dy on Mesh
| dx = dy | BodyNx | Est. IB Points | Resolution |
|---------|---------|----------------|------------|
| 0.00781 | 128     | ~730           | Coarse     |
| 0.00391 | 256     | ~2938          | Medium     |
| 0.00195 | 512     | ~11750         | Fine       |
| 0.00098 | 1024    | ~47000         | Very Fine  |

**Tradeoff:** Finer resolution → more accurate IB → higher computational cost

### Effect of thickness_scale
| Scale  | t_max  | Effective t/c | Appearance      |
|--------|--------|---------------|-----------------|
| 1.0    | 0.120  | 12%           | Thick (airfoil) |
| 0.5    | 0.060  | 6%            | Medium (fish)   |
| 0.2489 | 0.030  | 3%            | Thin (eel-like) |

---

## 9. Verification Checklist

Before running simulation:

- [ ] **Grid spacing:** dx = dy (isotropic) ✓
- [ ] **Lagrangian resolution:** ds ≈ dx ✓
- [ ] **Domain size:** Fish << Domain ✓
- [ ] **Clearance:** ≥1L in all directions ✓
- [ ] **AMR levels:** Finest level matches Nx,Ny
- [ ] **Vertex file:** Points = 2938 (carangiform)
- [ ] **Structure name:** "naca0012carangiform" in input2d
- [ ] **Kinematics:** Amplitude coefficients match mesh

---

## 10. Quick Reference

### Current Configuration
```
Domain:          8L × 4L
Grid (MATLAB):   2048 × 1024
Grid (input2d):  4096 × 2048 (finest AMR)  ⚠️ MISMATCH
Spacing:         dx = dy = 0.00391 (MATLAB)
                 dx = dy = 0.00195 (input2d)
Fish length:     L = 1
IB points:       2938
Cross-sections:  256
Mode:            Carangiform
```

### Swimming Kinematics (Carangiform)
```
Amplitude:   A(x/L) = 0.02 - 0.0825(x/L) + 0.1625(x/L)²
Centerline:  y(x,t) = A(x/L)*cos(2π(x/λ - ft))
Frequency:   f = 1 Hz
Wavelength:  λ = L = 1
```

---

## 11. Next Steps

1. **Resolve resolution mismatch:**
   - Either adjust MATLAB Nx,Ny to match input2d
   - Or adjust input2d N,MAX_LEVELS to match MATLAB

2. **Regenerate vertex file** with consistent parameters

3. **Verify in IBAMR:**
   ```bash
   ./main2d input2d
   ```

4. **Check outputs:**
   - Mesh deformation looks smooth
   - No IB points leaving domain
   - Forces/velocities physically reasonable

---

## References

- IBAMR eel2d example: `examples/ConstraintIB/eel2d/`
- NACA airfoil theory: Abbott & Von Doenhoff (1959)
- Carangiform kinematics: Khalid et al. (2016) *Phys. Rev. Fluids* 5:063104

---

**Generated:** 2025-11-13
**Configuration:** NACA0012 Carangiform, Lx=8L, Ly=4L, dx=dy=0.00391
