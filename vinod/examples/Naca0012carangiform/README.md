# NACA0012 Undulatory Swimming Simulation

Constraint-based immersed boundary (IB) simulation of undulatory swimming using NACA0012 airfoil profile to model swimmer bodies. Based on IBAMR's `eel2d` example with enhanced swimming kinematics.

## Overview

This simulation models aquatic undulatory locomotion where:
- **NACA0012 chord** represents the swimmer's spine at static equilibrium
- **Traveling waves** propagate along the body to generate thrust
- **Two swimming modes** are supported: anguilliform (eel-like) and carangiform (fish-like)

## Swimming Modes

### 1. Anguilliform (Eel-like Swimming)

**Biological Examples:** Eels, lampreys, sea snakes

**Characteristics:**
- Wave amplitude increases **gradually** from head to tail
- Large amplitude undulations along **entire body**
- Wavelength typically equals or exceeds body length (λ ≥ L)
- High maneuverability, moderate efficiency

**Amplitude Envelope:**
```
A(x/L) = 0.0367 + 0.0323(x/L) + 0.0310(x/L)²
```
- Head amplitude: **A(0) = 0.0367**
- Tail amplitude: **A(1) = 0.0700**
- All coefficients positive → monotonic increase

### 2. Carangiform (Fish-like Swimming)

**Biological Examples:** Tuna, mackerel, jacks, sailfish

**Characteristics:**
- Wave amplitude concentrated in **posterior half**
- Anterior body remains **relatively rigid**
- Wavelength approximately equals body length (λ ≈ L)
- High cruising speed and efficiency

**Amplitude Envelope (Khalid et al. 2016):**
```
A(x/L) = 0.02 - 0.0825(x/L) + 0.1625(x/L)²
```
- Head amplitude: **A(0) = 0.02** (small)
- Tail amplitude: **A(1) = 0.10** (5× larger)
- Negative linear term creates U-shaped envelope

## Mathematical Formulation

**Centerline Displacement:**
```
y(x,t) = A(x/L) × cos[2π(x/λ - ft)]
```

**Deformation Velocity:**
```
∂y/∂t = A(x/L) × 2πf × sin[2π(x/λ - ft)]
```

**Parameters:**
- `x` = chordwise position along spine [0, L]
- `L` = body length (chord length) = 1
- `λ` = wavelength of undulation = L
- `f` = swimming frequency = 1 Hz
- `t` = time (s)

## Directory Structure

```
Naca0012carangiform/
├── naca0012_swimmer_generator.m      # Unified mesh generator (RECOMMENDED)
├── naca_anguilliform.m                # Anguilliform-specific generator
├── naca_carangiform_clean.m           # Carangiform-specific generator
├── input2d                            # IBAMR input configuration
├── IBNACA0012Kinematics.cpp           # C++ kinematics implementation
├── IBNACA0012Kinematics.h             # Header file
├── example.cpp                        # Main simulation driver
├── CMakeLists.txt                     # Build configuration
├── README.md                          # This file
├── INSTALL.md                         # Installation instructions
├── CONTRIBUTING.md                    # Contribution guidelines
├── CODE_OF_CONDUCT.md                 # Community code of conduct
├── LICENSE                            # BSD 3-Clause license
├── CITATION.cff                       # Citation information
└── .github/                           # GitHub templates and workflows
    ├── workflows/c-cpp.yml            # CI/CD workflow
    ├── ISSUE_TEMPLATE/                # Issue templates
    └── PULL_REQUEST_TEMPLATE.md       # PR template
```

## Documentation

- **[INSTALL.md](INSTALL.md)** - Detailed installation instructions for IBAMR and dependencies
- **[CONTRIBUTING.md](CONTRIBUTING.md)** - Guidelines for contributing to the project
- **[MESH_PARAMETER_GUIDE.md](MESH_PARAMETER_GUIDE.md)** - Guide to mesh generation parameters
- **[FORMATTING_GUIDE.md](FORMATTING_GUIDE.md)** - Code formatting guidelines
- **[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)** - Community guidelines

## Quick Start

### Prerequisites

Before starting, ensure you have:
- IBAMR installed (see [INSTALL.md](INSTALL.md) for detailed instructions)
- MATLAB or GNU Octave for mesh generation
- MPI for parallel execution

### Method 1: Unified Generator (Recommended)

1. **Open** `naca0012_swimmer_generator.m`
2. **Set swimming mode** (line 20):
   ```matlab
   swimming_mode = 'carangiform';  % or 'anguilliform'
   ```
3. **Run** the script in MATLAB
4. **Outputs:**
   - `naca0012[mode].vertex` - Mesh file for IBAMR
   - Visualization plots (amplitude envelope, backbone motion, mesh)

### Method 2: Mode-Specific Generators

**For carangiform:**
```matlab
run naca_carangiform_clean.m
```

**For anguilliform:**
```matlab
run naca_anguilliform.m
```

## IBAMR Simulation Setup

### Step 1: Generate Mesh
```bash
cd /path/to/Naca0012carangiform
matlab -batch "swimming_mode='carangiform'; run('naca0012_swimmer_generator.m')"
```

### Step 2: Configure input2d

Edit `input2d` file:

**A. Update structure name (line 149, 237):**
```
structure_names = "naca0012carangiform"  # or "naca0012anguilliform"
```

**B. Set kinematics equations (lines 198-202 or 215-218):**

For **CARANGIFORM** (uncomment lines 198-202):
```
body_shape_equation = "(0.02 - 0.0825*X_0 + 0.1625*X_0^2) * cos(2*PI*X_0 - 2*PI*T)"
deformation_velocity_function_0 = "((0.02 - 0.0825*X_0 + 0.1625*X_0^2) * 2*PI * sin(2*PI*X_0 - 2*PI*T)) * N_0"
deformation_velocity_function_1 = "((0.02 - 0.0825*X_0 + 0.1625*X_0^2) * 2*PI * sin(2*PI*X_0 - 2*PI*T)) * N_1"
```

For **ANGUILLIFORM** (uncomment lines 215-218):
```
body_shape_equation = "(0.0367 + 0.0323*X_0 + 0.0310*X_0^2) * cos(2*PI*X_0 - 2*PI*T)"
deformation_velocity_function_0 = "((0.0367 + 0.0323*X_0 + 0.0310*X_0^2) * 2*PI * sin(2*PI*X_0 - 2*PI*T)) * N_0"
deformation_velocity_function_1 = "((0.0367 + 0.0323*X_0 + 0.0310*X_0^2) * 2*PI * sin(2*PI*X_0 - 2*PI*T)) * N_1"
```

### Step 3: Build and Run

```bash
# Configure and build
mkdir build && cd build
cmake ..
make

# Run simulation
mpirun -np 4 ./example input2d

# Visualize results
visit -o viz_naca0012_Str/dumps.visit
```

## Simulation Parameters

### Domain Configuration
- **Domain size:** 8L × 4L (where L = 1)
- **Grid resolution:** 2048 × 1024 (coarse level: 64 × 32)
- **AMR levels:** 3
- **Refinement ratio:** 4
- **Finest spacing:** dx = dy ≈ 0.00391

### Fluid Properties
- **Reynolds number:** Re = 5609
- **Viscosity:** μ = 0.785/Re
- **Density:** ρ = 1.0
- **Swimming frequency:** f = 1 Hz

### Boundary Conditions
- **All boundaries:** No-slip walls (Dirichlet)
- **Periodicity:** Both x and y directions (matches eel2d reference)

## Output Files

### Lagrangian Structure Outputs
Location: `./NACA0012Str/`
- `NACA0012_Drag.curve` - Drag force vs time
- `NACA0012_Torque.curve` - Torque vs time
- `NACA0012_CenterOfMass.curve` - Center of mass position
- `NACA0012_RigidTransVelocity.curve` - Rigid-body velocity

### Eulerian Field Outputs
Location: `./viz_naca0012_Str/`
- Velocity fields (u, v)
- Pressure field (p)
- Vorticity field (ω)
- Divergence of velocity (∇·u)

Visualization interval: Every 40 time steps

## Key Equations in Code

### Variables in Parser
- `X_0` = x/L (normalized chordwise position, 0 to 1)
- `T` = t (time)
- `N_0`, `N_1` = normal vector components (auto-computed)
- `PI` = π

### NACA0012 Thickness Distribution
```
y_t = 5t_max[0.2969√(x/L) - 0.1260(x/L) - 0.3516(x/L)² + 0.2843(x/L)³ - 0.1015(x/L)⁴]
```
where `t_max = 0.12 × 0.2489 ≈ 0.0299` (scaled for point density)

## References

1. **Khalid, M.S.U., et al.** (2020). "Flow transitions and mapping for undulating swimmers." *Physical Review Fluids*, 5(6):063104. [DOI:10.1103/PhysRevFluids.5.063104](https://doi.org/10.1103/PhysRevFluids.5.063104)

2. **Lighthill, M.J.** (1960). "Note on the swimming of slender fish." *Journal of Fluid Mechanics*, 9(2):305-317.

3. **Sfakiotakis, M., et al.** (1999). "Review of fish swimming modes for aquatic locomotion." *IEEE Journal of Oceanic Engineering*, 24(2):237-252.

4. **IBAMR:** https://github.com/IBAMR/IBAMR
   - Based on `examples/ConstraintIB/eel2d/`

## Switching Between Modes

### Quick Switch Checklist

| Step | Carangiform | Anguilliform |
|------|-------------|--------------|
| **1. Generate mesh** | `swimming_mode = 'carangiform'` | `swimming_mode = 'anguilliform'` |
| **2. Vertex file** | `naca0012carangiform.vertex` | `naca0012anguilliform.vertex` |
| **3. Structure name** | `"naca0012carangiform"` | `"naca0012anguilliform"` |
| **4. Shape equation** | `0.02 - 0.0825*X_0 + 0.1625*X_0^2` | `0.0367 + 0.0323*X_0 + 0.0310*X_0^2` |
| **5. Rebuild** | `make` | `make` |

## Troubleshooting

### Issue: Mesh file not found
**Solution:** Ensure vertex file name matches `structure_names` in `input2d`

### Issue: Simulation diverges
**Solution:**
- Check CFL condition (CFL_MAX = 0.3)
- Verify amplitude envelope doesn't create extreme deformations
- Reduce DT_MAX if needed

### Issue: Wrong swimming behavior
**Solution:** Verify equations in `input2d` match the mesh mode (anguilliform vs carangiform)

## Performance Notes

- **Typical runtime:** 2-4 hours for 10 time units on 4 cores
- **Memory usage:** ~8-16 GB depending on AMR refinement
- **Recommended processors:** 4-16 MPI ranks
- **Output size:** ~5-10 GB for full simulation

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Reporting bugs
- Suggesting enhancements
- Submitting pull requests
- Coding standards

Please read our [Code of Conduct](CODE_OF_CONDUCT.md) before contributing.

## License

This code is distributed under the 3-clause BSD license. See [LICENSE](LICENSE) file for details.

This project uses IBAMR, which is separately licensed under the BSD 3-Clause License.

## Citation

If you use this code in your research, please cite it using the information in [CITATION.cff](CITATION.cff), or cite:

**This software:**
```bibtex
@software{naca0012_swimming,
  author = {Thale, Vinod},
  title = {NACA0012 Undulatory Swimming Simulation},
  url = {https://github.com/vinodthale/Naca0012carangiform},
  license = {BSD-3-Clause}
}
```

**Key references:**
- IBAMR: https://github.com/IBAMR/IBAMR
- Khalid, M.S.U., et al. (2020). "Flow transitions and mapping for undulating swimmers." *Physical Review Fluids*, 5(6):063104. [DOI](https://doi.org/10.1103/PhysRevFluids.5.063104)

## Support

For questions or issues:
- Check the [documentation](#documentation) first
- Review [existing issues](https://github.com/vinodthale/Naca0012carangiform/issues)
- Open a new issue using the appropriate template
- See [INSTALL.md](INSTALL.md) for troubleshooting tips
