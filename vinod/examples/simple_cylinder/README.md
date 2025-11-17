# Simple Cylinder Example

## Description

This example simulates 2D flow past a circular cylinder using the immersed boundary method.

## Physics

- **Fluid**: Incompressible Navier-Stokes equations
- **Structure**: Rigid circular cylinder (represented by Lagrangian markers)
- **Coupling**: Standard IB method

## Parameters

- Reynolds number: Re = ρ U D / μ ≈ 20 (laminar flow)
- Domain: 4.0 × 4.0
- Cylinder radius: 0.2
- Inflow velocity: 1.0

## Files

- `main.cpp` - Main simulation code
- `input2d` - Input parameters
- `cylinder2d.vertex` - Cylinder geometry points
- `CMakeLists.txt` - Build configuration

## Building

### Using CMake (recommended)
```bash
mkdir build
cd build
cmake ..
make
```

### Using Makefile
```bash
make
```

## Running

```bash
# Serial run
./main2d input2d

# Parallel run with 4 processes
mpirun -np 4 ./main2d input2d
```

## Visualization

Use VisIt to visualize results:
```bash
visit -o viz_IB2d/dumps.visit
```

## Expected Results

- Steady flow pattern around cylinder (low Re)
- Symmetric pressure distribution
- No vortex shedding at this Reynolds number

## Modifications

To increase Reynolds number:
1. Decrease viscosity `MU` in input2d
2. Or increase velocity `U_MAX`
3. May need finer grid resolution for stability

To change cylinder size:
1. Modify `CYLINDER_RADIUS` in input2d
2. Regenerate geometry with modified parameters
