# IBAMR Examples Guide

## Available Examples

### 1. Simple IB Cylinder (examples/simple_cylinder)
A basic immersed boundary simulation of a cylinder in uniform flow.

**Features:**
- 2D flow past a circular cylinder
- Standard IB method with Lagrangian markers
- Demonstrates basic IB setup

**To Run:**
```bash
cd examples/simple_cylinder
make
mpirun -np 4 ./main2d input2d
```

### 2. Custom Force Application (examples/custom_force)
Demonstrates how to implement custom force functions for IB structures.

**Features:**
- Custom IBLagrangianForceStrategy implementation
- Time-dependent forcing
- Example of extending IBAMR classes

**To Run:**
```bash
cd examples/custom_force
make
mpirun -np 4 ./main2d input2d
```

## Visualization

All examples output data compatible with VisIt:
```bash
visit -o viz_IB2d/lag_data.visit
visit -o viz_IB2d/dumps.visit
```

## Modifying Examples

1. Edit the `.input` file to change simulation parameters
2. Modify the C++ source to change physics or geometry
3. Rebuild with `make clean && make`
4. Run with desired number of MPI processes

## Common Parameters

Key parameters in `.input` files:
- `N` - Grid resolution
- `DT` - Time step size
- `END_TIME` - Simulation end time
- `REGRID_INTERVAL` - AMR regridding frequency
- `IB_DELTA_FUNCTION` - IB kernel type (IB_4, IB_6, etc.)

## Troubleshooting

**Issue**: Simulation crashes with CFL error
- **Solution**: Reduce `DT` or increase `N`

**Issue**: Poor resolution of IB structure
- **Solution**: Increase Lagrangian point density or use finer grid

**Issue**: Parallel run fails
- **Solution**: Check MPI installation and process count compatibility
