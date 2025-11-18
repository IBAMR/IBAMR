# Installation Guide

This guide provides detailed instructions for installing dependencies and building the NACA0012 Carangiform swimming simulation.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installing IBAMR](#installing-ibamr)
- [Installing MATLAB (for mesh generation)](#installing-matlab)
- [Building the Project](#building-the-project)
- [Verification](#verification)
- [Troubleshooting](#troubleshooting)

## Prerequisites

### Required Software

- **C++ Compiler**: GCC 8.0+ or Clang 10.0+ with C++11 support
- **CMake**: Version 3.15.0 or higher
- **MPI**: OpenMPI, MPICH, or Intel MPI
- **IBAMR**: Version 0.10.0 or higher (installation instructions below)
- **MATLAB**: R2018b or higher (for mesh generation)
- **VisIt**: For visualization (optional but recommended)

### System Requirements

- **Memory**: Minimum 8 GB RAM (16 GB recommended)
- **Disk Space**: ~5-10 GB for dependencies and simulation output
- **Cores**: Multi-core processor recommended (4+ cores)

## Installing IBAMR

IBAMR (Immersed Boundary Adaptive Mesh Refinement) is the primary dependency for this project.

### Option 1: Install from Spack (Recommended)

Spack is the easiest way to install IBAMR with all dependencies:

```bash
# Install Spack
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
source spack/share/spack/setup-env.sh

# Add Spack to your shell profile
echo 'source /path/to/spack/share/spack/setup-env.sh' >> ~/.bashrc

# Install IBAMR (this will take 1-3 hours)
spack install ibamr@develop
spack load ibamr
```

### Option 2: Build from Source

#### Step 1: Install PETSc

```bash
# Download PETSc
git clone -b release https://gitlab.com/petsc/petsc.git
cd petsc

# Configure PETSc
./configure \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --download-hypre \
  --download-fblaslapack \
  --with-debugging=0 \
  COPTFLAGS='-O3 -march=native' \
  CXXOPTFLAGS='-O3 -march=native' \
  FOPTFLAGS='-O3 -march=native'

# Build and install
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt check

# Set environment variables
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
```

#### Step 2: Install libMesh

```bash
# Download libMesh
git clone https://github.com/libMesh/libmesh.git
cd libmesh

# Configure
./configure \
  --prefix=/usr/local/libmesh \
  --with-methods="opt" \
  --enable-petsc \
  --disable-eigen

# Build and install
make -j4
make install
```

#### Step 3: Install SAMRAI

```bash
# Download SAMRAI
git clone https://github.com/LLNL/SAMRAI.git
cd SAMRAI

# Configure
mkdir build && cd build
cmake \
  -DCMAKE_INSTALL_PREFIX=/usr/local/samrai \
  -DENABLE_HDF5=ON \
  -DHDF5_DIR=/path/to/hdf5 \
  ..

# Build and install
make -j4
make install
```

#### Step 4: Install IBAMR

```bash
# Clone IBAMR
git clone https://github.com/IBAMR/IBAMR.git
cd IBAMR

# Configure
mkdir build && cd build
cmake \
  -DCMAKE_INSTALL_PREFIX=/usr/local/ibamr \
  -DSAMRAI_ROOT=/usr/local/samrai \
  -DlibMesh_ROOT=/usr/local/libmesh \
  -DPETSC_ROOT=${PETSC_DIR} \
  ..

# Build and install
make -j4
make install

# Set environment variable
export IBAMR_ROOT=/usr/local/ibamr
```

### Ubuntu/Debian Quick Start

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  cmake \
  gfortran \
  libopenmpi-dev \
  openmpi-bin \
  libhdf5-openmpi-dev \
  libboost-all-dev \
  libeigen3-dev \
  git

# Then follow Spack installation (Option 1) for IBAMR
```

### macOS Quick Start

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install cmake gcc open-mpi hdf5 boost eigen

# Then follow Spack installation (Option 1) for IBAMR
```

## Installing MATLAB

MATLAB is required for generating the swimmer mesh.

### Option 1: MATLAB (Commercial)

1. Download MATLAB from [MathWorks](https://www.mathworks.com/downloads/)
2. Install with default settings
3. Ensure MATLAB is in your PATH:
   ```bash
   export PATH=$PATH:/usr/local/MATLAB/R2023a/bin
   ```

### Option 2: GNU Octave (Free Alternative)

Octave can run most MATLAB scripts:

```bash
# Ubuntu/Debian
sudo apt-get install octave

# macOS
brew install octave
```

**Note**: Some plotting features may differ between MATLAB and Octave.

## Building the Project

### Step 1: Clone the Repository

```bash
git clone https://github.com/vinodthale/Naca0012carangiform.git
cd Naca0012carangiform
```

### Step 2: Generate Swimmer Mesh

```bash
# Using MATLAB
matlab -batch "swimming_mode='carangiform'; run('naca0012_swimmer_generator.m')"

# Or using Octave
octave --eval "swimming_mode='carangiform'; naca0012_swimmer_generator"
```

This generates `naca0012carangiform.vertex`.

### Step 3: Configure Build

```bash
mkdir build
cd build

# Configure with CMake
cmake ..

# If IBAMR is not found, specify the path
cmake -DIBAMR_ROOT=/usr/local/ibamr ..
```

### Step 4: Build

```bash
make -j4
```

This creates the executable `main2d`.

### Step 5: Run Simulation

```bash
# Run with 4 MPI processes
mpirun -np 4 ./main2d ../input2d

# Or using the convenience script
cd ..
bash run_simulation.sh
```

## Verification

### Test Mesh Generation

```bash
matlab -batch "swimming_mode='carangiform'; run('naca0012_swimmer_generator.m')"
```

You should see:
- Console output showing mesh statistics
- Three plots: amplitude envelope, backbone motion, mesh visualization
- Output file: `naca0012carangiform.vertex`

### Test Compilation

```bash
cd build
make
```

Should complete without errors and produce `main2d` executable.

### Test Simulation (Short Run)

```bash
# Edit input2d to set shorter simulation time
# Change END_TIME to 0.1 for quick test
mpirun -np 2 ./main2d ../input2d
```

Should run without errors and produce output in:
- `NACA0012Str/` - Lagrangian structure data
- `viz_naca0012_Str/` - Eulerian field data

## Troubleshooting

### IBAMR Not Found

**Error**: `Could not find IBAMR`

**Solution**:
```bash
cmake -DIBAMR_ROOT=/path/to/ibamr/install ..
```

Or set environment variable:
```bash
export IBAMR_ROOT=/path/to/ibamr/install
```

### MPI Compilation Errors

**Error**: `mpicxx: command not found`

**Solution**:
```bash
# Ubuntu/Debian
sudo apt-get install libopenmpi-dev

# macOS
brew install open-mpi
```

### CMake Version Too Old

**Error**: `CMake 3.15 or higher is required`

**Solution**:
```bash
# Download latest CMake from https://cmake.org/download/
# Or use pip
pip install cmake --upgrade
```

### Simulation Divergence

**Error**: Simulation crashes or produces NaN values

**Solutions**:
1. Check grid resolution matches mesh density
2. Reduce time step: `DT_MAX = 0.0005` in `input2d`
3. Verify amplitude envelope is reasonable
4. Check CFL condition: `CFL_MAX = 0.3`

### Memory Issues

**Error**: `Out of memory` or slow performance

**Solutions**:
1. Reduce grid resolution in `input2d`:
   ```
   N = 1024, 512  # Instead of 2048, 1024
   ```
2. Reduce AMR levels:
   ```
   MAX_LEVELS = 2  # Instead of 3
   ```
3. Run on more cores with less memory per process

### Visualization Issues

**Problem**: Cannot open VisIt files

**Solution**:
1. Install VisIt from [visit.llnl.gov](https://visit.llnl.gov)
2. Open the `.visit` file:
   ```bash
   visit -o viz_naca0012_Str/dumps.visit
   ```

## Performance Optimization

### Compiler Flags

For better performance, configure with optimization flags:

```bash
cmake \
  -DCMAKE_CXX_FLAGS="-O3 -march=native" \
  -DCMAKE_C_FLAGS="-O3 -march=native" \
  ..
```

### MPI Tuning

For large simulations on clusters:

```bash
# Use processor binding
mpirun --bind-to core -np 16 ./main2d input2d

# For OpenMPI with Infiniband
mpirun --mca btl openib,self,sm -np 16 ./main2d input2d
```

## Getting Help

If you encounter issues not covered here:

1. Check the [README.md](README.md) for usage information
2. Review [IBAMR documentation](https://ibamr.github.io)
3. Open an issue on GitHub with:
   - Your system information
   - Complete error messages
   - Steps to reproduce

## Next Steps

After successful installation:

1. Read [README.md](README.md) for simulation usage
2. See [MESH_PARAMETER_GUIDE.md](MESH_PARAMETER_GUIDE.md) for mesh customization
3. Review [CONTRIBUTING.md](CONTRIBUTING.md) to contribute
4. Explore different swimming modes (carangiform vs anguilliform)
