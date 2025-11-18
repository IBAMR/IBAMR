# IBAMR Custom Examples

This directory contains custom IBAMR examples for testing and demonstration purposes.

## Directory Structure

```
vinod/examples/
├── README.md                              # This file
├── anguilliform/                          # Anguilliform swimming kinematics
├── carangiform/                           # Carangiform swimming kinematics
├── EelBAMRvinod/                         # Eel swimming simulation
├── Effect-of-shape-and-adaptive-kinematics/ # Fish shape and kinematics study
├── Four_fish_school/                     # Multi-fish schooling simulation
├── Naca0012carangiform/                  # NACA0012 airfoil with carangiform motion
├── Zhang_2018/                           # Zhang 2018 validation cases
├── simple-cmake/                         # CMake-based simple example
└── simple-makefile/                      # Makefile-based simple example
```

## Examples

### Basic Examples

#### simple-cmake
A simple IBAMR example demonstrating CMake-based build system integration. Shows:
- Basic IBAMR/IBTK initialization
- SAMRAI variable database usage
- Eigen library integration
- Optional libMesh support

See [simple-cmake/README.md](simple-cmake/README.md) for details.

#### simple-makefile
A simple IBAMR example using a traditional Makefile. Demonstrates:
- IBAMR initialization with Makefile build
- MPI parallel support
- SAMRAI grid geometry setup
- Parallel I/O operations

See [simple-makefile/README.md](simple-makefile/README.md) for details.

### Swimming Kinematics Examples

#### anguilliform
Anguilliform (eel-like) swimming kinematics implementation.
- Body wave propagation patterns
- High amplitude tail motion
- Efficient low-speed swimming

#### carangiform
Carangiform (tuna-like) swimming kinematics implementation.
- Posterior body flexure
- High-speed cruising
- Reduced drag coefficients

#### EelBAMRvinod
Complete eel swimming simulation with:
- Custom IBEELKinematics class
- 2D eel geometry generation
- Multiple swimming modes

See [EelBAMRvinod/README.md](EelBAMRvinod/README.md) for details.

### Advanced Examples

#### Effect-of-shape-and-adaptive-kinematics
Fish-inspired hydrodynamic study examining:
- Effects of body shape on swimming performance
- Adaptive kinematics optimization
- Reynolds number effects
- Parameter sweep studies

See [Effect-of-shape-and-adaptive-kinematics/README.md](Effect-of-shape-and-adaptive-kinematics/README.md) for details.

#### Four_fish_school
Multi-fish schooling simulation demonstrating:
- 4 swimming eels in rectangular formation
- Coupled Navier-Stokes + Immersed Boundary + Advection-Diffusion
- Odor plume transport and vortex dynamics
- Multi-structure force tracking
- Production-ready C++ test suite

See [Four_fish_school/README.md](Four_fish_school/README.md) for details.

#### Naca0012carangiform
NACA0012 airfoil with carangiform swimming motion:
- Airfoil geometry generation
- Undulatory propulsion mechanics
- Hydrodynamic performance analysis
- Thrust and efficiency calculations

See [Naca0012carangiform/README.md](Naca0012carangiform/README.md) for details.

### Validation Examples

#### Zhang_2018
Validation cases from Zhang et al. (2018):
- Undulating foil hydrodynamics
- Reynolds number effects (Re=50, Re=1000)
- Thickness ratio variations
- Performance metrics validation

See [Zhang_2018/README.md](Zhang_2018/README.md) for details.

## Building Examples

### Prerequisites

- IBAMR installed (either via autoibamr or manual build)
- CMake 3.15 or later (for CMake examples)
- MPI compiler (mpic++)
- Standard build tools (make, gcc/g++)

### Build Instructions

#### CMake-based Examples

```bash
cd <example-directory>
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/ibamr/install ..
make
```

#### Makefile-based Examples

```bash
cd <example-directory>
make IBAMR_PREFIX=/path/to/ibamr/install
```

## Running Examples

Most examples include input files (typically named `input2d` or `input3d`). Run with:

```bash
# Serial execution
./main2d input2d

# Parallel execution
mpirun -np 4 ./main2d input2d
```

## CI/CD Integration

These examples are automatically built and tested by the GitHub Actions workflow:
`.github/workflows/ibamr-full-build.yml`

The workflow:
1. Builds IBAMR from scratch using autoibamr
2. Compiles all examples in this directory
3. Runs static analysis (clang-tidy, cppcheck)
4. Uploads build artifacts and logs

## Adding New Examples

To add a new example:

1. Create a new subdirectory under `vinod/examples/`
2. Add either a `CMakeLists.txt` or `Makefile`
3. Include source files and a README.md
4. Document the example purpose and usage
5. The CI workflow will automatically detect and build it

## Troubleshooting

### CMake can't find IBAMR

Make sure to set `CMAKE_PREFIX_PATH` to your IBAMR installation directory:
```bash
cmake -DCMAKE_PREFIX_PATH=/path/to/ibamr/install ..
```

### Makefile build errors

Ensure `IBAMR_PREFIX` points to the correct installation:
```bash
export IBAMR_PREFIX=/path/to/ibamr/install
make
```

### MPI errors

Ensure your MPI installation is compatible with the MPI used to build IBAMR.

### Missing geometry files

Some examples require MATLAB/Octave to generate geometry files (.vertex files).
Check the example README for specific requirements.

## License

These examples are part of IBAMR and are distributed under the 3-clause BSD license.
See the COPYRIGHT file at the top level of the IBAMR repository.
