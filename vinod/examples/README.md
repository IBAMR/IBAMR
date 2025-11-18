# IBAMR Custom Examples

This directory contains custom IBAMR examples for testing and demonstration purposes.

## Directory Structure

```
vinod/examples/
├── README.md                  # This file
├── simple-cmake/              # CMake-based example
│   ├── CMakeLists.txt
│   ├── main.cpp
│   └── README.md
└── simple-makefile/           # Makefile-based example
    ├── Makefile
    ├── main.cpp
    └── README.md
```

## Examples

### simple-cmake

A simple IBAMR example demonstrating CMake-based build system integration. This example shows:
- Basic IBAMR/IBTK initialization
- SAMRAI variable database usage
- Eigen library integration
- Optional libMesh support

See [simple-cmake/README.md](simple-cmake/README.md) for build and run instructions.

### simple-makefile

A simple IBAMR example using a traditional Makefile. This example demonstrates:
- IBAMR initialization with Makefile build
- MPI parallel support
- SAMRAI grid geometry setup
- Parallel I/O operations

See [simple-makefile/README.md](simple-makefile/README.md) for build and run instructions.

## Building All Examples

### Prerequisites

- IBAMR installed (either via autoibamr or manual build)
- CMake 3.15 or later (for CMake examples)
- MPI compiler (mpic++)
- Standard build tools (make, gcc/g++)

### Build Instructions

#### CMake Example

```bash
cd simple-cmake
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/ibamr/install ..
make
./simple-cmake-example
```

#### Makefile Example

```bash
cd simple-makefile
make IBAMR_PREFIX=/path/to/ibamr/install
./simple-makefile-example
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
4. The CI workflow will automatically detect and build it

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

## License

These examples are part of IBAMR and are distributed under the 3-clause BSD license.
See the COPYRIGHT file at the top level of the IBAMR repository.
