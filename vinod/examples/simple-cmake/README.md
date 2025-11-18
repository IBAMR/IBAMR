# Simple CMake Example

This is a simple IBAMR example that demonstrates basic setup and initialization using CMake.

## Features

- IBAMR/IBTK initialization
- SAMRAI variable database access
- Eigen library integration
- Optional libMesh support

## Building

To build this example, you need to have IBAMR installed. Then:

```bash
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/ibamr/install ..
make
```

## Running

After building, run the example:

```bash
./simple-cmake-example
```

Or with MPI:

```bash
mpirun -np 4 ./simple-cmake-example
```

## Expected Output

The program should print:
- SAMRAI variable database information
- Eigen version
- A simple matrix trace calculation
- libMesh availability status
