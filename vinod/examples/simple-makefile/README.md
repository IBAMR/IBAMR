# Simple Makefile Example

This is a simple IBAMR example that demonstrates basic setup using a traditional Makefile build system.

## Features

- IBAMR/IBTK initialization
- MPI parallel support
- SAMRAI grid geometry setup
- Variable database access
- Parallel I/O operations

## Building

To build this example, you need to have IBAMR installed. Set the `IBAMR_PREFIX` environment variable to point to your IBAMR installation:

```bash
export IBAMR_PREFIX=/path/to/ibamr/install
make
```

Or specify it directly:

```bash
make IBAMR_PREFIX=/path/to/ibamr/install
```

## Running

After building, run the example:

```bash
make run
```

Or with multiple MPI processes:

```bash
make run-parallel
```

Or manually:

```bash
mpirun -np 4 ./simple-makefile-example
```

## Cleaning

To clean build artifacts:

```bash
make clean
```

## Expected Output

The program should print:
- MPI rank and process information
- SAMRAI dimension (NDIM)
- Variable database information
- Grid geometry setup details
- Periodic boundary condition settings
