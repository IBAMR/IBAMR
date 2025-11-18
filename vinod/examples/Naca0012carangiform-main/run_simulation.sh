#!/bin/bash
# Simple IBAMR simulation runner

# Number of MPI processes
NPROCS=6

# Input file
INPUT_FILE="input2d"

# Executable
EXECUTABLE="./build/main2d"

# Run simulation
echo "=========================================="
echo "Starting NACA0012 Carangiform Simulation"
echo "=========================================="
echo "Processors: $NPROCS"
echo "Input file: $INPUT_FILE"
echo "Time: $(date)"
echo "=========================================="

mpirun -np $NPROCS $EXECUTABLE $INPUT_FILE

echo ""
echo "Simulation completed at: $(date)"