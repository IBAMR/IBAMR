#!/bin/bash
#
# Parameter Study Runner for Undulatory Foil Simulations
# Runs simulations across different Reynolds numbers and thicknesses
#
# Usage: ./run_parameter_study.sh [num_processes]
#

# Number of MPI processes (default: 6)
NPROCS=${1:-6}

# Executable
EXEC="./build/main2d"

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: Executable $EXEC not found!"
    echo "Please compile the code first:"
    echo "  mkdir build && cd build && cmake .. && make -j4 && cd .."
    exit 1
fi

echo "========================================"
echo "Undulatory Foil Parameter Study"
echo "========================================"
echo "MPI Processes: $NPROCS"
echo "Executable: $EXEC"
echo "========================================"
echo ""

# List of input files for parameter study
INPUT_FILES=(
    "input2d_Re1000_h004"
    "input2d_Re5609_h006"
    "input2d_Re10000_h008"
)

# Run each simulation
for INPUT in "${INPUT_FILES[@]}"; do
    if [ -f "$INPUT" ]; then
        echo "========================================"
        echo "Running: $INPUT"
        echo "Started: $(date)"
        echo "========================================"

        # Run simulation
        mpirun -np $NPROCS $EXEC $INPUT

        STATUS=$?
        if [ $STATUS -eq 0 ]; then
            echo "✓ Completed successfully: $INPUT"
        else
            echo "✗ Failed with status $STATUS: $INPUT"
        fi

        echo "Finished: $(date)"
        echo ""
    else
        echo "Warning: Input file $INPUT not found, skipping..."
    fi
done

echo "========================================"
echo "All simulations completed!"
echo "========================================"
echo ""
echo "To analyze results, run:"
echo "  python analyze_performance.py"
echo ""
