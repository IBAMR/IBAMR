#!/bin/bash
# Script to build all IBAMR C++ tests

set -e  # Exit on error

echo "=========================================="
echo "Building IBAMR C++ Test Suite"
echo "=========================================="

# Check if IBAMR_ROOT is set
if [ -z "$IBAMR_ROOT" ]; then
    echo "ERROR: IBAMR_ROOT environment variable not set"
    echo "Please set IBAMR_ROOT to your IBAMR installation directory:"
    echo "  export IBAMR_ROOT=/path/to/ibamr"
    exit 1
fi

echo "IBAMR_ROOT: $IBAMR_ROOT"
echo ""

# Create build directory
if [ -d "build" ]; then
    echo "Build directory exists. Removing..."
    rm -rf build
fi

mkdir build
cd build

echo "Running CMake..."
cmake .. -DIBAMR_ROOT=$IBAMR_ROOT

echo ""
echo "Building tests..."
make -j4

echo ""
echo "=========================================="
echo "Build complete!"
echo "=========================================="
echo ""
echo "Executables created:"
ls -1 test* 2>/dev/null || echo "No test executables found"

echo ""
echo "To run a specific test:"
echo "  cd Test01_SmokeTest"
echo "  ../build/test01_smoke input2d"
echo ""
echo "Or use run_all_tests.sh to run all tests"
