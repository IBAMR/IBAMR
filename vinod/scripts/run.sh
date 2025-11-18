#!/bin/bash
# Run script for IBAMR simulations
# Usage: ./run.sh [example_name] [num_procs] [input_file]

set -e  # Exit on error

# Default values
EXAMPLE=${1:-"simple_cylinder"}
NUM_PROCS=${2:-4}
INPUT_FILE=${3:-"input2d"}

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}IBAMR Run Script${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "Example: ${YELLOW}$EXAMPLE${NC}"
echo -e "Processors: ${YELLOW}$NUM_PROCS${NC}"
echo -e "Input File: ${YELLOW}$INPUT_FILE${NC}"
echo ""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VINOD_DIR="$(dirname "$SCRIPT_DIR")"
EXAMPLE_DIR="$VINOD_DIR/examples/$EXAMPLE"

# Check if example exists
if [ ! -d "$EXAMPLE_DIR" ]; then
    echo -e "${RED}Error: Example directory not found: $EXAMPLE_DIR${NC}"
    exit 1
fi

cd "$EXAMPLE_DIR"

# Determine dimension from input file name
if [[ "$INPUT_FILE" == *"2d"* ]]; then
    EXECUTABLE="main2d"
    BUILD_DIR="build_2d"
elif [[ "$INPUT_FILE" == *"3d"* ]]; then
    EXECUTABLE="main3d"
    BUILD_DIR="build_3d"
else
    echo -e "${YELLOW}Warning: Cannot determine dimension from input file name${NC}"
    echo -e "${YELLOW}Assuming 2D...${NC}"
    EXECUTABLE="main2d"
    BUILD_DIR="build_2d"
fi

# Check if executable exists
if [ ! -f "$BUILD_DIR/$EXECUTABLE" ] && [ ! -f "$EXECUTABLE" ]; then
    echo -e "${RED}Error: Executable not found: $EXECUTABLE${NC}"
    echo -e "${YELLOW}Please build first using: ./scripts/build.sh $EXAMPLE${NC}"
    exit 1
fi

# Use executable from build directory if it exists, otherwise current directory
if [ -f "$BUILD_DIR/$EXECUTABLE" ]; then
    EXEC_PATH="$BUILD_DIR/$EXECUTABLE"
else
    EXEC_PATH="./$EXECUTABLE"
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

# Clean old output directories (optional, comment out if you want to keep them)
# rm -rf viz_* restart_* postproc_*

echo -e "${GREEN}Starting simulation...${NC}"
echo -e "${GREEN}Command: ${NC}mpirun -np $NUM_PROCS $EXEC_PATH $INPUT_FILE"
echo ""

# Run simulation
mpirun -np $NUM_PROCS $EXEC_PATH $INPUT_FILE 2>&1 | tee simulation_output.log

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Simulation complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Output log: ${YELLOW}simulation_output.log${NC}"
echo -e "Visualization data: ${YELLOW}viz_*/${NC}"
echo ""
echo -e "To visualize: ${YELLOW}visit -o viz_*/dumps.visit${NC}"
echo ""
