#!/bin/bash
# Build script for IBAMR examples and custom code
# Usage: ./build.sh [example_name] [2d|3d] [Debug|Release]

set -e  # Exit on error

# Default values
EXAMPLE=${1:-"simple_cylinder"}
DIMENSION=${2:-"2d"}
BUILD_TYPE=${3:-"Release"}

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}IBAMR Build Script${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Example: ${YELLOW}$EXAMPLE${NC}"
echo -e "Dimension: ${YELLOW}$DIMENSION${NC}"
echo -e "Build Type: ${YELLOW}$BUILD_TYPE${NC}"
echo ""

# Determine NDIM
if [ "$DIMENSION" == "2d" ]; then
    NDIM=2
elif [ "$DIMENSION" == "3d" ]; then
    NDIM=3
else
    echo -e "${RED}Error: Invalid dimension '$DIMENSION'. Use '2d' or '3d'${NC}"
    exit 1
fi

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VINOD_DIR="$(dirname "$SCRIPT_DIR")"
EXAMPLE_DIR="$VINOD_DIR/examples/$EXAMPLE"

# Check if example exists
if [ ! -d "$EXAMPLE_DIR" ]; then
    echo -e "${RED}Error: Example directory not found: $EXAMPLE_DIR${NC}"
    echo "Available examples:"
    ls -1 "$VINOD_DIR/examples/"
    exit 1
fi

echo -e "${GREEN}Building in: ${NC}$EXAMPLE_DIR"
cd "$EXAMPLE_DIR"

# Create build directory
BUILD_DIR="build_$DIMENSION"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

echo -e "${GREEN}Configuring with CMake...${NC}"
cmake .. \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DIBAMR_DIMENSIONS=$NDIM \
    || { echo -e "${RED}CMake configuration failed${NC}"; exit 1; }

echo -e "${GREEN}Building...${NC}"
make -j$(nproc) || { echo -e "${RED}Build failed${NC}"; exit 1; }

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Build successful!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Executable: ${YELLOW}$BUILD_DIR/main${DIMENSION:0:1}d${NC}"
echo -e "To run: ${YELLOW}cd $EXAMPLE_DIR && mpirun -np 4 ./$BUILD_DIR/main${DIMENSION:0:1}d input${DIMENSION}${NC}"
echo ""
