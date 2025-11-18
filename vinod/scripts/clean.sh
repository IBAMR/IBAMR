#!/bin/bash
# Clean script to remove build artifacts and output data
# Usage: ./clean.sh [option]
# Options:
#   build    - Clean only build directories
#   output   - Clean only simulation output
#   all      - Clean everything (default)

OPTION=${1:-"all"}

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VINOD_DIR="$(dirname "$SCRIPT_DIR")"

echo -e "${YELLOW}========================================${NC}"
echo -e "${YELLOW}IBAMR Clean Script${NC}"
echo -e "${YELLOW}========================================${NC}"

clean_build() {
    echo -e "${GREEN}Cleaning build directories...${NC}"
    find "$VINOD_DIR/examples" -type d -name "build_*" -exec rm -rf {} + 2>/dev/null
    find "$VINOD_DIR/examples" -type f -name "main2d" -delete 2>/dev/null
    find "$VINOD_DIR/examples" -type f -name "main3d" -delete 2>/dev/null
    find "$VINOD_DIR/examples" -type f -name "*.o" -delete 2>/dev/null
    echo -e "${GREEN}Build directories cleaned${NC}"
}

clean_output() {
    echo -e "${GREEN}Cleaning simulation output...${NC}"
    find "$VINOD_DIR/examples" -type d -name "viz_*" -exec rm -rf {} + 2>/dev/null
    find "$VINOD_DIR/examples" -type d -name "restart_*" -exec rm -rf {} + 2>/dev/null
    find "$VINOD_DIR/examples" -type d -name "postproc_*" -exec rm -rf {} + 2>/dev/null
    find "$VINOD_DIR/examples" -type f -name "*.log" -delete 2>/dev/null
    find "$VINOD_DIR/examples" -type f -name "simulation_output.log" -delete 2>/dev/null
    find "$VINOD_DIR/examples" -type f -name "*.png" -delete 2>/dev/null
    echo -e "${GREEN}Simulation output cleaned${NC}"
}

case "$OPTION" in
    build)
        clean_build
        ;;
    output)
        clean_output
        ;;
    all)
        clean_build
        clean_output
        ;;
    *)
        echo -e "${RED}Error: Unknown option '$OPTION'${NC}"
        echo "Usage: $0 [build|output|all]"
        exit 1
        ;;
esac

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Cleaning complete!${NC}"
echo -e "${GREEN}========================================${NC}"
