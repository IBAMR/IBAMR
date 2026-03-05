#!/bin/bash
set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <build_dir> <vanka_repo>"
  exit 1
fi

BUILD_DIR="$1"
VANKA_REPO="$2"
BASE_DIR="$BUILD_DIR/tests/external-reference/cav-stokes-ib"
TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)"
GEN_ONE="$TOOLS_DIR/generate_single_level_case_01.sh"

bash "$GEN_ONE" "$BUILD_DIR" "$VANKA_REPO" \
  "$BASE_DIR/implicit_stokes_ib_reference_single_level_01.input" \
  "$BASE_DIR/reference-generated/single_level_case_01"

bash "$GEN_ONE" "$BUILD_DIR" "$VANKA_REPO" \
  "$BASE_DIR/implicit_stokes_ib_reference_single_level_01.multilevel_l2.input" \
  "$BASE_DIR/reference-generated/single_level_case_01.multilevel_l2"

bash "$GEN_ONE" "$BUILD_DIR" "$VANKA_REPO" \
  "$BASE_DIR/implicit_stokes_ib_reference_single_level_01.multilevel_l3.input" \
  "$BASE_DIR/reference-generated/single_level_case_01.multilevel_l3"
