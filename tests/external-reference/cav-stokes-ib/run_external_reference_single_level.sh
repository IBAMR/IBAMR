#!/bin/bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <build_dir>"
  exit 1
fi

BUILD_DIR="$1"
TEST_DIR="$BUILD_DIR/tests/IB"
INPUT_FILE="$BUILD_DIR/tests/external-reference/cav-stokes-ib/implicit_stokes_ib_reference_single_level_01.input"
REF_DIR="$BUILD_DIR/tests/external-reference/cav-stokes-ib/reference-generated/single_level_case_01"
REF_FILE="$REF_DIR/X.ref"

if [ ! -f "$REF_FILE" ]; then
  echo "Missing local reference artifacts:"
  echo "  $REF_FILE"
  echo "Generate them first with:"
  echo "  tests/external-reference/cav-stokes-ib/reference-tools/generate_single_level_case_01.sh <build_dir> <path_to_implicit_ib_coupling_aware_vanka>"
  exit 1
fi

cd "$TEST_DIR"
TMP_INPUT="${INPUT_FILE%.input}.run_external_reference.input"
cp "$INPUT_FILE" "$TMP_INPUT"

TMP_INPUT_ENV="$TMP_INPUT" REF_DIR_ENV="$REF_DIR" python3 - <<'PY'
import os
import re
from pathlib import Path
p = Path(os.environ["TMP_INPUT_ENV"])
s = p.read_text()
s = re.sub(r"write_reference\s*=\s*TRUE", "write_reference = FALSE", s)
s = re.sub(r'reference_dir\s*=\s*"[^"]*"', f'reference_dir = "{os.environ["REF_DIR_ENV"]}"', s)
p.write_text(s)
PY

./implicit_stokes_ib_operator_chain_parity_01 "$TMP_INPUT"
rm -f "$TMP_INPUT"
