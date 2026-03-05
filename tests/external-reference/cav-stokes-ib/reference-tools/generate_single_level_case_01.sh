#!/bin/bash
set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <build_dir> <vanka_repo> [input_file] [reference_dir]"
  exit 1
fi

BUILD_DIR="$1"
VANKA_REPO="$2"
INPUT_FILE="${3:-$BUILD_DIR/tests/external-reference/cav-stokes-ib/implicit_stokes_ib_reference_single_level_01.input}"
REF_DIR="${4:-$BUILD_DIR/tests/external-reference/cav-stokes-ib/reference-generated/single_level_case_01}"
TEST_DIR="$BUILD_DIR/tests/IB"
TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)"

cd "$TEST_DIR"
mkdir -p "$REF_DIR"

TMP_INPUT="${INPUT_FILE%.input}.write_reference.input"
cp "$INPUT_FILE" "$TMP_INPUT"

TMP_INPUT_ENV="$TMP_INPUT" REF_DIR_ENV="$REF_DIR" python3 - <<'PY'
import os
import re
from pathlib import Path
p = Path(os.environ["TMP_INPUT_ENV"])
s = p.read_text()
s = re.sub(r"write_reference\s*=\s*FALSE", "write_reference = TRUE", s)
s = re.sub(r'reference_dir\s*=\s*"[^"]*"', f'reference_dir = "{os.environ["REF_DIR_ENV"]}"', s)
p.write_text(s)
PY

./implicit_stokes_ib_operator_chain_parity_01 "$TMP_INPUT"
rm -f "$TMP_INPUT"

octave --quiet --eval "addpath('$TOOLS_DIR'); generate_single_level_case_01('$VANKA_REPO', '$REF_DIR');"

echo "Reference files regenerated from Octave baseline under $REF_DIR"
