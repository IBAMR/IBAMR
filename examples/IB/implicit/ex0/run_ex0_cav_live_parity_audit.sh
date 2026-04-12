#!/usr/bin/env bash
set -euo pipefail

IBAMR_BIN=${IBAMR_BIN:-/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/main2d}
IBAMR_INPUT_REF=${IBAMR_INPUT_REF:-/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d.reference}
IBAMR_PETSC_REF=${IBAMR_PETSC_REF:-/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/petsc_options.reference.dat}
MATLAB_BIN=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
MATLAB_SCRIPT=${MATLAB_SCRIPT:-/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_export_matlab.m}
COMPARE_SCRIPT=${COMPARE_SCRIPT:-/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_compare.py}

OUTPUT_ROOT=${OUTPUT_ROOT:-/tmp/ibamr_cav_live_parity}
CASE_ID=${CASE_ID:-ex0_live_parity}
ABS_TOL=${ABS_TOL:-1.0e-12}
REL_TOL=${REL_TOL:-1.0e-10}
DT_FACTOR=${DT_FACTOR:-0.5}
ELASTIC_K=${ELASTIC_K:-1.0e8}
USE_MATRIX_BASED_SAJ=${USE_MATRIX_BASED_SAJ:-}
STAGE_D_PRESSURE_GAUGE_PROJECTED=${STAGE_D_PRESSURE_GAUGE_PROJECTED:-true}
STAGE_D_COND_SAFETY_FACTOR=${STAGE_D_COND_SAFETY_FACTOR:-1.0}
STAGE_D_COND_THRESHOLD_FLOOR=${STAGE_D_COND_THRESHOLD_FLOOR:-1.0e-12}
STAGE_D_COND2_BACKEND=${STAGE_D_COND2_BACKEND:-matlab}
STAGE_D_COND2_MATLAB_BIN=${STAGE_D_COND2_MATLAB_BIN:-$MATLAB_BIN}

COARSE_N=${COARSE_N:-8}
MAX_LEVELS=${MAX_LEVELS:-3}
N_FINE=$((COARSE_N * (2 ** (MAX_LEVELS - 1))))
if (( N_FINE < 32 )); then
  echo "error: finest grid must be at least 32x32 (got N_FINE=$N_FINE from COARSE_N=$COARSE_N MAX_LEVELS=$MAX_LEVELS)" >&2
  exit 2
fi
DX_FINE=$(python3 - <<PY
N=$N_FINE
print(f"{1.0/float(N):.17g}")
PY
)

mkdir -p "$OUTPUT_ROOT/reports"
echo "grid_config: COARSE_N=$COARSE_N MAX_LEVELS=$MAX_LEVELS N_FINE=$N_FINE ELASTIC_K=$ELASTIC_K"
if [[ -n "$USE_MATRIX_BASED_SAJ" ]]; then
  echo "solver_config: USE_MATRIX_BASED_SAJ=$USE_MATRIX_BASED_SAJ"
fi

prepare_policy_input() {
  local input_file=$1
  cp "$IBAMR_INPUT_REF" "$input_file"
  sed -i '' 's/#.*$//' "$input_file"
  sed -i '' -E "s/^MAX_LEVELS[[:space:]]*=.*/MAX_LEVELS = $MAX_LEVELS/" "$input_file"
  sed -i '' -E "s/^N[[:space:]]*=.*/N = $COARSE_N/" "$input_file"
  sed -i '' -E "s/^K[[:space:]]*=.*/K = $ELASTIC_K/" "$input_file"
  if [[ -n "$USE_MATRIX_BASED_SAJ" ]]; then
    if rg -q '^USE_MATRIX_BASED_SAJ[[:space:]]*=' "$input_file"; then
      sed -i '' -E "s/^USE_MATRIX_BASED_SAJ[[:space:]]*=.*/USE_MATRIX_BASED_SAJ = $USE_MATRIX_BASED_SAJ/" "$input_file"
    else
      printf "\nUSE_MATRIX_BASED_SAJ = %s\n" "$USE_MATRIX_BASED_SAJ" >> "$input_file"
    fi
  fi
}

write_parity_block() {
  local input_file=$1
  local policy=$2
  local export_preconditioned=$3

  cat >> "$input_file" <<EOF2

parity_audit_db {
  enabled = TRUE
  output_dir = "$OUTPUT_ROOT/ibamr"
  case_id = "$CASE_ID"
  closure_policy = "$policy"
  export_preconditioned_operator = $export_preconditioned
}
EOF2
}

run_ibamr_export() {
  local policy=$1
  local export_preconditioned=$2
  local input_file=$3
  local ibamr_dir=$4

  prepare_policy_input "$input_file"
  write_parity_block "$input_file" "$policy" "$export_preconditioned"
  "$IBAMR_BIN" "$input_file" > "$ibamr_dir/run_ibamr.stdout" 2>&1
}

run_matlab_export() {
  local policy=$1
  local export_preconditioned=$2
  local matlab_dir=$3

  "$MATLAB_BIN" -batch "\
setenv('CAV_MATLAB_PARITY_DIR','$matlab_dir'); \
setenv('CAV_CLOSURE_POLICY','$policy'); \
setenv('CAV_CASE_ID','$CASE_ID'); \
setenv('CAV_COARSE_N','$COARSE_N'); \
setenv('CAV_MAX_LEVELS','$MAX_LEVELS'); \
setenv('CAV_X_CENTER','0.5'); \
setenv('CAV_Y_CENTER','0.5'); \
setenv('CAV_X_RADIUS','0.2'); \
setenv('CAV_Y_RADIUS','0.2'); \
setenv('CAV_DS','$DX_FINE'); \
setenv('CAV_DT_FACTOR','$DT_FACTOR'); \
setenv('CAV_K','$ELASTIC_K'); \
setenv('CAV_NU1','1'); \
setenv('CAV_NU2','1'); \
setenv('CAV_BOX_SIZE','4'); \
setenv('CAV_OVERLAP','2'); \
setenv('CAV_RELAX_ALPHA','1.0'); \
setenv('CAV_EXPORT_PRECONDITIONED_OPERATOR','$export_preconditioned'); \
run('$MATLAB_SCRIPT');" > "$matlab_dir/run_matlab.stdout" 2>&1
}

run_policy() {
  local policy=$1
  local policy_lower
  policy_lower=$(echo "$policy" | tr 'A-Z' 'a-z')

  local input_dir="$OUTPUT_ROOT/inputs/$CASE_ID/policy=$policy"
  local ibamr_dir="$OUTPUT_ROOT/ibamr/$CASE_ID/policy=$policy"
  local matlab_dir="$OUTPUT_ROOT/matlab/$CASE_ID/policy=$policy"
  local report_json="$OUTPUT_ROOT/reports/${CASE_ID}.${policy_lower}.report.json"
  local report_txt="$OUTPUT_ROOT/reports/${CASE_ID}.${policy_lower}.report.txt"
  local report_json_stage_c="$OUTPUT_ROOT/reports/${CASE_ID}.${policy_lower}.stage_c.report.json"
  local report_txt_stage_c="$OUTPUT_ROOT/reports/${CASE_ID}.${policy_lower}.stage_c.report.txt"
  local input_file="$input_dir/input2d.live_parity"

  mkdir -p "$input_dir" "$ibamr_dir" "$matlab_dir"
  cp "$IBAMR_PETSC_REF" "$input_dir/petsc_options.reference.dat"

  run_ibamr_export "$policy" "FALSE" "$input_file" "$ibamr_dir"
  run_matlab_export "$policy" "FALSE" "$matlab_dir"
  if ! "$COMPARE_SCRIPT" \
    --ibamr-dir "$ibamr_dir" \
    --matlab-dir "$matlab_dir" \
    --max-stage "C" \
    --abs-tol "$ABS_TOL" \
    --rel-tol "$REL_TOL" \
    --report-json "$report_json_stage_c" \
    --report-text "$report_txt_stage_c"; then
    :
  fi

  local stage_c_status="FAIL"
  if [[ -f "$report_json_stage_c" ]]; then
    stage_c_status=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "$report_json_stage_c")
  fi
  if [[ "$stage_c_status" != "PASS" ]]; then
    if [[ -f "$report_json_stage_c" ]]; then
      cp "$report_json_stage_c" "$report_json"
    fi
    if [[ -f "$report_txt_stage_c" ]]; then
      cp "$report_txt_stage_c" "$report_txt"
    fi
    echo "policy=$policy report=$report_txt"
    return 1
  fi

  run_ibamr_export "$policy" "TRUE" "$input_file" "$ibamr_dir"
  run_matlab_export "$policy" "TRUE" "$matlab_dir"
  if ! "$COMPARE_SCRIPT" \
    --ibamr-dir "$ibamr_dir" \
    --matlab-dir "$matlab_dir" \
    --max-stage "D" \
    --abs-tol "$ABS_TOL" \
    --rel-tol "$REL_TOL" \
    --stage-d-pressure-gauge-projected "$STAGE_D_PRESSURE_GAUGE_PROJECTED" \
    --stage-d-cond-safety-factor "$STAGE_D_COND_SAFETY_FACTOR" \
    --stage-d-cond-threshold-floor "$STAGE_D_COND_THRESHOLD_FLOOR" \
    --stage-d-cond2-backend "$STAGE_D_COND2_BACKEND" \
    --stage-d-cond2-matlab-bin "$STAGE_D_COND2_MATLAB_BIN" \
    --report-json "$report_json" \
    --report-text "$report_txt"; then
    echo "policy=$policy report=$report_txt"
    return 1
  fi

  echo "policy=$policy report=$report_txt"
}

status=0
run_policy RELAXED || status=1
run_policy STRICT || status=1

echo "live parity audit complete: $OUTPUT_ROOT"
exit "$status"
