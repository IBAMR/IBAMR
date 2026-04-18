#!/usr/bin/env bash
set -euo pipefail

IBAMR_BIN=${IBAMR_BIN:-/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/main2d}
IBAMR_INPUT_REF=${IBAMR_INPUT_REF:-/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/input2d.reference}
IBAMR_PETSC_REF=${IBAMR_PETSC_REF:-/Users/boyceg/code/ibamr-objs-opt/examples/IB/implicit/ex0/petsc_options.reference.dat}
MATLAB_BIN=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
MATLAB_SCRIPT=${MATLAB_SCRIPT:-/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_export_matlab.m}
MATLAB_REFERENCE_REPO=${MATLAB_REFERENCE_REPO:-/Users/boyceg/code/implicit_ib_coupling_aware_vanka}
COMPARE_SCRIPT=${COMPARE_SCRIPT:-/Users/boyceg/code/IBAMR/examples/IB/implicit/ex0/run_ex0_cav_live_parity_compare.py}

OUTPUT_ROOT=${OUTPUT_ROOT:-/tmp/ibamr_cav_live_parity}
CASE_ID=${CASE_ID:-ex0_live_parity}
REL_TOL=${REL_TOL:-1.0e-10}
DT_FACTOR=${DT_FACTOR:-0.5}
ELASTIC_K=${ELASTIC_K:-1.0e8}
USE_MATRIX_BASED_SAJ=${USE_MATRIX_BASED_SAJ:-}
STAGE_D_PRESSURE_GAUGE_PROJECTED=${STAGE_D_PRESSURE_GAUGE_PROJECTED:-true}
STAGE_D_COND_SAFETY_FACTOR=${STAGE_D_COND_SAFETY_FACTOR:-5.0}
STAGE_D_COND_THRESHOLD_FLOOR=${STAGE_D_COND_THRESHOLD_FLOOR:-1.0e-12}
STAGE_D_FIXED_THRESHOLD=${STAGE_D_FIXED_THRESHOLD:-1.0e-10}
STAGE_D_COND2_BACKEND=${STAGE_D_COND2_BACKEND:-matlab}
STAGE_D_COND2_MATLAB_BIN=${STAGE_D_COND2_MATLAB_BIN:-$MATLAB_BIN}
STAGE_D_COND2_STRICT_BACKEND=${STAGE_D_COND2_STRICT_BACKEND:-true}
STAGE_D_REQUIRE_COARSE_SEMANTICS=${STAGE_D_REQUIRE_COARSE_SEMANTICS:-true}
STAGE_D_REQUIRE_SMOOTHER_SEMANTICS=${STAGE_D_REQUIRE_SMOOTHER_SEMANTICS:-true}

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

IBAMR_GIT_SHA=$(git -C /Users/boyceg/code/IBAMR rev-parse HEAD 2>/dev/null || echo "unknown")
MATLAB_GIT_SHA=$(git -C "$MATLAB_REFERENCE_REPO" rev-parse HEAD 2>/dev/null || echo "unknown")
MATLAB_V_CYCLE_OVERRIDE_USED="false"
MATLAB_V_CYCLE_OVERRIDE_PATH=""
MATLAB_V_CYCLE_OVERRIDE_SHA256=""
if [[ "$(basename "$MATLAB_SCRIPT")" == "run_ex0_cav_live_parity_export_matlab.m" ]]; then
  script_dir=$(cd "$(dirname "$MATLAB_SCRIPT")" && pwd)
  if [[ -f "$script_dir/v_cycle.m" ]]; then
    MATLAB_V_CYCLE_OVERRIDE_USED="true"
    MATLAB_V_CYCLE_OVERRIDE_PATH="$script_dir/v_cycle.m"
    MATLAB_V_CYCLE_OVERRIDE_SHA256=$(shasum -a 256 "$MATLAB_V_CYCLE_OVERRIDE_PATH" | awk '{print $1}')
  fi
fi

write_runner_failure_report() {
  local report_json=$1
  local report_txt=$2
  local failed_stage=$3
  local message=$4

  python3 - "$report_json" "$report_txt" "$failed_stage" "$message" <<'PY'
import json
import os
import sys

report_json, report_txt, failed_stage, message = sys.argv[1:5]
report = {
    "status": "FAIL",
    "failed_stage": failed_stage,
    "message": message,
    "stages": [{"stage": failed_stage, "ok": False, "message": message}],
}

os.makedirs(os.path.dirname(report_json), exist_ok=True)
with open(report_json, "w", encoding="utf-8") as f:
    json.dump(report, f, indent=2, sort_keys=True)

summary = "\n".join(
    [
        "Status: FAIL",
        f"Failed stage: {failed_stage}",
        f"Message: {message}",
        f"Stage {failed_stage}: FAIL - {message}",
    ]
) + "\n"
with open(report_txt, "w", encoding="utf-8") as f:
    f.write(summary)
PY
}

verify_bundle_artifacts() {
  local stage_letter=$1
  local policy=$2
  local ibamr_dir=$3
  local matlab_dir=$4
  local report_json=$5
  local report_txt=$6

  local validation_output=""
  if ! validation_output=$(python3 - "$ibamr_dir" "$matlab_dir" "$stage_letter" <<'PY'
import glob
import json
import os
import re
import sys

ibamr_dir, matlab_dir, stage_letter = sys.argv[1:4]
errors = []

def parse_levels(root):
    levels = set()
    for path in glob.glob(os.path.join(root, "dof_map_level*.json")):
        match = re.search(r"dof_map_level(\d+)\.json$", os.path.basename(path))
        if match:
            levels.add(int(match.group(1)))
    return sorted(levels)

def inspect_bundle(root, role):
    role_errors = []
    metadata_path = os.path.join(root, "metadata.json")
    finest_level = None
    if not os.path.exists(metadata_path):
        role_errors.append(f"{role}: missing metadata file: {metadata_path}")
    else:
        try:
            with open(metadata_path, "r", encoding="utf-8") as f:
                metadata = json.load(f)
            finest_level = int(metadata["finest_level"])
        except Exception as ex:
            role_errors.append(f"{role}: invalid metadata JSON {metadata_path}: {ex}")

    discovered_levels = parse_levels(root)
    if finest_level is not None:
        expected_levels = list(range(finest_level + 1))
        marker_path = os.path.join(root, f"markers_level{finest_level}.json")
        if not os.path.exists(marker_path):
            role_errors.append(f"{role}: missing marker file: {marker_path}")
    else:
        expected_levels = discovered_levels

    for level in expected_levels:
        dof_path = os.path.join(root, f"dof_map_level{level}.json")
        A_path = os.path.join(root, f"A_level{level}.mtx")
        subdomain_path = os.path.join(root, f"subdomains_level{level}.json")
        if not os.path.exists(dof_path):
            role_errors.append(f"{role}: missing DOF map: {dof_path}")
        if not os.path.exists(A_path):
            role_errors.append(f"{role}: missing level matrix: {A_path}")
        if not os.path.exists(subdomain_path):
            role_errors.append(f"{role}: missing subdomain JSON: {subdomain_path}")
            continue
        try:
            with open(subdomain_path, "r", encoding="utf-8") as f:
                subdomains = json.load(f)
        except Exception as ex:
            role_errors.append(f"{role}: invalid subdomain JSON {subdomain_path}: {ex}")
            continue
        overlap = subdomains.get("overlap")
        if not isinstance(overlap, list):
            role_errors.append(f"{role}: subdomain JSON missing list field 'overlap': {subdomain_path}")
            continue
        for k in range(len(overlap)):
            sub_mat_path = os.path.join(root, f"A_subdomain_level{level}_k{k}.mtx")
            if not os.path.exists(sub_mat_path):
                role_errors.append(f"{role}: missing subdomain matrix: {sub_mat_path}")

    if stage_letter in {"C", "D"} and expected_levels:
        for level in expected_levels[:-1]:
            p_path = os.path.join(root, f"P_level{level}.mtx")
            r_path = os.path.join(root, f"R_level{level}.mtx")
            if not os.path.exists(p_path):
                role_errors.append(f"{role}: missing prolongation matrix: {p_path}")
            if not os.path.exists(r_path):
                role_errors.append(f"{role}: missing restriction matrix: {r_path}")

    if stage_letter == "D":
        for name in [
            "preconditioned_apply_input_level_fine.mtx",
            "preconditioned_apply_output_level_fine.mtx",
            "preconditioned_apply_coarse_rhs_level0.mtx",
            "preconditioned_apply_coarse_correction_level0.mtx",
        ]:
            path = os.path.join(root, name)
            if not os.path.exists(path):
                role_errors.append(f"{role}: missing Stage-D vector artifact: {path}")

    if stage_letter == "D" and role == "ibamr" and expected_levels:
        coarsest_level = min(expected_levels)
        smoother_levels = [level for level in expected_levels if level != coarsest_level]
        for level in smoother_levels:
            first_subdomain_path = os.path.join(root, f"subdomains_first_sweep_level{level}.json")
            if not os.path.exists(first_subdomain_path):
                role_errors.append(
                    f"{role}: missing required first-sweep subdomain JSON for Stage-D semantics: {first_subdomain_path}"
                )
                continue
            try:
                with open(first_subdomain_path, "r", encoding="utf-8") as f:
                    first_subdomains = json.load(f)
            except Exception as ex:
                role_errors.append(f"{role}: invalid first-sweep subdomain JSON {first_subdomain_path}: {ex}")
                continue
            first_overlap = first_subdomains.get("overlap")
            if not isinstance(first_overlap, list):
                role_errors.append(
                    f"{role}: first-sweep subdomain JSON missing list field 'overlap': {first_subdomain_path}"
                )
                continue
            for k in range(len(first_overlap)):
                first_sub_mat_path = os.path.join(root, f"A_subdomain_first_sweep_level{level}_k{k}.mtx")
                if not os.path.exists(first_sub_mat_path):
                    role_errors.append(
                        f"{role}: missing required first-sweep subdomain matrix for Stage-D semantics: {first_sub_mat_path}"
                    )

    return role_errors

errors.extend(inspect_bundle(ibamr_dir, "ibamr"))
errors.extend(inspect_bundle(matlab_dir, "matlab"))
if errors:
    sys.stdout.write("\n".join(errors))
    sys.exit(1)
PY
  ); then
    local ibamr_log="$ibamr_dir/run_ibamr.stdout"
    local matlab_log="$matlab_dir/run_matlab.stdout"
    local failure_message="policy=$policy stage=$stage_letter artifact validation failed before comparator: $validation_output (ibamr_log=$ibamr_log matlab_log=$matlab_log)"
    write_runner_failure_report "$report_json" "$report_txt" "$stage_letter" "$failure_message"
    echo "$failure_message" >&2
    return 1
  fi
}

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
setenv('CAV_MATLAB_REFERENCE_REPO','$MATLAB_REFERENCE_REPO'); \
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

  rm -rf "$ibamr_dir" "$matlab_dir"
  mkdir -p "$input_dir" "$ibamr_dir" "$matlab_dir"
  cp "$IBAMR_PETSC_REF" "$input_dir/petsc_options.reference.dat"

  if ! run_ibamr_export "$policy" "FALSE" "$input_file" "$ibamr_dir"; then
    local failure_message="policy=$policy stage=C IBAMR export failed (ibamr_log=$ibamr_dir/run_ibamr.stdout)"
    write_runner_failure_report "$report_json_stage_c" "$report_txt_stage_c" "C" "$failure_message"
    cp "$report_json_stage_c" "$report_json"
    cp "$report_txt_stage_c" "$report_txt"
    echo "$failure_message" >&2
    echo "policy=$policy report=$report_txt"
    return 1
  fi
  if ! run_matlab_export "$policy" "FALSE" "$matlab_dir"; then
    local failure_message="policy=$policy stage=C MATLAB export failed (matlab_log=$matlab_dir/run_matlab.stdout)"
    write_runner_failure_report "$report_json_stage_c" "$report_txt_stage_c" "C" "$failure_message"
    cp "$report_json_stage_c" "$report_json"
    cp "$report_txt_stage_c" "$report_txt"
    echo "$failure_message" >&2
    echo "policy=$policy report=$report_txt"
    return 1
  fi
  if ! verify_bundle_artifacts "C" "$policy" "$ibamr_dir" "$matlab_dir" "$report_json_stage_c" "$report_txt_stage_c"; then
    cp "$report_json_stage_c" "$report_json"
    cp "$report_txt_stage_c" "$report_txt"
    echo "policy=$policy report=$report_txt_stage_c"
    return 1
  fi
  if ! "$COMPARE_SCRIPT" \
    --ibamr-dir "$ibamr_dir" \
    --matlab-dir "$matlab_dir" \
    --max-stage "C" \
    --rel-tol "$REL_TOL" \
    --normalize-pressure-row-sign "true" \
    --metadata-ibamr-git-sha "$IBAMR_GIT_SHA" \
    --metadata-matlab-git-sha "$MATLAB_GIT_SHA" \
    --metadata-matlab-v-cycle-override-used "$MATLAB_V_CYCLE_OVERRIDE_USED" \
    --metadata-matlab-v-cycle-override-path "$MATLAB_V_CYCLE_OVERRIDE_PATH" \
    --metadata-matlab-v-cycle-override-sha256 "$MATLAB_V_CYCLE_OVERRIDE_SHA256" \
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

  rm -rf "$ibamr_dir" "$matlab_dir"
  mkdir -p "$ibamr_dir" "$matlab_dir"
  if ! run_ibamr_export "$policy" "TRUE" "$input_file" "$ibamr_dir"; then
    local failure_message="policy=$policy stage=D IBAMR export failed (ibamr_log=$ibamr_dir/run_ibamr.stdout)"
    write_runner_failure_report "$report_json" "$report_txt" "D" "$failure_message"
    echo "$failure_message" >&2
    echo "policy=$policy report=$report_txt"
    return 1
  fi
  if ! run_matlab_export "$policy" "TRUE" "$matlab_dir"; then
    local failure_message="policy=$policy stage=D MATLAB export failed (matlab_log=$matlab_dir/run_matlab.stdout)"
    write_runner_failure_report "$report_json" "$report_txt" "D" "$failure_message"
    echo "$failure_message" >&2
    echo "policy=$policy report=$report_txt"
    return 1
  fi
  if ! verify_bundle_artifacts "D" "$policy" "$ibamr_dir" "$matlab_dir" "$report_json" "$report_txt"; then
    echo "policy=$policy report=$report_txt"
    return 1
  fi
  if ! "$COMPARE_SCRIPT" \
    --ibamr-dir "$ibamr_dir" \
    --matlab-dir "$matlab_dir" \
    --max-stage "D" \
    --rel-tol "$REL_TOL" \
    --normalize-pressure-row-sign "true" \
    --stage-d-pressure-gauge-projected "$STAGE_D_PRESSURE_GAUGE_PROJECTED" \
    --stage-d-cond-safety-factor "$STAGE_D_COND_SAFETY_FACTOR" \
    --stage-d-cond-threshold-floor "$STAGE_D_COND_THRESHOLD_FLOOR" \
    --stage-d-fixed-threshold "$STAGE_D_FIXED_THRESHOLD" \
    --stage-d-cond2-backend "$STAGE_D_COND2_BACKEND" \
    --stage-d-cond2-matlab-bin "$STAGE_D_COND2_MATLAB_BIN" \
    --stage-d-cond2-strict-backend "$STAGE_D_COND2_STRICT_BACKEND" \
    --stage-d-require-coarse-semantics "$STAGE_D_REQUIRE_COARSE_SEMANTICS" \
    --stage-d-require-smoother-semantics "$STAGE_D_REQUIRE_SMOOTHER_SEMANTICS" \
    --metadata-ibamr-git-sha "$IBAMR_GIT_SHA" \
    --metadata-matlab-git-sha "$MATLAB_GIT_SHA" \
    --metadata-matlab-v-cycle-override-used "$MATLAB_V_CYCLE_OVERRIDE_USED" \
    --metadata-matlab-v-cycle-override-path "$MATLAB_V_CYCLE_OVERRIDE_PATH" \
    --metadata-matlab-v-cycle-override-sha256 "$MATLAB_V_CYCLE_OVERRIDE_SHA256" \
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
