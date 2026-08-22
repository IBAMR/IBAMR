#!/usr/bin/env bash

set -euo pipefail

readonly oracle_sha=5b77344db6746269f8c77695c99e9043907ba74b

if [[ $# -lt 5 || $# -gt 7 ]]; then
    echo "usage: $0 CANDIDATE_WORKTREE CANDIDATE_BUILD CANDIDATE_SHA ORACLE_WORKTREE MATLAB [OUTPUT_PARENT] [N0|N1|N2]" >&2
    exit 2
fi

candidate_worktree=$1
candidate_build=$2
candidate_sha=$3
oracle_worktree=$4
matlab=$5
output_parent=${6:-/private/tmp}
scope=${7:-N1}
tools_dir=$(cd "$(dirname "$0")" && pwd -P)

if [[ $scope != N0 && $scope != N1 && $scope != N2 ]]; then
    echo "scope must be N0, N1, or N2" >&2
    exit 2
fi
if [[ ! $candidate_sha =~ ^[0-9a-f]{40}$ ]]; then
    echo "candidate SHA must be a full lowercase 40-character Git SHA" >&2
    exit 2
fi

resolved_candidate_sha=$(git -C "$candidate_worktree" rev-parse HEAD)
resolved_oracle_sha=$(git -C "$oracle_worktree" rev-parse HEAD)
candidate_dirty=$(git -C "$candidate_worktree" status --porcelain=v1 --untracked-files=all)
oracle_dirty=$(git -C "$oracle_worktree" status --porcelain=v1 --untracked-files=all)

if [[ $resolved_candidate_sha != "$candidate_sha" ]]; then
    echo "candidate worktree is at $resolved_candidate_sha, expected $candidate_sha" >&2
    exit 1
fi
if [[ $resolved_oracle_sha != "$oracle_sha" ]]; then
    echo "oracle worktree is at $resolved_oracle_sha, expected $oracle_sha" >&2
    exit 1
fi
if [[ -n $candidate_dirty && ${CAV_ALLOW_DIRTY_CANDIDATE:-0} != 1 ]]; then
    echo "candidate worktree is dirty" >&2
    exit 1
fi
if [[ -n $oracle_dirty ]]; then
    echo "oracle worktree is dirty" >&2
    exit 1
fi

cache_file="$candidate_build/CMakeCache.txt"
candidate_executable="$candidate_build/tests/IB/implicit_stokes_ib_solver_components_01"
candidate_template="$candidate_worktree/tests/IB/implicit_stokes_ib_solver_components_01.max_levels=2.pressure_cav_K=1.input"
candidate_petsc_options="$candidate_worktree/tests/IB/PetscOptions.implicit_stokes_ib_solver_components_01.pressure_cav.dat"

if [[ ! -f $cache_file || ! -x $candidate_executable || ! -f $candidate_template || ! -f $candidate_petsc_options ]]; then
    echo "candidate build or live pressure-CAV fixture is incomplete" >&2
    exit 1
fi
configured_source=$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$cache_file")
configured_build_type=$(sed -n 's/^CMAKE_BUILD_TYPE:[^=]*=//p' "$cache_file")
configured_cxx_compiler=$(sed -n 's/^CMAKE_CXX_COMPILER:[^=]*=//p' "$cache_file")
configured_fortran_compiler=$(sed -n 's/^CMAKE_Fortran_COMPILER:[^=]*=//p' "$cache_file")
if [[ $configured_source != "$candidate_worktree" ]]; then
    echo "candidate build was configured from $configured_source, not $candidate_worktree" >&2
    exit 1
fi
if [[ ! -x $matlab ]]; then
    echo "MATLAB executable is not executable: $matlab" >&2
    exit 1
fi

if [[ $scope == N0 ]]; then
    stiffnesses=(1)
else
    stiffnesses=(1 1e2 1e4)
fi

if [[ $scope == N2 ]]; then
    base_grid_cells=8
    finest_grid_cells=16
    patch_count=256
    closure_policies=(RELAXED)
    native_ib_collection_tolerance=1.0e-6
else
    base_grid_cells=4
    finest_grid_cells=8
    patch_count=64
    closure_policies=(RELAXED STRICT)
    native_ib_collection_tolerance=1.0e-10
fi

mkdir -p "$output_parent"
scope_lower=$(printf '%s' "$scope" | tr '[:upper:]' '[:lower:]')
run_root=$(mktemp -d "$output_parent/cav-f2c-$scope_lower.XXXXXX")
candidate_root="$run_root/candidate"
oracle_root="$run_root/oracle"
mkdir -p "$candidate_root" "$oracle_root"

for stiffness in "${stiffnesses[@]}"; do
    case $stiffness in
        1)
            spring_stiffness=8.0e0
            fd_tolerance=1.0e-6
            ;;
        1e2)
            spring_stiffness=8.0e2
            fd_tolerance=1.0e-6
            ;;
        1e4)
            spring_stiffness=8.0e4
            fd_tolerance=2.0e-6
            ;;
    esac
    if [[ $scope == N2 && $stiffness == 1e4 ]]; then
        fd_tolerance=3.0e-6
    fi

    for policy in "${closure_policies[@]}"; do
        candidate_output="$candidate_root/K$stiffness/$policy"
        mkdir -p "$candidate_output"
        if [[ $scope == N2 ]]; then
            spring_stiffness=$(awk -v value="$spring_stiffness" 'BEGIN { printf "%.1e", 2.0*value }')
        fi
        sed -e "s/^N                                 = 4$/N                                 = $base_grid_cells/" \
            -e "s/SPRING_STIFFNESS                  = 8.0e0/SPRING_STIFFNESS                  = $spring_stiffness/" \
            -e "s/FD_REL_TOL                        = 1.0e-6/FD_REL_TOL                        = $fd_tolerance/" \
            -e "s/^EXPECTED_PRESSURE_CAV_PATCH_COUNT = 64$/EXPECTED_PRESSURE_CAV_PATCH_COUNT = $patch_count/" \
            -e "s/^IB_COUPLING_RELATIVE_RESIDUAL_TOL = 1.0e-10$/IB_COUPLING_RELATIVE_RESIDUAL_TOL = $native_ib_collection_tolerance/" \
            -e "s/coupling_aware_asm_closure_policy = \"RELAXED\"/coupling_aware_asm_closure_policy = \"$policy\"/" \
            "$candidate_template" > "$candidate_output/input"
        if [[ $scope == N2 ]]; then
            awk '
                /coarse_solver_db \{/ { in_coarse_solver = 1 }
                !in_coarse_solver && /ksp_type = "preonly"/ { sub(/"preonly"/, "\"richardson\"") }
                in_coarse_solver && /pc_type = "lu"/ { sub(/"lu"/, "\"svd\"") }
                { print }
                /num_post_sweeps = 1/ {
                    print "   U_prolongation_method = \"RT0_REFINE\""
                    print "   P_prolongation_method = \"LINEAR_REFINE\""
                    print "   U_restriction_method = \"RT0_COARSEN\""
                    print "   P_restriction_method = \"CONSERVATIVE_COARSEN\""
                }
            ' "$candidate_output/input" > "$candidate_output/input.with-explicit-transfers"
            mv "$candidate_output/input.with-explicit-transfers" "$candidate_output/input"
        fi
        cp "$candidate_petsc_options" "$candidate_output/$(basename "$candidate_petsc_options")"
        if [[ $scope == N2 ]]; then
            sed -i.bak \
                's/-stokes_ib_pc_level_ksp_type preonly/-stokes_ib_pc_level_ksp_type richardson/' \
                "$candidate_output/$(basename "$candidate_petsc_options")"
            rm "$candidate_output/$(basename "$candidate_petsc_options").bak"
        fi
        (
            cd "$candidate_output"
            "$candidate_executable" input > candidate.stdout 2> candidate.stderr
        )
        if ! rg -q '^test_failures = 0$' "$candidate_output/candidate.stdout"; then
            echo "candidate $policy K=$stiffness live diagnostic failed" >&2
            exit 1
        fi
    done

    oracle_output="$oracle_root/K$stiffness"
    mkdir -p "$oracle_output"
    env CAV_ORACLE_ROOT="$oracle_worktree" \
        CAV_ORACLE_OUTPUT="$oracle_output" \
        CAV_PHYSICAL_K="$stiffness" \
        CAV_ORACLE_FINE_N="$finest_grid_cells" \
        "$matlab" -batch "addpath('$tools_dir'); export_oracle_multiplicative" \
        > "$run_root/oracle-K$stiffness.stdout" 2> "$run_root/oracle-K$stiffness.stderr"
done

replay_status=0
if env CAV_CANDIDATE_ROOT="$candidate_root" \
    CAV_ORACLE_ROOT="$oracle_root" \
    CAV_ORACLE_WORKTREE="$oracle_worktree" \
    CAV_REPLAY_SCOPE="$scope" \
    CAV_REPLAY_REPORT="$run_root/multiplicative-paired-replay-report.json" \
    "$matlab" -batch "addpath('$tools_dir'); replay_multiplicative_common_arithmetic" \
    > "$run_root/replay.stdout" 2> "$run_root/replay.stderr"; then
    replay_status=0
else
    replay_status=$?
fi

{
    echo "schema cav-f2c-multiplicative-provenance-v2"
    echo "scope $scope"
    echo "candidate_sha $resolved_candidate_sha"
    echo "candidate_dirty $([[ -n $candidate_dirty ]] && echo 1 || echo 0)"
    echo "candidate_worktree $candidate_worktree"
    echo "candidate_build $candidate_build"
    echo "candidate_build_type $configured_build_type"
    echo "candidate_cxx_compiler $configured_cxx_compiler"
    echo "candidate_fortran_compiler $configured_fortran_compiler"
    echo "candidate_executable $candidate_executable"
    echo "candidate_template $candidate_template"
    echo "oracle_sha $resolved_oracle_sha"
    echo "oracle_dirty 0"
    echo "oracle_worktree $oracle_worktree"
    echo "matlab $matlab"
    echo "stiffnesses ${stiffnesses[*]}"
    echo "base_grid_cells_per_axis $base_grid_cells"
    echo "finest_grid_cells_per_axis $finest_grid_cells"
    echo "closure_policies ${closure_policies[*]}"
    echo "native_ib_artifact_collection_tolerance $native_ib_collection_tolerance"
    echo "replay_ib_acceptance_tolerance 1.0e-10"
    echo "n2_K1e4_fd_artifact_collection_tolerance 3.0e-6"
    echo "n2_velocity_prolongation_method RT0_REFINE"
    echo "n2_velocity_restriction_method RT0_COARSEN"
    echo "n2_pressure_prolongation_method LINEAR_REFINE"
    echo "n2_pressure_restriction_method CONSERVATIVE_COARSEN"
    echo "n2_pressure_nullspace_coarse_solver svd"
    echo "n2_level_smoother_ksp_type richardson"
} > "$run_root/provenance.txt"

find "$candidate_root" "$oracle_root" -type f -print0 | sort -z | xargs -0 shasum -a 256 > "$run_root/RAW_ARTIFACTS.sha256"
if [[ -f $run_root/multiplicative-paired-replay-report.json ]]; then
    shasum -a 256 "$run_root/multiplicative-paired-replay-report.json" "$run_root/provenance.txt" > "$run_root/REPORTS.sha256"
else
    shasum -a 256 "$run_root/provenance.txt" > "$run_root/REPORTS.sha256"
fi

echo "cav_replay_run_root=$run_root"
echo "cav_replay_report=$run_root/multiplicative-paired-replay-report.json"
if [[ $replay_status -eq 0 ]]; then
    echo "cav_replay_pass=true"
else
    echo "cav_replay_pass=false"
fi
exit "$replay_status"
