#!/usr/bin/env bash

set -euo pipefail

readonly oracle_sha=5b77344db6746269f8c77695c99e9043907ba74b

if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "usage: $0 CANDIDATE_WORKTREE CANDIDATE_BUILD CANDIDATE_SHA ORACLE_WORKTREE MATLAB [OUTPUT_PARENT]" >&2
    exit 2
fi

candidate_worktree=$1
candidate_build=$2
candidate_sha=$3
oracle_worktree=$4
matlab=$5
output_parent=${6:-/private/tmp}
tools_dir=$(cd "$(dirname "$0")" && pwd -P)

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
candidate_input="$candidate_worktree/tests/IB/implicit_stokes_ib_solver_components_01.max_levels=2.pressure_cav_K=1.input"
candidate_petsc_options="$candidate_worktree/tests/IB/PetscOptions.implicit_stokes_ib_solver_components_01.pressure_cav.dat"

if [[ ! -f $cache_file || ! -x $candidate_executable || ! -f $candidate_input || ! -f $candidate_petsc_options ]]; then
    echo "candidate build or N0 fixture is incomplete" >&2
    exit 1
fi
configured_source=$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$cache_file")
if [[ $configured_source != "$candidate_worktree" ]]; then
    echo "candidate build was configured from $configured_source, not $candidate_worktree" >&2
    exit 1
fi
if [[ ! -x $matlab ]]; then
    echo "MATLAB executable is not executable: $matlab" >&2
    exit 1
fi

mkdir -p "$output_parent"
run_root=$(mktemp -d "$output_parent/cav-f2c-n0.XXXXXX")
candidate_output="$run_root/candidate"
oracle_output="$run_root/oracle"
mkdir -p "$candidate_output" "$oracle_output"

cp "$candidate_input" "$candidate_output/input"
cp "$candidate_petsc_options" "$candidate_output/$(basename "$candidate_petsc_options")"

(
    cd "$candidate_output"
    "$candidate_executable" input > candidate.stdout 2> candidate.stderr
)

env CAV_N0_ORACLE_ROOT="$oracle_worktree" \
    CAV_N0_ORACLE_OUTPUT="$oracle_output" \
    "$matlab" -batch "addpath('$tools_dir'); export_oracle_n0" \
    > "$run_root/oracle.stdout" 2> "$run_root/oracle.stderr"

env CAV_N0_CANDIDATE_OUTPUT="$candidate_output" \
    CAV_N0_ORACLE_OUTPUT="$oracle_output" \
    CAV_N0_REPORT="$run_root/n0-paired-replay-report.json" \
    "$matlab" -batch "addpath('$tools_dir'); replay_n0_common_arithmetic" \
    > "$run_root/replay.stdout" 2> "$run_root/replay.stderr"

{
    echo "schema cav-f2c-n0-provenance-v1"
    echo "candidate_sha $resolved_candidate_sha"
    echo "candidate_dirty $([[ -n $candidate_dirty ]] && echo 1 || echo 0)"
    echo "candidate_worktree $candidate_worktree"
    echo "candidate_build $candidate_build"
    echo "candidate_executable $candidate_executable"
    echo "oracle_sha $resolved_oracle_sha"
    echo "oracle_dirty 0"
    echo "oracle_worktree $oracle_worktree"
    echo "matlab $matlab"
} > "$run_root/provenance.txt"

find "$candidate_output" "$oracle_output" -type f -print0 | sort -z | xargs -0 shasum -a 256 > "$run_root/RAW_ARTIFACTS.sha256"
shasum -a 256 "$run_root/n0-paired-replay-report.json" "$run_root/provenance.txt" > "$run_root/REPORTS.sha256"

echo "cav_n0_run_root=$run_root"
echo "cav_n0_report=$run_root/n0-paired-replay-report.json"
echo "cav_n0_pass=true"
