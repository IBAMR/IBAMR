#!/bin/bash
# Run only in the separately approved disposable Actions checkout/container.
set -euo pipefail
task_source=$(realpath "$1")
task_results=$(realpath -m "$2")
task_package=$(realpath "$3")
task_head=e70a6b6fdb4af99e1cf569daa5279f16e2c25cf3
task_merge=c8c4ff66730c3ce8daec8eac59b600c7cefa1faa
task_target=tests-IB_implicit_stokes_ib_solver_components_01
task_selector='^IB/implicit_stokes_ib_solver_components_01\.(operators\.(backward_euler|midpoint|trapezoidal)|max_levels=2\.galerkin_borrowing)\.input$'
mkdir "$task_results"
exec > >(tee "$task_results/driver.log") 2>&1
cd "$task_package"
sha256sum -c SHA256SUMS
task_statuses="$task_results/status.tsv"
record() {
    local label=$1
    shift
    set +e
    "$@" > "$task_results/$label.log" 2>&1
    local status=$?
    set -e
    printf '%s\t%s\n' "$label" "$status" >> "$task_statuses"
    cat "$task_results/$label.log"
    return "$status"
}
git config --global --add safe.directory "$task_source"
test "$(git -C "$task_source" rev-parse HEAD)" = "$task_head"
test -z "$(git -C "$task_source" status --porcelain)"
record merge-fetch git -C "$task_source" fetch --depth=1 https://github.com/IBAMR/IBAMR.git "$task_merge"
test "$(git -C "$task_source" rev-parse 'HEAD^{tree}')" = "$(git -C "$task_source" rev-parse "$task_merge^{tree}")"
git -C "$task_source" show -s --format=fuller HEAD "$task_merge" > "$task_results/source.txt"
git -C "$task_source" rev-parse 'HEAD^{tree}' "$task_merge^{tree}" >> "$task_results/source.txt"
printf '%s\n' 'wellsd2/ibamr@sha256:78ae125b6152a3f1108dd60c9620ad6561bc7c3d826f9f8ba019aca0104e8642' > "$task_results/image.txt"
record compiler c++ --version
record fortran gfortran --version
record mpi mpiexec --version
record numdiff /numdiff/bin/numdiff --version
cp /petsc/include/petscversion.h /petsc/include/petscconf.h "$task_results/"
cp /petsc/lib/petsc/conf/petscvariables "$task_results/"
cp "$task_package/sidecar.patch" "$task_results/"
sha256sum "$task_package/sidecar.patch" "$task_source/tests/IB/implicit_stokes_ib_solver_components_01.cpp" > "$task_results/input-sha256.txt"
task_build="$task_source/build"
mkdir "$task_build"
cd "$task_build"
record configure cmake -DSILO_ROOT=/petsc -DHYPRE_ROOT=/petsc -DPETSC_ROOT=/petsc -DSAMRAI_ROOT=/samrai \
    -DNUMDIFF_ROOT=/numdiff/ -GNinja -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_COMPILER_LAUNCHER="$(command -v ccache)" \
    -DCMAKE_CXX_FLAGS='-O1 -Wall -Wextra -Wpedantic -Werror -Wno-deprecated-declarations -D_GLIBCXX_ASSERTIONS -ftrivial-auto-var-init=pattern -fuse-ld=mold -D_FORTIFY_SOURCE=3' \
    -DCMAKE_Fortran_FLAGS='-O3 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -Wno-compare-reals' "$task_source"
cp CMakeCache.txt compile_commands.json "$task_results/"
if [[ -f CMakeFiles/CMakeConfigureLog.yaml ]]; then cp CMakeFiles/CMakeConfigureLog.yaml "$task_results/"; fi
record build-original cmake --build . --target "$task_target" -j4
bash "$task_source/tests/link-test-files.sh" "$task_source/tests/IB" "$task_build/tests/IB"
task_executable="$task_build/tests/IB/implicit_stokes_ib_solver_components_01"
capture_output() {
    test ! -e "$task_executable.real"
    mv "$task_executable" "$task_executable.real"
    cp "$task_package/capture-launcher.sh" "$task_executable"
    chmod a+rx "$task_executable"
}
capture_output
# attest writes both successful and failed case work directories below TMPDIR.
chmod a+rx "$task_results"
run_attest() {
    local stage=$1
    mkdir -p "$task_results/$stage/work"
    chown -R build "$task_results/$stage"
    record "$stage" runuser -u build -- env TMPDIR="$task_results/$stage/work" \
        "$task_source/attest" --verbose --keep-work-directories --test-timeout=120 -j1 -R "$task_selector"
}
record discovery runuser -u build -- "$task_source/attest" -N -R "$task_selector"
task_original_status=0
run_attest original || task_original_status=$?
if ! grep -q 'DIFF FAILED' "$task_results/original.log"; then
    printf '%s\n' 'Original numerical DIFF not reproduced; no diagnostic patch applied.' > "$task_results/disposition.txt"
    exit "$task_original_status"
fi
record patch-check git -C "$task_source" apply --check "$task_package/sidecar.patch"
git -C "$task_source" apply "$task_package/sidecar.patch"
git -C "$task_source" diff > "$task_results/applied.patch"
# Restore the original binary before Ninja rebuilds it from the diagnostic TU.
mv "$task_executable.real" "$task_executable"
record build-sidecar cmake --build . --target "$task_target" -j4
capture_output
task_diagnostic_status=0
run_attest sidecar-1 || task_diagnostic_status=$?
run_attest sidecar-2 || task_diagnostic_status=$?
git -C "$task_source" status --porcelain > "$task_results/ephemeral-source-status.txt"
sha256sum "$task_source/tests/IB/implicit_stokes_ib_solver_components_01.cpp" > "$task_results/patched-sha256.txt"
# Leave the job failed when original comparisons fail; do not mask them with successful diagnostics.
if ((task_original_status != 0)); then exit "$task_original_status"; fi
exit "$task_diagnostic_status"
