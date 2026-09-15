#!/usr/bin/env bash
set -euo pipefail
experiment_source=$(cd "$(dirname "$0")/../.." && pwd)
experiment_mode=${1:-Debug}
case "$experiment_mode" in
  Debug) package_mode=dbg; mesh_method=devel; cxx_opt=-O1; fortran_opt=-O2; cpu_flags=; build_name=Debug ;;
  Release) package_mode=opt; mesh_method=opt; cxx_opt=-O3; fortran_opt=-O3; cpu_flags=-mcpu=native; build_name=Release-native ;;
  *) echo "Usage: $0 [Debug|Release]" >&2; exit 1 ;;
esac
experiment_packages=/Users/boyceg/code/autoibamr/$package_mode/packages
if [ "$experiment_mode" = Release ]; then
  source /Users/boyceg/code/autoibamr/opt/configuration/enable.sh
  experiment_samrai=${EXPERIMENT_SAMRAI_ROOT:-$SAMRAI_DIR}
  experiment_boost=$BOOST_DIR
else
  experiment_samrai=/Users/boyceg/sfw/samrai/IBSAMRAI2/darwin-clang-dbg
  experiment_boost=/opt/homebrew
fi
export CCACHE_DIR=${CCACHE_DIR:-/Users/boyceg/Library/Caches/ccache}
export CCACHE_TEMPDIR="$experiment_source/.cache/ibkernels-ccache-tmp"
export CCACHE_BASEDIR="$experiment_source"
mkdir -p "$CCACHE_TEMPDIR"
cmake -S "$experiment_source" -B "$experiment_source/build/ibkernels-matrix-free/$build_name" \
  -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE="$experiment_mode" \
  -DCMAKE_C_COMPILER="$(xcrun --find clang)" -DCMAKE_CXX_COMPILER="$(xcrun --find clang++)" \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_OSX_SYSROOT="$(xcrun --show-sdk-path)" \
  -DCMAKE_C_COMPILER_LAUNCHER=/opt/homebrew/bin/ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=/opt/homebrew/bin/ccache \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_FLAGS="$cxx_opt $cpu_flags -fno-fast-math" \
  -DCMAKE_CXX_FLAGS="$cxx_opt $cpu_flags -fno-fast-math -Wall -Wextra -Wpedantic -Werror -Wmost -Wmove -Wunused -isystem $experiment_packages/libmesh-1.7.8/include" \
  -DCMAKE_Fortran_FLAGS="$fortran_opt $cpu_flags -fno-fast-math -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -Wno-compare-reals" \
  -DIBAMR_ENABLE_TESTING=ON -DIBAMR_ENABLE_DOCUMENTATION=OFF \
  -DBOOST_ROOT="$experiment_boost" \
  -DSAMRAI_ROOT="$experiment_samrai" \
  -DPETSC_ROOT="$experiment_packages/petsc-3.23.3" \
  -DHYPRE_ROOT="$experiment_packages/petsc-3.23.3" \
  -DHDF5_ROOT="$experiment_packages/hdf5-1.12.2" \
  -DLIBMESH_ROOT="$experiment_packages/libmesh-1.7.8" -DLIBMESH_METHOD="$mesh_method" \
  -DSILO_ROOT="$experiment_packages/silo-4.11-bsd" \
  -DNUMDIFF_ROOT="$experiment_packages/numdiff-5.9.0"
