# IBAMR ccache setup (Codex sandbox-friendly)

This note documents an IBAMR-only ccache workflow that keeps build products and
cache outside the source tree and inside writable temporary paths.

## Paths

- Source: `/Users/boyceg/code/IBAMR`
- Build: `/tmp/ibamr-objs-dbg-ccache`
- ccache dir: `/tmp/ibamr-ccache`
- ccache temp dir: `/tmp/ibamr-ccache-tmp`

## Environment

```bash
export CCACHE_DIR=/tmp/ibamr-ccache
export CCACHE_TEMPDIR=/tmp/ibamr-ccache-tmp
export CCACHE_COMPILERCHECK=content
```

## Initial setup

```bash
rm -rf /tmp/ibamr-objs-dbg-ccache /tmp/ibamr-ccache /tmp/ibamr-ccache-tmp
mkdir -p /tmp/ibamr-objs-dbg-ccache /tmp/ibamr-ccache /tmp/ibamr-ccache-tmp

source /Users/boyceg/code/autoibamr/dbg/configuration/enable.sh

ccache -z
ccache -M 20G

cmake -S /Users/boyceg/code/IBAMR -B /tmp/ibamr-objs-dbg-ccache \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER="$(which mpicc)" \
  -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
  -DCMAKE_Fortran_COMPILER="$(which mpif90)" \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_INSTALL_PREFIX="" \
  -DHDF5_ROOT="$HDF5_DIR" \
  -DHYPRE_ROOT="$PETSC_DIR" \
  -DIBAMR_FORCE_BUNDLED_BOOST=ON \
  -DIBAMR_FORCE_BUNDLED_EIGEN3=ON \
  -DIBAMR_FORCE_BUNDLED_MUPARSER=ON \
  -DLIBMESH_ROOT="$LIBMESH_DIR" \
  -DLIBMESH_METHOD=devel \
  -DNUMDIFF_ROOT="$NUMDIFF_DIR" \
  -DPETSC_ROOT="$PETSC_DIR" \
  -DSAMRAI_ROOT=/Users/boyceg/sfw/samrai/IBSAMRAI2/darwin-clang-dbg \
  -DSILO_ROOT="$SILO_DIR" \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_Fortran_COMPILER_LAUNCHER=ccache
```

## Build

```bash
source /Users/boyceg/code/autoibamr/dbg/configuration/enable.sh
export CCACHE_DIR=/tmp/ibamr-ccache
export CCACHE_TEMPDIR=/tmp/ibamr-ccache-tmp
export CCACHE_COMPILERCHECK=content
cmake --build /tmp/ibamr-objs-dbg-ccache -j8
```

## Check cache stats

```bash
ccache -s
```
