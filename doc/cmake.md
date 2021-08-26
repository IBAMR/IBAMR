# The CMake-based IBAMR build system

## Introduction

IBAMR 0.8.0 introduced a complete rewrite of its build system (that is, the set
of scripts that detect IBAMR's dependencies and then compile and install the
library). The new build system requires CMake version 3.15 or newer.

CMake is a build system generator - i.e., unlike traditional autotools usage,
users will use the `cmake` executable to generate a build system instead of
relying on a prepackaged `configure` script. Hence users must have CMake
installed on their computer to compile IBAMR.

## Advantages of the CMake build system

A few things are now possible with CMake that were not possible before:
1. `make install` works correctly.
2. In the future, IBAMR will go one step further than installation and offer
   downloadable Mac packages.
3. IBAMR can now be compiled with either dynamicly or staticly linked libraries.
4. IBAMR now exports its targets - i.e., IBAMR can be added as a dependency to
   your own projects with the standard `FIND_PACKAGE` CMake macro, which sets up
   the complete header and link interface.
5. The build system itself is much easier to maintain and modify since it is written
   in CMake code instead of shell and m4. Boilerplate files like `aclocal.m4` no
   longer exist and there is no need to rerun `autoreconf`, `autoheader`, or any
   other such program (CMake will do this automatically if scripts are updated).
6. The CMake build system integrates well with most IDEs.

## How to use the CMake build system

By default, IBAMR will compile itself with shared libraries. *This requires that
SAMRAI be compiled with position-independent code: i.e., SAMRAI must be compiled
with `-fPIC`*. If you do not compile SAMRAI this way then you must use static
libraries (see the text below on static linkage).

This is similar to how autotools is used with a separate build directory
`build`. IBAMR's configuration scripts expect all dependencies to have their
paths provided in the standard way for CMake (i.e., we pass in the root path of
the directory for things installed outside the system search path). Here, for
example, `PETSC_ROOT` is typically `$PETSC_DIR/$PETSC_ARCH`. IBAMR's configuration
script looks for each package in `PKG_ROOT`, where `PKG` is the all-caps version
of the name: for example, CMake will look for boost in the directory specified
by `BOOST_ROOT` and for muParser in the directory specified by `MUPARSER_ROOT`.

Same packages have ambiguous capitalizations - IBAMR expects package names to be
in `ALL_CAPS` (e.g., `PETSC_ROOT` and not `Petsc_ROOT`). Some aliases are also
defined, mostly so that one can use the common CMake package names during
configuration:
- `BOOST_ROOT` and `Boost_ROOT` are equivalent
- `EIGEN3_ROOT` and `Eigen3_ROOT` are equivalent
- `LIBMESH_ROOT` and `libMesh_ROOT` are equivalent
- `MUPARSER_ROOT` and `muParser_ROOT` are equivalent
- `PETSC_ROOT` and `PETSc_ROOT` are equivalent

though inside IBAMR's build system the `ALL_CAPS` names are used.

```
mkdir build
cd build
cmake -DCMAKE_C_FLAGS="-O3 -march=native"                     \
      -DCMAKE_CXX_FLAGS="-O3 -march=native"                   \
      -DCMAKE_Fortran_FLAGS="-O3 -march=native"               \
      -DCMAKE_INSTALL_PREFIX=$HOME/Applications/ibamr         \
      -DIBAMR_ENABLE_TESTING=OFF                              \
      -DSAMRAI_ROOT=$HOME/Applications/samrai-2.4.4           \
      -DLIBMESH_ROOT=$HOME/Applications/libmesh-dev           \
      -DLIBMESH_METHOD=OPT                                    \
      -DPETSC_ROOT=$HOME/Applications/petsc-3.13.0/x86_64-opt \
      -DMUPARSER_ROOT=$HOME/Applications/muParser-2.3.2/      \
      ../
make -j6
make -j6 install
```
IBAMR requires
- Boost,
- Eigen3,
- HDF5,
- Hypre,
- MPI,
- muParser,
- PETSc, and
- SAMRAI.

Of these, Boost, Eigen3, and muParser are also bundled with IBAMR and IBAMR will
build its own copies of those libraries if it cannot find working externally
available versions. If you wish to use the bundled version of one of these
libraries instead of one found at either system or specified search locations,
then you can specify that as well by passing any of
- `-DIBAMR_FORCE_BUNDLED_BOOST=ON`
- `-DIBAMR_FORCE_BUNDLED_EIGEN3=ON`
- `-DIBAMR_FORCE_BUNDLED_MUPARSER=ON`
as options to CMake.

The CMake build system will attempt to find these with default search paths and
also in the directory provided by `PKG_ROOT` (i.e., for Boost, the build system
will search in `BOOST_ROOT`).

IBAMR's optional dependencies are
- GSL (which is used in some examples),
- Silo, and
- libMesh.

These are only searched for if a directory path is provided (i.e., by passing
`-DSILO_ROOT=/usr/`). libMesh also requires that the caller pass
`LIBMESH_METHOD` to CMake to indicate which libMesh variant we should use (e.g.,
above we use the optimized version by passing `OPT`: `DBG` and `DEVELOP` are
alternative options).

Subsequent runs of CMake are much faster than the initial run since CMake will
cache some computed values (such as the MPI configuration). If you run into
problems with CMake, a good way to get out of trouble is to delete
`CMakeCache.txt`, which will delete the cache, and the `CMakeFiles` directory,
which contains all of the files generated by CMake. Deleting the cache is
harmless - the subsequent initial configuration run will just take a few more
seconds.

### Configuring the build
- The example call to CMake above disables tests. If you wish to compile the
  tests then pass `-DIBAMR_ENABLE_TESTING=ON` instead.
- If you want to build IBAMR with static libraries then pass the argument
  `-DBUILD_SHARED_LIBS=OFF` to the initial call to `cmake`. IBAMR defaults to
  building shared libraries.
- At the current time the build system does not support compiling with CMake's
  own debug or release modes since this is not compatible with the way things
  were done with the autotools based build system. If you want to compile a
  debug build or an optimized build you will have to provide the flags yourself,
  like what was done at the top of the previous section.

## Debugging configuration issues
- If you need to change a parameter passed to CMake and things don't work then
  try deleting the `CMakeCache.txt` file.
- If CMake cannot find your MPI implementation then try explicitly passing in
  the MPI compiler wrappers - e.g.,
 ```sh
     cmake -DCMAKE_C_COMPILER=/path/to/mpicc \
           -DCMAKE_CXX_COMPILER=/path/to/mpicxx \
           -DCMAKE_Fortran_COMPILER=/path/to/mpifort \
           -DPETSC_ROOT=/path/to/petsc
     # etc.
 ```
- The CMake build system is new and may still contain bugs. If things don't seem
  to work then post a message on the IBAMR mailing list.
- If all else fails use the autotools-based build system, which is currently
  known to work well on a wide variety of computers.

## How to use IBAMR in your own project

If you have installed IBAMR correctly then it is easy to add it as a dependency
to your own project using standard CMake tools. You can configure your project
(here assumed to consist of a single file, `main.cpp`, which will be compiled in
2D) with the following `CMakeLists.txt` script:
```cmake
PROJECT(main2d)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

ADD_EXECUTABLE(main2d main.cpp)

FIND_PACKAGE(IBAMR 0.8.0 REQUIRED)
TARGET_LINK_LIBRARIES(main2d IBAMR::IBAMR2d)
```
Run cmake as
```
cmake -DIBAMR_ROOT=/path/to/ibamr/installation ./
```
where `/path/to/ibamr/installation` should be replaced by the directory path
(either relative or absolute) to the installation location of IBAMR or IBAMR's
build directory.

## Build system internals

### General setup

IBAMR is now compiled as four libraries: 2D and 3D IBTK libraries and 2D and 3D
IBAMR libraries. Since IBAMR depends on IBTK, many features in the library are
actually set up for IBTK and then copied into IBAMR proper. IBAMR and IBTK have
their dependencies configured in the top-level `CMakeLists.txt` file. Source
files for each library are listed in `src/CMakeLists.txt` and
`ibtk/src/CMakeLists.txt`, respectively. Similarly, all test executables are
compiled in `tests/CMakeLists.txt`.

Since there is a wide variety of necessary setup for examples, IBAMR includes
the macro `IBAMR_ADD_EXAMPLE` (defined in `examples/CMakeLists.txt`) which is
used to set up a single example program (which may have multiple source and
input files). New examples should use this macro as well.

### Configuration headers

At the current time IBAMR now generates a single configuration header:
`ibtk/config.h`. All macros defined in this file with the prefix `IBTK_` are
redefined in `ibamr/config.h` with the prefix `IBAMR_`. The old top-level
configuration files `IBAMR_config.h` and `IBTK_config.h` do not have equivalents
in the CMake build system.
