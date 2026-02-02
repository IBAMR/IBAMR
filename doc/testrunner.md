# Using `attest`, IBAMR's test suite runner

## Introduction

IBAMR includes a test suite consisting of both unit tests (tests of
individual functionality implemented in IBAMR) and integration tests
(tests of codes written using IBAMR). These tests range from simply
printing some output on an 8x8 grid to running a three-dimensional
problem with MPI. The goal is to ultimately implement unit or
integration tests for all features in IBAMR.

## Design Goals

The goal of the test suite is to design a system that lets IBAMR
developers easily write unit tests for small pieces of functionality.
The tests should be compiled by the build system like other example
programs. The test runner should be completely separate from the build
system and should be able to run tests in parallel and, when requested,
use multiple MPI processes for each test.

## Requirements

The test runner, `attest`, requires Python 3.5 or newer and numdiff. The
location of the numdiff installation is provided to cmake via `-DNUMDIFF_ROOT`.
CMake extracts the paths to the numdiff and mpi executables as well as some other
basic parameters to populate `attest.conf`, which usually depends on the current
build environment.

Tests are compiled via the `tests` target (e.g., via `make -j4 tests`).

## Files

The test suite consists of three components:

1.  Test executables: these are executables compiled and linked against IBAMR by
    the build system. They are implemented in nearly the same way as the example
    programs. All source files for the test executables are located in
    subdirectories of `tests/` in the top directory.
2.  The test input and output files: IBAMR's test suite assumes that each test
    corresponds to exactly one pair of input and output files. The file names
    for the input and output file starts with the name of the required test
    executable followed by a period, with additional data separated by periods.
    For example: the test executable `interpolate_velocity_01` is used in
    multiple tests since there are several input and output files starting with
    `interpolate_velocity_01` (though each has different text after the period).
    The number of MPI processes is encoded into the test input and output files
    by writing, e.g., ``.mpirun=42.` in the filename. The build system will
    create symbolic links to these files inside the build directory.
3.  The test runner: this is the python script `attest` in the top-level
    directory. The build system creates a symbolic link to this script and
    generates a configuration file, `attest.conf`, inside the build directory.

## attest

The program that finds and executes tests is the test runner: `attest`. `attest`
is very similar, aside from being decoupled from the build system, to `ctest`,
and supports many of the same options. All options can be displayed by running
`./attest --help`: a few common options include

1.  `-j,--jobs N Running attest -j4` or `attest --jobs 4` will run four job
    processes at once. This could be one test with 4 MPI processes, two tests
    each using two MPI processes, etc.
2.  `-R,--tests-regex TEST_REGEX` Running `attest -R interp` will filter the
    tests and only run the tests with `interp` in their names.
3.  `-N,--show-only` Running `attest -N` will disable the actual execution of
    tests, but will print the names of all tests that would have been run.

Since the build system compiles the executables and links input and output files
into the build directory, `attest` simply looks for tests in the `tests/`
directory in which it is run. Hence, running

      make -j8 tests
      ./attest -j8

from the top build directory will compile and run the whole test suite while
using eight processors.
