# Working with IBAMR's Continuous Integration (CI) setup

## Introduction

Changes proposed to IBAMR need to not introduce new bugs or compilation
problems. To this end, IBAMR presently uses GitHub actions to, for each pull
request, compile the library, all tests and examples, and run a subset of the
test suite. Since IBAMR supports both autotools and CMake we use our own test
runner, attest, which has some preliminary integration with ctest and cdash
for monitoring builds.

## How it works

In order to compile IBAMR we first need a basic computing environment (MPI,
compilers, a shell, etc.) as well as IBAMR's dependencies. To this end, all of
IBAMR's build dependencies are packaged inside one of several docker containers.
These containers provide a portable environment in which we can compile the
library with specified precompiled versions of PETSc, SAMRAI, etc. Most of our
containers are based on recent releases of Fedora: one is based on CENTOS 7 to
verify that the library still works with GCC 4.8 and minimal C++11 support.

The file specifying what GitHub does is `.github/workflows/push-pull.yml`. This
YAML file specifies what each virtual machine (hosted by GitHub) should do when
a new pull request is created. This file is responsible for setting up the
environment for compiling the library: in particular, it pulls down David Wells'
Docker images (either Fedora 33 or CENTOS 7) from DockerHub to get copies of all
recent dependencies.

Scripts for compiling docker images are available at

https://github.com/IBAMR/docker-files

## Improving performance with sccache and ccache

The Docker setup we use was initially developed by Kitware and uses scache; for
all intents and purposes, however, sccache and ccache work in the same way so I
will collectively refer to them as ccache. ccache is the more standard tool and
we should ultimately migrate away from sccache to using ccache instead. The
additional features of sccache are not relevant to our workflow on GitHub.

IBAMR takes upwards of an hour to compile on one of GitHub's virtual cloud
machines. Most of this time is spent compiling the same files between runs. To
accelerate this process, we use GitHub's cache action (`actions/cache@v2`) which
stores build artifacts between runs. To speed up compilation, ccache will run
the C preprocessor on each file and hash the result: if the hash matches a cache
entry, then ccache returns that value (an object file) instead of running the
compiler itself. This lowers compilation time by about 80%.

GitHub's cache infrastructure was not intended to do this. To make it work
better with our needs, we assign each cache, at the end of each run, a key equal
to the present GitHub action number and some additional information about the
particular compilation run (e.g., whether or not we used CMake). Subsequent runs
will look for compilation caches with numbers close to their own: this way we
are constantly cycling through cached dependencies and will attempt to use the
newest caches. GitHub will automatically delete caches that have been unused
for seven days, or caches over the 10 GB limit.

## Actually compiling the library and running the test suite

At the present time we use ctest to compile the library and run the test suite.
In particular, the ctest scripts in `.github/cmake/` have enough information to
set up the library with image-specific settings, compile the library, run the
test suite, and report the results. Each setup uses a different configuration
file, which in turn configures the library in a different way.
