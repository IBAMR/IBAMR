# Working with IBAMR's Continuous Integration (CI) setup

## Introduction

Changes proposed to IBAMR need to not introduce new bugs or compilation
problems. To this end, IBAMR presently uses GitHub actions to, for each pull
request, compile the library, all tests and examples, and run a subset of the
test suite. For historical reasons we use our own test runner, attest, which
understands IBAMR-specific things (like running with MPI and restart files).

## How it works

In order to compile IBAMR we first need a basic computing environment (MPI,
compilers, a shell, etc.) as well as IBAMR's dependencies. To this end, all of
IBAMR's build dependencies are packaged inside one of several docker containers.
These containers provide a portable environment in which we can compile the
library with specified precompiled versions of PETSc, SAMRAI, etc. We maintain
three different CI configurations:
- A newer configuration (based on Arch Linux) with reasonably recent versions
  of dependencies
- The oldest supported version of Ubuntu
- macOS

The file specifying what GitHub does is `.github/workflows/push-pull.yml`. This
YAML file specifies what each virtual machine (hosted by GitHub) should do when
a new pull request is created. This file is responsible for setting up the
environment for compiling the library: in particular, it pulls down David Wells'
Docker images from DockerHub to get copies of all recent dependencies.

Scripts for compiling docker images are available at

https://github.com/IBAMR/docker-files

## Improving performance with ccache

IBAMR takes upwards of an hour to compile on one of GitHub's virtual cloud
machines. Most of this time is spent compiling the same files between runs. To
accelerate this process, we use GitHub's cache action (`actions/cache`) which
stores build artifacts between runs. To speed up compilation, ccache will run
the C preprocessor on each file and hash the result: if the hash matches a cache
entry, then ccache returns that value (an object file) instead of running the
compiler itself. This lowers compilation time by about 80%.

GitHub's cache infrastructure was not intended to do this. To make things work,
we rebuild the cache on the master branch after merges and pull that into
feature branches. Feature branches cannot access caches from other feature
branches.

## Actually compiling the library and running the test suite

For more information see the actual Dockerfiles: we use them, plus the copies of
ccache and ninja in the container, to run attest.
