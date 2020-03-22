#! /bin/bash
## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2019 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------


REBUILD_LIBRARIES=0
PROJECT_ROOT=$PWD

SETUP_PETSC=0
PETSC_DIR=$HOME/sfw/petsc/petsc-maint
PETSC_ARCH=darwin-dbg

SETUP_SAMRAI=0
SAMRAI_BUILD_DIR=$HOME/sfw/samrai/2.4.4/objs-dbg

SETUP_LIBMESH=0
LIBMESH_BUILD_DIR=$HOME/sfw/libmesh/master/build-dbg

SETUP_IBAMR=1
IBAMR_BUILD_DIR=$HOME/code/ibamr-objs-dbg
IBAMR_PROJECT_ROOT=$HOME/code/IBAMR

intercept_make="make CC=intercept-cc CXX=intercept-c++ MPICH_CC=intercept-cc MPICH_CXX=intercept-cxx OMPI_CC=intercept-cc OMPI_CXX=intercept-cxx"

CXX_COMMAND=/usr/local/opt/llvm/bin/clang++

if [ $SETUP_PETSC -eq 1 ]; then
  echo "setting up petsc"
  cd $PETSC_DIR
  echo "working directory is: $PWD"
  echo "PETSC_ARCH is: $PETSC_ARCH"
  make PETSC_ARCH=$PETSC_ARCH clean
  rm $PWD/compile_commands.json
  intercept-build --override-compiler $intercept_make PETSC_ARCH=$PETSC_ARCH all
  compdb -p . list > compile_commands-new.json
  sed -i -e "s|\"command\": \"c++|\"command\": \"$CXX_COMMAND|" compile_commands-new.json
  cp $PWD/compile_commands-new.json $PROJECT_ROOT/petsc-compile_commands.json
  if [ $REBUILD_LIBRARIES -eq 1 ]; then
    make PETSC_ARCH=$PETSC_ARCH all
  fi
  cd $PROJECT_ROOT
fi

if [ $SETUP_SAMRAI -eq 1 ]; then
  echo "setting up samrai"
  cd $SAMRAI_BUILD_DIR
  echo "working directory is: $PWD"
  make clean
  rm $PWD/compile_commands.json
  intercept-build --override-compiler $intercept_make
  compdb -p . list > compile_commands-new.json
  sed -i -e "s|\"command\": \"c++|\"command\": \"$CXX_COMMAND|" compile_commands-new.json
  cp $PWD/compile_commands-new.json $PROJECT_ROOT/samrai-compile_commands.json
  if [ $REBUILD_LIBRARIES -eq 1 ]; then
    make all
  fi
  cd $PROJECT_ROOT
fi

if [ $SETUP_LIBMESH -eq 1 ]; then
  echo "setting up libmesh"
  cd $LIBMESH_BUILD_DIR
  echo "working directory is: $PWD"
  make clean
  rm $PWD/compile_commands.json
  intercept-build --override-compiler $intercept_make
  compdb -p . list > compile_commands-new.json
  sed -i -e "s|\"command\": \"c++|\"command\": \"$CXX_COMMAND|" compile_commands-new.json
  cp $PWD/compile_commands-new.json $PROJECT_ROOT/libmesh-compile_commands.json
  if [ $REBUILD_LIBRARIES -eq 1 ]; then
    make all
  fi
  cd $PROJECT_ROOT
fi

if [ $SETUP_IBAMR -eq 1 ]; then
  echo "setting up ibamr"
  cd $IBAMR_BUILD_DIR
  echo "working directory is: $PWD"
  make clean
  rm $PWD/compile_commands.json
  intercept-build --override-compiler $intercept_make
  cp compile_commands.json compile_commands.json.bak
  compdb -p . list > compile_commands-new.json
  sed -i -e "s|\"command\": \"c++|\"command\": \"$CXX_COMMAND|" compile_commands-new.json
  cp $PWD/compile_commands-new.json $PROJECT_ROOT/ibamr-compile_commands.json
  if [ $REBUILD_LIBRARIES -eq 1 ]; then
    make all
  fi
  cd $PROJECT_ROOT
fi

# Use the IBAMR compile_commands.json for now.
#
# TODO: Should we merge the compile_commands.json files for all libraries together?
ln -s $PWD/ibamr-compile_commands.json $IBAMR_PROJECT_ROOT/compile_commands.json
