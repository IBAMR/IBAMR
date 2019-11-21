## ---------------------------------------------------------------------
##
## Copyright (c) 2011 - 2017 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = CCLaplace CCPoisson PETScOps PhysBdryOps SCLaplace SCPoisson VCLaplace VCViscousSolver

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

## Standard make targets.
examples:
	@(cd CCLaplace   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd CCPoisson   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd PETScOps    && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd PhysBdryOps && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd SCLaplace   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd SCPoisson   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd VCLaplace   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd VCViscousSolver   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;


gtest:
	@(cd CCLaplace   && make gtest) || exit 1;
	@(cd CCPoisson   && make gtest) || exit 1;
	@(cd PETScOps    && make gtest) || exit 1;
##	@(cd PhysBdryOps && make gtest) || exit 1;
	@(cd SCLaplace   && make gtest) || exit 1;
	@(cd SCPoisson   && make gtest) || exit 1;
	@(cd VCLaplace   && make gtest) || exit 1;
	@(cd VCViscousSolver   && make gtest) || exit 1;

gtest-long:
	@(cd CCLaplace   && make gtest) || exit 1;
	@(cd CCPoisson   && make gtest) || exit 1;
	@(cd PETScOps    && make gtest) || exit 1;
##	@(cd PhysBdryOps && make gtest-long) || exit 1;
	@(cd SCLaplace   && make gtest-long) || exit 1;
	@(cd SCPoisson   && make gtest-long) || exit 1;
	@(cd VCLaplace   && make gtest-long) || exit 1;
	@(cd VCViscousSolver   && make gtest-long) || exit 1;

gtest-all:
	@(cd CCLaplace   && make gtest-long) || exit 1;
	@(cd CCPoisson   && make gtest-long) || exit 1;
	@(cd PETScOps    && make gtest-long) || exit 1;
##	@(cd PhysBdryOps && make gtest-long) || exit 1;
	@(cd SCLaplace   && make gtest-long) || exit 1;
	@(cd SCPoisson   && make gtest-long) || exit 1;
	@(cd VCLaplace   && make gtest-long) || exit 1;
	@(cd VCViscousSolver   && make gtest-long) || exit 1;
