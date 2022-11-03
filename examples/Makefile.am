## ---------------------------------------------------------------------
##
## Copyright (c) 2009 - 2022 by the IBAMR developers
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
SUBDIRS = CIB ConstraintIB DLM IB IBFE IBLevelSet IIM IMP adv_diff advect complex_fluids fe_mechanics navier_stokes level_set multiphase_flow vc_navier_stokes wave_tank

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

COMPLETE_EXAMPLES = CIB ConstraintIB DLM IB IBFE IBLevelSet IIM adv_diff advect complex_fluids fe_mechanics navier_stokes level_set multiphase_flow vc_navier_stokes wave_tank
# Not all examples can be compiled at the present time
INCOMPLETE_EXAMPLES = IMP
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

