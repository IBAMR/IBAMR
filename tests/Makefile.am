## ---------------------------------------------------------------------
##
## Copyright (c) 2015 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

include $(top_srcdir)/config/Make-rules
TEST_DIRECTORIES = adv_diff advect CIB ConstraintIB coarsen complex_fluids \
fe_mechanics IB IBFE IIM IBTK interpolate level_set multiphase_flow \
navier_stokes physical_boundary refine spread vc_navier_stokes wave_tank

tests: $(TEST_DIRECTORIES)

$(TEST_DIRECTORIES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) tests
.PHONY: tests $(TEST_DIRECTORIES)
