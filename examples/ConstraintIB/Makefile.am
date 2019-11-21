## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2018 by the IBAMR developers
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
SUBDIRS = eel2d eel3d

COMPLETE_EXAMPLES = eel2d eel3d falling_sphere flow_past_cylinder flow_past_cylinder_HF knifefish impulsively_started_cylinder moving_plate oscillating_rigid_cylinder stokes_first_problem
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

## Make targets for gtests
gtest:
	@(cd eel2d                        && make gtest) || exit 1;
	@(cd eel3d                        && make gtest) || exit 1;
	@(cd falling_sphere               && make gtest) || exit 1;
	@(cd flow_past_cylinder           && make gtest) || exit 1;
	@(cd flow_past_cylinder_HF           && make gtest) || exit 1;
	@(cd knifefish                    && make gtest) || exit 1;
	@(cd impulsively_started_cylinder && make gtest) || exit 1;
	@(cd moving_plate                 && make gtest) || exit 1;
	@(cd oscillating_rigid_cylinder   && make gtest) || exit 1;
	@(cd stokes_first_problem         && make gtest) || exit 1;

gtest-long:
	@(cd eel2d                        && make gtest) || exit 1;
	@(cd eel3d                        && make gtest) || exit 1;
	@(cd falling_sphere               && make gtest) || exit 1;
	@(cd flow_past_cylinder           && make gtest) || exit 1;
	@(cd flow_past_cylinder_HF           && make gtest) || exit 1;
	@(cd knifefish                    && make gtest) || exit 1;
	@(cd impulsively_started_cylinder && make gtest) || exit 1;
	@(cd moving_plate                 && make gtest) || exit 1;
	@(cd oscillating_rigid_cylinder   && make gtest-long) || exit 1;
	@(cd stokes_first_problem         && make gtest) || exit 1;

gtest-all:
	@(cd eel2d                        && make gtest-long) || exit 1;
	@(cd eel3d                        && make gtest-long) || exit 1;
	@(cd falling_sphere               && make gtest-long) || exit 1;
	@(cd flow_past_cylinder           && make gtest-long) || exit 1;
	@(cd flow_past_cylinder_HF           && make gtest-long) || exit 1;
	@(cd knifefish                    && make gtest-long) || exit 1;
	@(cd impulsively_started_cylinder && make gtest-long) || exit 1;
	@(cd moving_plate                 && make gtest-long) || exit 1;
	@(cd oscillating_rigid_cylinder   && make gtest-long) || exit 1;
	@(cd stokes_first_problem         && make gtest-long) || exit 1;
