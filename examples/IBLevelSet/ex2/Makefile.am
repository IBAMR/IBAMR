## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2022 by the IBAMR developers
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

## Dimension-dependent testers
EXAMPLE_DRIVER = example.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidDensity.h SetFluidGasSolidViscosity.cpp SetFluidGasSolidViscosity.h LSLocateGasInterface.cpp LSLocateGasInterface.h  TagLSRefinementCells.cpp TagLSRefinementCells.h SetLSProperties.cpp SetLSProperties.h GravityForcing.cpp GravityForcing.h FlowGravityForcing.cpp FlowGravityForcing.h LevelSetInitialCondition.h LevelSetInitialCondition.cpp
EXTRA_DIST     = input2d 

EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)


examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[2]d* $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input[2]d ; \
	fi ;


