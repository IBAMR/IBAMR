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

## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = example.cpp SetFluidProperties.cpp SetLSProperties.cpp LSLocateColumnInterface.cpp GravityForcing.cpp TagLSRefinementCells.cpp
EXTRA_DIST     = input2d
EXTRA_DIST    += input2d_CaseA input2d_CaseA_IrregularWave input2d_CaseB input2d_CaseC input2d_5thOrderStokes
EXTRA_DIST    += input3d
EXTRA_DIST    += compareIrregularWave.m compareStokesFirstOrder.m compareStokesSecondOrder.m compareStokesFifthOrder.m

EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)



examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d* $(PWD) ; \
	  cp -f $(srcdir)/*.m $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input[23]d               ; \
	  rm -f $(builddir)/input[2]d_CaseA          ; \
	  rm -f $(builddir)/input[2]d_CaseB          ; \
	  rm -f $(builddir)/input[2]d_CaseC          ; \
	  rm -f $(builddir)/input[2]d_5thOrderStokes ; \
	fi ;
