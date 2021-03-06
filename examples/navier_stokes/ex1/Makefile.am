## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2019 by the IBAMR developers
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
EXAMPLE_DRIVER = example.cpp
EXTRA_DIST     = input2d
EXTRA_DIST    += input3d

EXAMPLES =
CONVERGENCE_TESTERS =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
CONVERGENCE_TESTERS += convergence_tester2d
EXTRA_PROGRAMS += $(EXAMPLES) $(CONVERGENCE_TESTERS)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d
CONVERGENCE_TESTERS += convergence_tester3d
EXTRA_PROGRAMS += $(EXAMPLES) $(CONVERGENCE_TESTERS)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)


convergence_tester2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
convergence_tester2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
convergence_tester2d_SOURCES = convergence_tester.cpp

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)


convergence_tester3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
convergence_tester3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
convergence_tester3d_SOURCES = convergence_tester.cpp

examples: $(EXAMPLES) $(CONVERGENCE_TESTERS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input[23]d ; \
	fi ;
