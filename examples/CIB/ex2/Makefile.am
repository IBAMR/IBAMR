## ---------------------------------------------------------------------
##
## Copyright (c) 2015 - 2019 by the IBAMR developers
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
EXTRA_DIST  = input2d.plate input2d.plate.amr plate2d.vertex plate2d.m petsc_options.dat
EXTRA_DIST +=

EXTRA_PROGRAMS =
EXAMPLES = 
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)


examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d* $(PWD) ; \
	  cp -f $(srcdir)/*.vertex $(PWD) ; \
	  cp -f $(srcdir)/*.m $(PWD) ; \
	  cp -f $(srcdir)/*.dat $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d.plate ; \
	  rm -f $(builddir)/input2d.plate.amr ; \
          rm -f $(builddir)/plate2d.vertex ; \
          rm -f $(builddir)/plate2d.m ; \
          rm -f $(builddir)/petsc_options.dat ; \
	fi ;
