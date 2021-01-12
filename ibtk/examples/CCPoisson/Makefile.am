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
EXTRA_DIST     = README input2d input3d

EXTRA_PROGRAMS =
EXAMPLES =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d
endif
EXTRA_PROGRAMS += $(EXAMPLES)

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBTK_LDFLAGS) $(IBTK2d_LIBS) $(IBTK_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBTK_LDFLAGS) $(IBTK3d_LIBS) $(IBTK_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/{README,input[23]d} $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/{README,input[23]d} ; \
	fi ;
