## ---------------------------------------------------------------------
##
## Copyright (c) 2016 - 2020 by the IBAMR developers
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
SOURCES = main.cpp generate_membrane2d.cpp
EXTRA_DIST  = input2d input2d_tune automate.sh run_IPD_comparisons.sh
EXTRA_DIST += read_file_ibamr2d.m

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d membrane2d_generator
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = main.cpp

membrane2d_generator_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
membrane2d_generator_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
membrane2d_generator_SOURCES = generate_membrane2d.cpp

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	  cp -f $(srcdir)/generate_membrane2d.cpp $(PWD) ; \
	  cp -f $(srcdir)/read_file_ibamr2d.m $(PWD) ; \
	  cp -f $(srcdir)/automate.sh $(PWD) ; \
	  cp -f $(srcdir)/run_IPD_comparisons.sh $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/generate_membrane2d.cpp ; \
	  rm -f $(builddir)/read_file_ibamr2d.m ; \
	  rm -f $(builddir)/automate.sh ; \
	  rm -f $(builddir)/run_IPD_comparisons.sh ; \
	fi ;