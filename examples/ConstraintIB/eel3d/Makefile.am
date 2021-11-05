## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2021 by the IBAMR developers
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
EXTRA_DIST  = input3d
EXTRA_DIST += eelgenerator3d.cpp


EXTRA_PROGRAMS =
EXAMPLES =
if SAMRAI3D_ENABLED
EXAMPLES += main3d eelgenerator3d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = IBEELKinematics3d.cpp $(EXAMPLE_DRIVER)

eelgenerator3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
eelgenerator3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
eelgenerator3d_SOURCES = eelgenerator3d.cpp


examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/eel3d.vertex $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/eel3d.vertex ;\
	fi ;
