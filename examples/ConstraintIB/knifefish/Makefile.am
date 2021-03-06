## ---------------------------------------------------------------------
##
## Copyright (c) 2016 - 2019 by the IBAMR developers
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
EXAMPLE_DRIVER = KnifeFishKinematics.cpp example.cpp
EXTRA_DIST  = input3d
EXTRA_DIST += KnifeFish.vertex RadiusAmp.dat

EXTRA_PROGRAMS =
EXAMPLES =
if SAMRAI3D_ENABLED
EXAMPLES += main3d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)


examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/KnifeFish.vertex $(PWD) ; \
	  cp -f $(srcdir)/RadiusAmp.dat $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/KnifeFish.vertex ; \
	  rm -f $(builddir)/RadiusAmp.dat ; \
	fi ;
