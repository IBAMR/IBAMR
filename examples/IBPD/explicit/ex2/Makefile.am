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
SOURCES = main.cpp generate_failed_sheet2d.cpp generate_failed_sheet3d.cpp
EXTRA_DIST  = input2d.failure input3d.failure
EXTRA_DIST += generate_failed_sheet2d.cpp generate_failed_sheet3d.cpp 

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d failed_sheet2d_generator
endif
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d failed_sheet3d_generator
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = main.cpp

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = main.cpp

failed_sheet2d_generator_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
failed_sheet2d_generator_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
failed_sheet2d_generator_SOURCES = generate_failed_sheet2d.cpp

failed_sheet3d_generator_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
failed_sheet3d_generator_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
failed_sheet3d_generator_SOURCES = generate_failed_sheet3d.cpp

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.failure $(PWD) ; \
	  cp -f $(srcdir)/generate_failed_sheet2d.cpp $(PWD) ; \
	  cp -f $(srcdir)/input3d.failure $(PWD) ; \
	  cp -f $(srcdir)/generate_failed_sheet3d.cpp $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d.failure ; \
	  rm -f $(builddir)/generate_failed_sheet2d.cpp ; \
	  rm -f $(builddir)/input3d.failure ; \
	  rm -f $(builddir)/generate_failed_sheet3d.cpp ; \
	fi ;
