## ---------------------------------------------------------------------
##
## Copyright (c) 2025 - 2025 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = array_01 array_02

array_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
array_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
array_01_SOURCES = array_01.cpp

array_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
array_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
array_02_SOURCES = array_02.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
