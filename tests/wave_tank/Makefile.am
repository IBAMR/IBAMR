## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2020 by the IBAMR developers
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

EXTRA_PROGRAMS = nwt
if LIBMESH_ENABLED
EXTRA_PROGRAMS += nwt_cylinder 
endif

nwt_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 
nwt_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
nwt_SOURCES = nwt.cpp

if LIBMESH_ENABLED
# this test needs some extra input files, so make SOURCE_DIR available:
nwt_cylinder_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
nwt_cylinder_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
nwt_cylinder_SOURCES = nwt_cylinder.cpp
endif

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
