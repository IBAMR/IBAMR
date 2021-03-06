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

EXTRA_PROGRAMS =
if LIBMESH_ENABLED
EXTRA_PROGRAMS += spread_01_2d spread_01_3d spread_02_2d spread_02_3d
endif

if LIBMESH_ENABLED
spread_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
spread_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
spread_01_2d_SOURCES = spread_01.cpp
endif

if LIBMESH_ENABLED
spread_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
spread_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
spread_01_3d_SOURCES = spread_01.cpp
endif

if LIBMESH_ENABLED
spread_02_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
spread_02_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
spread_02_2d_SOURCES = spread_02.cpp
endif

if LIBMESH_ENABLED
spread_02_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
spread_02_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
spread_02_3d_SOURCES = spread_02.cpp
endif


tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
