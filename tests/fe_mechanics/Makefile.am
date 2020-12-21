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
EXTRA_PROGRAMS += fe_mechanics_ex0_2d fe_mechanics_ex0_3d
endif

if LIBMESH_ENABLED
fe_mechanics_ex0_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
fe_mechanics_ex0_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
fe_mechanics_ex0_2d_SOURCES = fe_mechanics_ex0.cpp
endif

if LIBMESH_ENABLED
fe_mechanics_ex0_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
fe_mechanics_ex0_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
fe_mechanics_ex0_3d_SOURCES = fe_mechanics_ex0.cpp
endif

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
