## ---------------------------------------------------------------------
##
## Copyright (c) 2015 - 2021 by the IBAMR developers
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
EXTRA_PROGRAMS += poiseuille_flow hagen_poiseuille_flow \
taylor_couette_2d taylor_couette_3d flow_past_cylinder \
flow_past_sphere
endif

if LIBMESH_ENABLED
poiseuille_flow_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
poiseuille_flow_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
poiseuille_flow_SOURCES = poiseuille_flow.cpp
endif

if LIBMESH_ENABLED
hagen_poiseuille_flow_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
hagen_poiseuille_flow_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
hagen_poiseuille_flow_SOURCES = hagen_poiseuille_flow.cpp
endif


if LIBMESH_ENABLED
taylor_couette_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
taylor_couette_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
taylor_couette_2d_SOURCES = taylor_couette_2d.cpp
endif

if LIBMESH_ENABLED
taylor_couette_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
taylor_couette_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
taylor_couette_3d_SOURCES = taylor_couette_3d.cpp
endif

if LIBMESH_ENABLED
flow_past_cylinder_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
flow_past_cylinder_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
flow_past_cylinder_SOURCES = flow_past_cylinder.cpp
endif

if LIBMESH_ENABLED
flow_past_sphere_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
flow_past_sphere_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
flow_past_sphere_SOURCES = flow_past_sphere.cpp
endif


tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
