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

EXTRA_PROGRAMS = fo_acoustic_streaming_solver_2d fo_acoustic_streaming_solver_3d \
acoustic_streaming_hier_integrator_2d acoustic_streaming_hier_integrator_3d \
so_boundary_conditions_2d sd_mass_src_2d

fo_acoustic_streaming_solver_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
fo_acoustic_streaming_solver_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
fo_acoustic_streaming_solver_2d_SOURCES = fo_acoustic_streaming_solver.cpp

fo_acoustic_streaming_solver_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
fo_acoustic_streaming_solver_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
fo_acoustic_streaming_solver_3d_SOURCES = fo_acoustic_streaming_solver.cpp

acoustic_streaming_hier_integrator_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
acoustic_streaming_hier_integrator_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
acoustic_streaming_hier_integrator_2d_SOURCES = acoustic_streaming_hier_integrator.cpp

acoustic_streaming_hier_integrator_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
acoustic_streaming_hier_integrator_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
acoustic_streaming_hier_integrator_3d_SOURCES = acoustic_streaming_hier_integrator.cpp

so_boundary_conditions_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
so_boundary_conditions_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
so_boundary_conditions_2d_SOURCES = so_boundary_conditions.cpp

sd_mass_src_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
sd_mass_src_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
sd_mass_src_2d_SOURCES = sd_mass_src.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
