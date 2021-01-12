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

EXTRA_PROGRAMS = free_falling_cyl_cib rotating_barge_cib check_hydro_force high_density_droplet

# this test needs some extra input files, so make SOURCE_DIR available:
free_falling_cyl_cib_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
free_falling_cyl_cib_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
free_falling_cyl_cib_SOURCES = free_falling_cyl_cib.cpp

rotating_barge_cib_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
rotating_barge_cib_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
rotating_barge_cib_SOURCES = rotating_barge_cib.cpp

check_hydro_force_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
check_hydro_force_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
check_hydro_force_SOURCES = check_hydro_force.cpp

# this test needs some extra input files, so make SOURCE_DIR available:
high_density_droplet_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
high_density_droplet_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
high_density_droplet_SOURCES = high_density_droplet.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
