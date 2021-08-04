## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2021 by the IBAMR developers
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

EXTRA_PROGRAMS = explicit_ex0 explicit_ex1 ib_body_force ib_body_force_kirchhoff

explicit_ex0_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
explicit_ex0_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex0_SOURCES = explicit_ex0.cpp

explicit_ex1_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
explicit_ex1_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex1_SOURCES = explicit_ex1.cpp

ib_body_force_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ib_body_force_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ib_body_force_SOURCES = ib_body_force.cpp

ib_body_force_kirchhoff_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
ib_body_force_kirchhoff_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
ib_body_force_kirchhoff_SOURCES = ib_body_force_kirchhoff.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
