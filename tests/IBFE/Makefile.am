## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2022 by the IBAMR developers
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

# All programs in this subdirectory depend on libMesh.
EXTRA_PROGRAMS =

if LIBMESH_ENABLED
EXTRA_PROGRAMS += instrument_panel_01_3d interpolate_velocity_01_2d interpolate_velocity_01_3d \
interpolate_velocity_02 explicit_ex0_2d explicit_ex1_2d explicit_ex2_3d explicit_ex4_2d \
explicit_ex4_3d explicit_ex5_2d explicit_ex5_3d explicit_ex8_2d \
ib_partitioning_01_2d ib_partitioning_01_3d ib_partitioning_02_2d \
ib_partitioning_02_3d zero_exterior_values_2d zero_exterior_values_3d

instrument_panel_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
instrument_panel_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
instrument_panel_01_3d_SOURCES = instrument_panel_01.cpp

interpolate_velocity_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
interpolate_velocity_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
interpolate_velocity_01_2d_SOURCES = interpolate_velocity_01.cpp

interpolate_velocity_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
interpolate_velocity_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
interpolate_velocity_01_3d_SOURCES = interpolate_velocity_01.cpp

interpolate_velocity_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
interpolate_velocity_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
interpolate_velocity_02_SOURCES = interpolate_velocity_02.cpp

explicit_ex0_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
explicit_ex0_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex0_2d_SOURCES = explicit_ex0.cpp

explicit_ex1_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
explicit_ex1_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex1_2d_SOURCES = explicit_ex1.cpp

explicit_ex2_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
explicit_ex2_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
explicit_ex2_3d_SOURCES = explicit_ex2.cpp

# this test needs some extra input files, so make SOURCE_DIR available:
explicit_ex4_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
explicit_ex4_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex4_2d_SOURCES = explicit_ex4.cpp

# add it here too, even though we don't use it directly in 3D:
explicit_ex4_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3 -DSOURCE_DIR=\"$(abs_srcdir)\"
explicit_ex4_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
explicit_ex4_3d_SOURCES = explicit_ex4.cpp

explicit_ex5_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
explicit_ex5_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex5_2d_SOURCES = explicit_ex5.cpp

explicit_ex5_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
explicit_ex5_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
explicit_ex5_3d_SOURCES = explicit_ex5.cpp

explicit_ex8_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
explicit_ex8_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex8_2d_SOURCES = explicit_ex8.cpp

ib_partitioning_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ib_partitioning_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ib_partitioning_01_2d_SOURCES = ib_partitioning_01.cpp

ib_partitioning_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
ib_partitioning_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
ib_partitioning_01_3d_SOURCES = ib_partitioning_01.cpp

ib_partitioning_02_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ib_partitioning_02_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ib_partitioning_02_2d_SOURCES = ib_partitioning_02.cpp

ib_partitioning_02_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
ib_partitioning_02_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
ib_partitioning_02_3d_SOURCES = ib_partitioning_02.cpp

zero_exterior_values_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
zero_exterior_values_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
zero_exterior_values_2d_SOURCES = zero_exterior_values.cpp

zero_exterior_values_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
zero_exterior_values_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
zero_exterior_values_3d_SOURCES = zero_exterior_values.cpp
endif

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
