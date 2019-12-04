## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2019 by the IBAMR developers
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
EXAMPLE_DRIVER = example.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateCircularInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
EXTRA_DIST     = input2d.cyl
EXTRA_DIST    += input3d


EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d
EXTRA_PROGRAMS += $(EXAMPLES) sphereGen
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

sphereGen_CXXFLAGS = $(AM_CXXFLAGS)
sphereGen_LDADD =
sphereGen_SOURCES = sphereGen.cpp



sphere3d.vertex: sphereGen
	./sphereGen

examples: $(EXAMPLES) sphere3d.vertex
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d* $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d.vertex $(PWD) ; \
	  cp -f $(srcdir)/Cylinder2d_cubic.m $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d.cyl ; \
	  rm -f $(builddir)/input3d.sphere ; \
	  rm -f $(builddir)/cylinder2d.vertex ; \
	  rm -f $(builddir)/Cylinder2d_cubic.m ; \
	  rm -f $(builddir)/sphereGen.cpp ; \
	fi ;
