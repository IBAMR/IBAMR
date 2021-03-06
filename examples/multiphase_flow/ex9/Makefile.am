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
EXAMPLE_DRIVER = example.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateCircularInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp FlowGravityForcing.cpp
EXTRA_DIST     = input2d.flow input2d.full

EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES) cylinderGen
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)


cylinderGen_CXXFLAGS = $(AM_CXXFLAGS)
cylinderGen_LDADD =
cylinderGen_SOURCES = cylinderGen.cpp

cylinder2d.vertex: cylinderGen
	./cylinderGen

examples: $(EXAMPLES) cylinder2d.vertex
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[2]d* $(PWD) ; \
	  cp -f $(srcdir)/Cylinder2d_cubic.m $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d.flow ; \
	  rm -f $(builddir)/input2d.full ; \
	  rm -f $(builddir)/cylinder2d.vertex ; \
	  rm -f $(builddir)/cylinderGen ; \
	  rm -f $(builddir)/Cylinder2d_cubic.m ; \
	fi ;
