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

## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = example.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateStructureInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp FlowGravityForcing.cpp LevelSetGasInitialCondition.cpp LevelSetSolidInitialCondition.cpp
EXTRA_DIST     = input3d

EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXAMPLES += main3d
EXTRA_PROGRAMS += $(EXAMPLES) rectangleGen
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)


rectangleGen_CXXFLAGS = $(AM_CXXFLAGS)
rectangleGen_LDADD =
rectangleGen_SOURCES = rectangleGen.cpp

rectangle3d.vertex: rectangleGen
	./rectangleGen

examples: $(EXAMPLES) rectangle3d.vertex
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[3]d* $(PWD) ; \
	  cp -f $(srcdir)/Rectangle3d.m $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/rectangle3d.vertex ; \
	  rm -f $(builddir)/rectangleGen ; \
	  rm -f $(builddir)/Rectangle3d.m ; \
	fi ;
