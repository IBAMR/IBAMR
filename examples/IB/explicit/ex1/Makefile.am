## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2019 by the IBAMR developers
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
EXAMPLE_DRIVER = example.cpp
EXTRA_DIST  = input2d input2d.shell input2d.shell_circum_fibers
EXTRA_DIST += input3d
EXTRA_DIST += generate_curve2d.m generate_shell2d.m

EXTRA_DIST += curve2d_64.vertex curve2d_128.vertex curve2d_256.vertex curve2d_512.vertex 
EXTRA_DIST += curve2d_64.spring curve2d_128.spring curve2d_256.spring curve2d_512.spring 
EXTRA_DIST += shell2d_64.vertex shell2d_128.vertex shell2d_256.vertex ##shell2d_512.vertex
EXTRA_DIST += shell2d_64.spring shell2d_128.spring shell2d_256.spring ##shell2d_512.spring
EXTRA_DIST += shell2d_radial_64.vertex shell2d_radial_128.vertex shell2d_radial_256.vertex ##shell2d_radial_512.vertex
EXTRA_DIST += shell2d_radial_64.spring shell2d_radial_128.spring shell2d_radial_256.spring ##shell2d_radial_512.spring
EXTRA_DIST += sphere3d.vertex sphere3d.spring

EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)



examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d{,.shell,.shell_circum_fibers} $(PWD) ; \
	  cp -f $(srcdir)/curve2d_{64,128,256,512}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{64,128,256}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_radial_{64,128,256}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/sphere3d.{vertex,spring} $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d{,.shell,.shell_circum_fibers} ; \
	  rm -f $(builddir)/curve2d_{64,128,256,512}.{vertex,spring} ; \
	  rm -f $(builddir)/shell2d_{64,128,256}.{vertex,spring} ; \
	  rm -f $(builddir)/shell2d_radial_{64,128,256}.{vertex,spring} ; \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/sphere3d.{vertex,spring} ; \
	fi ;
