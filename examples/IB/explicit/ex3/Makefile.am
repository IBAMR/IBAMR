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
EXTRA_DIST  = input2d
EXTRA_DIST +=  fila_256.vertex  fila_256.spring  fila_256.beam  fila_256.target
EXTRA_DIST +=  fila_512.vertex  fila_512.spring  fila_512.beam  fila_512.target
EXTRA_DIST += fila_1024.vertex fila_1024.spring fila_1024.beam fila_1024.target

EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
EXTRA_PROGRAMS += $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)



examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d* $(PWD) ; \
	  cp -f $(srcdir)/fila_{256,512,1024}.{vertex,spring,beam,target} $(PWD) ; \
	fi ;


clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/fila_{256,512,1024}.{vertex,spring,beam,target} ; \
	fi ;
