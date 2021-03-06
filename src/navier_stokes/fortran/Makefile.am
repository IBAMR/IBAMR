## ---------------------------------------------------------------------
##
## Copyright (c) 2006 - 2019 by the IBAMR developers
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

EXTRA_DIST =                                                                                 \
  navier_stokes_bdryop2d.f.m4                     navier_stokes_bdryop3d.f.m4                \
  navier_stokes_divsource2d.f.m4                  navier_stokes_divsource3d.f.m4             \
  navier_stokes_staggered_derivatives2d.f.m4      navier_stokes_staggered_derivatives3d.f.m4 \
  navier_stokes_staggered_helpers2d.f.m4          navier_stokes_staggered_helpers3d.f.m4     \
  navier_stokes_stabledt2d.f.m4                   navier_stokes_stabledt3d.f.m4              \
  navier_stokes_stochastic_forcing2d.f.m4         navier_stokes_stochastic_forcing3d.f.m4    \
  navier_stokes_surface_tension_forcing2d.f.m4    navier_stokes_surface_tension_forcing2d.f.m4
BUILT_SOURCES =                                                                              \
  navier_stokes_bdryop2d.f                        navier_stokes_bdryop3d.f                   \
  navier_stokes_divsource2d.f                     navier_stokes_divsource3d.f                \
  navier_stokes_staggered_derivatives2d.f         navier_stokes_staggered_derivatives3d.f    \
  navier_stokes_staggered_helpers2d.f             navier_stokes_staggered_helpers3d.f        \
  navier_stokes_stabledt2d.f                      navier_stokes_stabledt3d.f                 \
  navier_stokes_stochastic_forcing2d.f            navier_stokes_stochastic_forcing3d.f       \
  navier_stokes_surface_tension_forcing3d.f       navier_stokes_surface_tension_forcing3d.f
CLEANFILES = ${BUILT_SOURCES}
