## ---------------------------------------------------------------------
##
## Copyright (c) 2006 - 2020 by the IBAMR developers
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

EXTRA_DIST =                                                             \
  advect_helpers2d.f.m4               advect_helpers3d.f.m4              \
  advect_centered_derivatives2d.f.m4  advect_centered_derivatives3d.f.m4 \
  advect_detect2d.f.m4                advect_detect3d.f.m4               \
  advect_diff2d.f.m4                  advect_diff3d.f.m4                 \
  advect_stable2d.f.m4                advect_stable3d.f.m4               \
  advect_predictors2d.f.m4            advect_predictors3d.f.m4
BUILT_SOURCES =                                                   \
  advect_helpers2d.f              advect_helpers3d.f              \
  advect_centered_derivatives2d.f advect_centered_derivatives3d.f \
  advect_detect2d.f               advect_detect3d.f               \
  advect_diff2d.f                 advect_diff3d.f                 \
  advect_stable2d.f               advect_stable3d.f               \
  advect_predictors2d.f           advect_predictors3d.f
CLEANFILES = ${BUILT_SOURCES}
