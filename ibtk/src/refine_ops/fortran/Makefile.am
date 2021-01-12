## ---------------------------------------------------------------------
##
## Copyright (c) 2011 - 2020 by the IBAMR developers
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

EXTRA_DIST = cart_side_helpers.f.m4                     \
  cart_side_refine2d.f.m4    cart_side_refine3d.f.m4    \
  divpreservingrefine2d.f.m4 divpreservingrefine3d.f.m4
BUILT_SOURCES = cart_side_helpers.f               \
  cart_side_refine2d.f    cart_side_refine3d.f    \
  divpreservingrefine2d.f divpreservingrefine3d.f
CLEANFILES = ${BUILT_SOURCES}
