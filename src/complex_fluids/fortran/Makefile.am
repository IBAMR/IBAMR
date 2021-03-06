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

EXTRA_DIST =                                                       \
  div_tensor2d.f.m4                div_tensor3d.f.m4               \
  log_upper_convective_op2d.f.m4   log_upper_convective_op3d.f.m4  \
  sqrt_upper_convective_op2d.f.m4  sqrt_upper_convective_op3d.f.m4 \
  upper_convective_op2d.f.m4       upper_convective_op3d.f.m4
BUILT_SOURCES =                                                    \
  div_tensor2d.f                   div_tensor3d.f                  \
  log_upper_convective_op2d.f      log_upper_convective_op3d.f     \
  sqrt_upper_convective_op2d.f     sqrt_upper_convective_op3d.f    \
  upper_convective_op2d.f          upper_convective_op3d.f
CLEANFILES = ${BUILT_SOURCES}
