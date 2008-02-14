## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = limitertypes.i                          \
  godunov_helpers.f.m4                                \
  advect_detect2d.f.m4      advect_detect3d.f.m4      \
  advect_diff2d.f.m4        advect_diff3d.f.m4        \
  advect_stable2d.f.m4      advect_stable3d.f.m4      \
  godunov_predictors2d.f.m4 godunov_predictors3d.f.m4
BUILT_SOURCES =                                 \
  godunov_helpers.f                             \
  advect_detect2d.f      advect_detect3d.f      \
  advect_diff2d.f        advect_diff3d.f        \
  advect_stable2d.f      advect_stable3d.f      \
  godunov_predictors2d.f godunov_predictors3d.f
CLEANFILES = ${BUILT_SOURCES}

noinst_LIBRARIES = libIBAMR_advect_fortran.a
if SAMRAI2D_ENABLED
noinst_LIBRARIES += libIBAMR_advect_fortran2d.a
endif
if SAMRAI3D_ENABLED
noinst_LIBRARIES += libIBAMR_advect_fortran3d.a
endif

## Dimension-independent library
libIBAMR_advect_fortran_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_advect_fortran_a_SOURCES = godunov_helpers.f

## Dimension-dependent libraries
libIBAMR_advect_fortran2d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_advect_fortran2d_a_SOURCES = advect_detect2d.f advect_diff2d.f advect_stable2d.f godunov_predictors2d.f

libIBAMR_advect_fortran3d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_advect_fortran3d_a_SOURCES = advect_detect3d.f advect_diff3d.f advect_stable3d.f godunov_predictors3d.f
