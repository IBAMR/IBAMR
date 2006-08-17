## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = limitertypes.i \
  godunov_helpers.f.m4 \
  advect_detect1d.f.m4      advect_detect2d.f.m4      advect_detect3d.f.m4     \
  advect_diff1d.f.m4        advect_diff2d.f.m4        advect_diff3d.f.m4       \
  advect_stable1d.f.m4      advect_stable2d.f.m4      advect_stable3d.f.m4     \
  godunov_predictors1d.f.m4 godunov_predictors2d.f.m4 godunov_predictors3d.f.m4

BUILT_SOURCES = \
  godunov_helpers.f \
  advect_detect1d.f      advect_detect2d.f      advect_detect3d.f     \
  advect_diff1d.f        advect_diff2d.f        advect_diff3d.f       \
  advect_stable1d.f      advect_stable2d.f      advect_stable3d.f     \
  godunov_predictors1d.f godunov_predictors2d.f godunov_predictors3d.f
CLEANFILES = ${BUILT_SOURCES}

noinst_LTLIBRARIES = libIBAMR_advect_fortran.la
libIBAMR_advect_fortran_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_advect_fortran_la_SOURCES = \
  godunov_helpers.f

noinst_LTLIBRARIES += libIBAMR_advect_fortran2d.la
libIBAMR_advect_fortran2d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_advect_fortran2d_la_SOURCES = \
  advect_detect2d.f      \
  advect_diff2d.f        \
  advect_stable2d.f      \
  godunov_predictors2d.f

noinst_LTLIBRARIES += libIBAMR_advect_fortran3d.la
libIBAMR_advect_fortran3d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_advect_fortran3d_la_SOURCES = \
  advect_detect3d.f      \
  advect_diff3d.f        \
  advect_stable3d.f      \
  godunov_predictors3d.f