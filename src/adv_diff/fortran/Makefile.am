## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = \
  adv_diff_consdiff1d.f.m4 adv_diff_consdiff2d.f.m4 adv_diff_consdiff3d.f.m4

BUILT_SOURCES = \
  adv_diff_consdiff1d.f adv_diff_consdiff2d.f adv_diff_consdiff3d.f
CLEANFILES = ${BUILT_SOURCES}

## Dimension-independent library
noinst_LTLIBRARIES = #libIBAMR_adv_diff_fortran.la
#libIBAMR_adv_diff_fortran_la_F77FLAGS = $(AM_F77FLAGS)
#libIBAMR_adv_diff_fortran_la_SOURCES =

## Dimension-dependent libraries
noinst_LTLIBRARIES += libIBAMR_adv_diff_fortran2d.la
libIBAMR_adv_diff_fortran2d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_adv_diff_fortran2d_la_SOURCES = adv_diff_consdiff2d.f

noinst_LTLIBRARIES += libIBAMR_adv_diff_fortran3d.la
libIBAMR_adv_diff_fortran3d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_adv_diff_fortran3d_la_SOURCES = adv_diff_consdiff3d.f
