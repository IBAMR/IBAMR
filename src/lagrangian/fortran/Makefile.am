## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = \
  lagrangian_delta.f.m4 lagrangian_interaction2d.f.m4 lagrangian_interaction3d.f.m4

BUILT_SOURCES = \
  lagrangian_delta.f lagrangian_interaction2d.f lagrangian_interaction3d.f
CLEANFILES = ${BUILT_SOURCES}

## Dimension-independent library
noinst_LTLIBRARIES = libIBAMR_lagrangian_fortran.la
libIBAMR_lagrangian_fortran_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_lagrangian_fortran_la_SOURCES = lagrangian_delta.f

## Dimension-dependent libraries
noinst_LTLIBRARIES += libIBAMR_lagrangian_fortran2d.la
libIBAMR_lagrangian_fortran2d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_lagrangian_fortran2d_la_SOURCES = \
  lagrangian_interaction2d.f

noinst_LTLIBRARIES += libIBAMR_lagrangian_fortran3d.la
libIBAMR_lagrangian_fortran3d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_lagrangian_fortran3d_la_SOURCES = \
  lagrangian_interaction3d.f
