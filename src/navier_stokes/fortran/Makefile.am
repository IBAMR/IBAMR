## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = \
  navier_stokes_divsource2d.f.m4 navier_stokes_divsource3d.f.m4 \
  navier_stokes_stabledt2d.f.m4  navier_stokes_stabledt3d.f.m4

BUILT_SOURCES = \
  navier_stokes_divsource2d.f navier_stokes_divsource3d.f \
  navier_stokes_stabledt2d.f  navier_stokes_stabledt3d.f
CLEANFILES = ${BUILT_SOURCES}

## Dimension-independent library
noinst_LTLIBRARIES = #libIBAMR_navier_stokes_fortran.la
#libIBAMR_navier_stokes_fortran_la_F77FLAGS = $(AM_F77FLAGS)
#libIBAMR_navier_stokes_fortran_la_SOURCES = # intentionally blank

## Dimension-dependent libraries
noinst_LTLIBRARIES += libIBAMR_navier_stokes_fortran2d.la
libIBAMR_navier_stokes_fortran2d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_navier_stokes_fortran2d_la_SOURCES = \
  navier_stokes_divsource2d.f  navier_stokes_stabledt2d.f

noinst_LTLIBRARIES += libIBAMR_navier_stokes_fortran3d.la
libIBAMR_navier_stokes_fortran3d_la_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_navier_stokes_fortran3d_la_SOURCES = \
  navier_stokes_divsource3d.f  navier_stokes_stabledt3d.f
