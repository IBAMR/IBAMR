## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = lagrangian_delta.f.m4 lagrangian_interaction2d.f.m4 lagrangian_interaction3d.f.m4
BUILT_SOURCES = lagrangian_delta.f lagrangian_interaction2d.f lagrangian_interaction3d.f
CLEANFILES = ${BUILT_SOURCES}

noinst_LIBRARIES = libIBAMR_lagrangian_fortran.a libIBAMR_lagrangian_fortran2d.a libIBAMR_lagrangian_fortran3d.a

## Dimension-independent library
libIBAMR_lagrangian_fortran_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_lagrangian_fortran_a_SOURCES = lagrangian_delta.f

## Dimension-dependent libraries
libIBAMR_lagrangian_fortran2d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_lagrangian_fortran2d_a_SOURCES = lagrangian_interaction2d.f

libIBAMR_lagrangian_fortran3d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_lagrangian_fortran3d_a_SOURCES = lagrangian_interaction3d.f
