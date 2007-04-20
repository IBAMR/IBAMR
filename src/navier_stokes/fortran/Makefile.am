## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST =                                                    \
  navier_stokes_divsource2d.f.m4 navier_stokes_divsource3d.f.m4 \
  navier_stokes_stabledt2d.f.m4  navier_stokes_stabledt3d.f.m4
BUILT_SOURCES =                                           \
  navier_stokes_divsource2d.f navier_stokes_divsource3d.f \
  navier_stokes_stabledt2d.f  navier_stokes_stabledt3d.f
CLEANFILES = ${BUILT_SOURCES}

noinst_LIBRARIES = libIBAMR_navier_stokes_fortran2d.a libIBAMR_navier_stokes_fortran3d.a

## Dimension-dependent libraries
libIBAMR_navier_stokes_fortran2d_a_F77FLAGS = $(AM_F77FLAGS)
if SAMRAI2D_ENABLED
libIBAMR_navier_stokes_fortran2d_a_SOURCES = navier_stokes_divsource2d.f  navier_stokes_stabledt2d.f
else
libIBAMR_navier_stokes_fortran2d_a_SOURCES = # intentionally blank
endif

libIBAMR_navier_stokes_fortran3d_a_F77FLAGS = $(AM_F77FLAGS)
if SAMRAI3D_ENABLED
libIBAMR_navier_stokes_fortran3d_a_SOURCES = navier_stokes_divsource3d.f  navier_stokes_stabledt3d.f
else
libIBAMR_navier_stokes_fortran3d_a_SOURCES = # intentionally blank
endif
