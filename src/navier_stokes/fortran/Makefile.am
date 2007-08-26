## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST =                                                    \
  navier_stokes_bdryop2d.f.m4    navier_stokes_bdryop3d.f.m4    \
  navier_stokes_divsource2d.f.m4 navier_stokes_divsource3d.f.m4 \
  navier_stokes_stabledt2d.f.m4  navier_stokes_stabledt3d.f.m4
BUILT_SOURCES =                                           \
  navier_stokes_bdryop2d.f    navier_stokes_bdryop3d.f    \
  navier_stokes_divsource2d.f navier_stokes_divsource3d.f \
  navier_stokes_stabledt2d.f  navier_stokes_stabledt3d.f
CLEANFILES = ${BUILT_SOURCES}

noinst_LIBRARIES =
if SAMRAI2D_ENABLED
noinst_LIBRARIES += libIBAMR_navier_stokes_fortran2d.a
endif
if SAMRAI3D_ENABLED
noinst_LIBRARIES += libIBAMR_navier_stokes_fortran3d.a
endif

## Dimension-dependent libraries
libIBAMR_navier_stokes_fortran2d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_navier_stokes_fortran2d_a_SOURCES = navier_stokes_bdryop2d.f navier_stokes_divsource2d.f navier_stokes_stabledt2d.f

libIBAMR_navier_stokes_fortran3d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_navier_stokes_fortran3d_a_SOURCES = navier_stokes_bdryop2d.f navier_stokes_divsource3d.f navier_stokes_stabledt3d.f
