## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST =                                                                                 \
  navier_stokes_bdryop2d.f.m4                     navier_stokes_bdryop3d.f.m4                \
  navier_stokes_divsource2d.f.m4                  navier_stokes_divsource3d.f.m4             \
  navier_stokes_staggered_derivatives2d.f.m4      navier_stokes_staggered_derivatives3d.f.m4 \
  navier_stokes_staggered_helpers2d.f.m4          navier_stokes_staggered_helpers3d.f.m4     \
  navier_stokes_stabledt2d.f.m4                   navier_stokes_stabledt3d.f.m4              \
  navier_stokes_stochastic_forcing2d.f.m4         navier_stokes_stochastic_forcing3d.f.m4    \
  navier_stokes_surface_tension_forcing2d.f.m4    navier_stokes_surface_tension_forcing2d.f
BUILT_SOURCES =                                                                              \
  navier_stokes_bdryop2d.f                        navier_stokes_bdryop3d.f                   \
  navier_stokes_divsource2d.f                     navier_stokes_divsource3d.f                \
  navier_stokes_staggered_derivatives2d.f         navier_stokes_staggered_derivatives3d.f    \
  navier_stokes_staggered_helpers2d.f             navier_stokes_staggered_helpers3d.f        \
  navier_stokes_stabledt2d.f                      navier_stokes_stabledt3d.f                 \
  navier_stokes_stochastic_forcing2d.f            navier_stokes_stochastic_forcing3d.f       \
  navier_stokes_surface_tension_forcing3d.f.m4    navier_stokes_surface_tension_forcing3d.f
CLEANFILES = ${BUILT_SOURCES}
