## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran . # tests

## Dimension-independent library
pkg_include_HEADERS = # intentionally blank

noinst_LTLIBRARIES = #libIBAMR_navier_stokes.la
#libIBAMR_navier_stokes_la_CXXFLAGS = $(AM_CXXFLAGS)
#libIBAMR_navier_stokes_la_LIBADD = fortran/libIBAMR_navier_stokes_fortran.la
#libIBAMR_navier_stokes_la_SOURCES =

## Dimension-dependent libraries
SOURCES = HierarchyProjector.C INSHierarchyIntegrator.C
pkg_include_HEADERS += HierarchyProjector.h INSHierarchyIntegrator.h

noinst_LTLIBRARIES += libIBAMR_navier_stokes2d.la
libIBAMR_navier_stokes2d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_navier_stokes2d_la_LIBADD = fortran/libIBAMR_navier_stokes_fortran2d.la
libIBAMR_navier_stokes2d_la_SOURCES = $(SOURCES)

noinst_LTLIBRARIES += libIBAMR_navier_stokes3d.la
libIBAMR_navier_stokes3d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_navier_stokes3d_la_LIBADD = fortran/libIBAMR_navier_stokes_fortran3d.la
libIBAMR_navier_stokes3d_la_SOURCES = $(SOURCES)
