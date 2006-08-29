## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran . tests

## Dimension-independent library
pkg_include_HEADERS = # intentionally blank

noinst_LTLIBRARIES = libIBAMR_advect.la
libIBAMR_advect_la_CXXFLAGS = $(AM_CXXFLAGS)
libIBAMR_advect_la_LIBADD = fortran/libIBAMR_advect_fortran.la
libIBAMR_advect_la_SOURCES =

## Dimension-dependent libraries
SOURCES = GodunovAdvector.C GodunovHypPatchOps.C
pkg_include_HEADERS += GodunovAdvector.h GodunovHypPatchOps.h

noinst_LTLIBRARIES += libIBAMR_advect2d.la
libIBAMR_advect2d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_advect2d_la_LIBADD = fortran/libIBAMR_advect_fortran2d.la
libIBAMR_advect2d_la_SOURCES = $(SOURCES)

noinst_LTLIBRARIES += libIBAMR_advect3d.la
libIBAMR_advect3d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_advect3d_la_LIBADD = fortran/libIBAMR_advect_fortran3d.la
libIBAMR_advect3d_la_SOURCES = $(SOURCES)
