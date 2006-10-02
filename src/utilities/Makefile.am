## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

## Dimension-independent library
pkg_include_HEADERS = # intentionally blank

noinst_LTLIBRARIES = #libIBAMR_utilities.la
#libIBAMR_utilities_la_CXXFLAGS = $(AM_CXXFLAGS)
#libIBAMR_utilities_la_LIBADD =
#libIBAMR_utilities_la_SOURCES =

## Dimension-dependent libraries
SOURCES              = ConvergenceMonitor.C PhysicalBCDataStrategy.C SetDataStrategy.C Stashable.C StashableFactory.C StashableManager.C StashableStream.C
pkg_include_HEADERS += ConvergenceMonitor.h PhysicalBCDataStrategy.h SetDataStrategy.h Stashable.h StashableFactory.h StashableManager.h StashableStream.h
pkg_include_HEADERS +=                                                                                                StashableManager.I StashableStream.I

noinst_LTLIBRARIES += libIBAMR_utilities2d.la
libIBAMR_utilities2d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_utilities2d_la_LIBADD =
libIBAMR_utilities2d_la_SOURCES = $(SOURCES)

noinst_LTLIBRARIES += libIBAMR_utilities3d.la
libIBAMR_utilities3d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_utilities3d_la_LIBADD =
libIBAMR_utilities3d_la_SOURCES = $(SOURCES)
