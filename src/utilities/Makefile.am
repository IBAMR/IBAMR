## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

noinst_LIBRARIES = libIBAMR_utilities2d.a libIBAMR_utilities3d.a

## Dimension-dependent libraries
SOURCES              = ConvergenceMonitor.C PhysicalBCDataStrategy.C SetDataStrategy.C Stashable.C StashableFactory.C StashableManager.C StashableStream.C get_version.C
pkg_include_HEADERS  = ConvergenceMonitor.h PhysicalBCDataStrategy.h SetDataStrategy.h Stashable.h StashableFactory.h StashableManager.h StashableStream.h get_version.h
pkg_include_HEADERS +=                                                                                                StashableManager.I StashableStream.I

libIBAMR_utilities2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_utilities2d_a_SOURCES = $(SOURCES)

libIBAMR_utilities3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_utilities3d_a_SOURCES = $(SOURCES)
