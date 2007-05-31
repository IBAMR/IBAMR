## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

noinst_LIBRARIES =
if SAMRAI2D_ENABLED
noinst_LIBRARIES += libIBAMR_utilities2d.a
endif
if SAMRAI3D_ENABLED
noinst_LIBRARIES += libIBAMR_utilities3d.a
endif

## Dimension-dependent libraries
SOURCES = \
Stashable.C \
StashableFactory.C \
StashableManager.C \
StashableStream.C \
get_version.C

pkg_include_HEADERS = \
Stashable.h \
StashableFactory.h \
StashableManager.h \
StashableManager.I \
StashableStream.h \
StashableStream.I \
get_version.h

libIBAMR_utilities2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_utilities2d_a_SOURCES = $(SOURCES)

libIBAMR_utilities3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_utilities3d_a_SOURCES = $(SOURCES)

# Dimension-independent executables
bin_PROGRAMS  = spring2hdf

spring2hdf_CXXFLAGS= $(AM_CXXFLAGS)
spring2hdf_SOURCES = spring2hdf.C

# Dimension-dependent executables
bin_PROGRAMS += hdf2vertex2d hdf2vertex3d vertex2hdf2d vertex2hdf3d

hdf2vertex2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
hdf2vertex2d_SOURCES = hdf2vertex.C

hdf2vertex3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
hdf2vertex3d_SOURCES = hdf2vertex.C

vertex2hdf2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
vertex2hdf2d_SOURCES = vertex2hdf.C

vertex2hdf3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
vertex2hdf3d_SOURCES = vertex2hdf.C
