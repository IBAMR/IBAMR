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
