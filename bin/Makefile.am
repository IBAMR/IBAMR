## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

bin_PROGRAMS = hdf2vertex2d hdf2vertex3d spring2hdf vertex2hdf2d vertex2hdf3d

# Dimension-independent executables
spring2hdf_CXXFLAGS= $(AM_CXXFLAGS)
spring2hdf_SOURCES = spring2hdf.C

# Dimension-dependent executables
hdf2vertex2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
hdf2vertex2d_SOURCES = hdf2vertex.C

hdf2vertex3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
hdf2vertex3d_SOURCES = hdf2vertex.C

vertex2hdf2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
vertex2hdf2d_SOURCES = vertex2hdf.C

vertex2hdf3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
vertex2hdf3d_SOURCES = vertex2hdf.C

## Copy all scripts except for includes.pl to the build directory.
all-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(top_srcdir)/bin/*.pl $(top_builddir)/bin ; \
	  rm -f $(top_builddir)/bin/includes.pl ; \
	fi ;
