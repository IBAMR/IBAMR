## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

bin_PROGRAMS = hdf2vertex2d hdf2vertex3d m3D_hdf5_fiber_converter m3D_hdf5_marker_converter spring2hdf vertex2hdf2d vertex2hdf3d

# Dimension-independent executables
spring2hdf_CXXFLAGS= $(AM_CXXFLAGS)
spring2hdf_SOURCES = spring2hdf.C

# Dimension-dependent executables
hdf2vertex2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
hdf2vertex2d_SOURCES = hdf2vertex.C

hdf2vertex3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
hdf2vertex3d_SOURCES = hdf2vertex.C

m3D_hdf5_fiber_converter_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
m3D_hdf5_fiber_converter_SOURCES = m3D_hdf5_fiber_converter.C

m3D_hdf5_marker_converter_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
m3D_hdf5_marker_converter_SOURCES = m3D_hdf5_marker_converter.C

vertex2hdf2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
vertex2hdf2d_SOURCES = vertex2hdf.C

vertex2hdf3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
vertex2hdf3d_SOURCES = vertex2hdf.C
