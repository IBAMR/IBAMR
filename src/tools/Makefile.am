## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

bin_PROGRAMS = ascii2hdf2d ascii2hdf3d m3D_hdf5_fiber_converter m3D_hdf5_marker_converter m3D_hdf5_strain_analyzer marker_flooder3d

# Dimension-dependent executables
ascii2hdf2d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=2
ascii2hdf2d_SOURCES = ascii2hdf.C

ascii2hdf3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
ascii2hdf3d_SOURCES = ascii2hdf.C

m3D_hdf5_fiber_converter_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
m3D_hdf5_fiber_converter_SOURCES = m3D_hdf5_fiber_converter.C

m3D_hdf5_marker_converter_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
m3D_hdf5_marker_converter_SOURCES = m3D_hdf5_marker_converter.C

m3D_hdf5_strain_analyzer_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
m3D_hdf5_strain_analyzer_SOURCES = m3D_hdf5_strain_analyzer.C

marker_flooder3d_CXXFLAGS= $(AM_CXXFLAGS) -DNDIM=3
marker_flooder3d_SOURCES = marker_flooder.C
