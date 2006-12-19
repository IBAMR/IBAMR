## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

noinst_LIBRARIES = libIBAMR_advect2d.a libIBAMR_advect3d.a

## Dimension-dependent libraries
SOURCES              = GodunovAdvector.C GodunovHypPatchOps.C
pkg_include_HEADERS  = GodunovAdvector.h GodunovHypPatchOps.h

libIBAMR_advect2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_advect2d_a_SOURCES = $(SOURCES)

libIBAMR_advect3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_advect3d_a_SOURCES = $(SOURCES)
