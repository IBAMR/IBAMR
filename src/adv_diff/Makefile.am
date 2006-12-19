## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

noinst_LIBRARIES = libIBAMR_adv_diff2d.a libIBAMR_adv_diff3d.a

## Dimension-dependent libraries
SOURCES              = AdvDiffHierarchyIntegrator.C AdvDiffHypPatchOps.C
pkg_include_HEADERS  = AdvDiffHierarchyIntegrator.h AdvDiffHypPatchOps.h

libIBAMR_adv_diff2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_adv_diff2d_a_SOURCES = $(SOURCES)

libIBAMR_adv_diff3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_adv_diff3d_a_SOURCES = $(SOURCES)
