## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

## Dimension-independent library
pkg_include_HEADERS = # intentionally blank

noinst_LTLIBRARIES = #libIBAMR_adv_diff.la
#libIBAMR_adv_diff_la_CXXFLAGS = $(AM_CXXFLAGS)
#libIBAMR_adv_diff_la_LIBADD = fortran/libIBAMR_adv_diff_fortran.la
#libIBAMR_adv_diff_la_SOURCES =

## Dimension-dependent libraries
SOURCES = AdvDiffHierarchyIntegrator.C AdvDiffHypPatchOps.C
pkg_include_HEADERS += AdvDiffHierarchyIntegrator.h AdvDiffHypPatchOps.h

noinst_LTLIBRARIES += libIBAMR_adv_diff2d.la
libIBAMR_adv_diff2d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_adv_diff2d_la_LIBADD = fortran/libIBAMR_adv_diff_fortran2d.la
libIBAMR_adv_diff2d_la_SOURCES = $(SOURCES)

noinst_LTLIBRARIES += libIBAMR_adv_diff3d.la
libIBAMR_adv_diff3d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_adv_diff3d_la_LIBADD = fortran/libIBAMR_adv_diff_fortran3d.la
libIBAMR_adv_diff3d_la_SOURCES = $(SOURCES)
