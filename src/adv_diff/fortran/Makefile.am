## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = adv_diff_consdiff1d.f.m4 adv_diff_consdiff2d.f.m4 adv_diff_consdiff3d.f.m4
BUILT_SOURCES = adv_diff_consdiff1d.f adv_diff_consdiff2d.f adv_diff_consdiff3d.f
CLEANFILES = ${BUILT_SOURCES}

noinst_LIBRARIES = libIBAMR_adv_diff_fortran2d.a libIBAMR_adv_diff_fortran3d.a

## Dimension-dependent libraries
libIBAMR_adv_diff_fortran2d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_adv_diff_fortran2d_a_SOURCES = adv_diff_consdiff2d.f

libIBAMR_adv_diff_fortran3d_a_F77FLAGS = $(AM_F77FLAGS)
libIBAMR_adv_diff_fortran3d_a_SOURCES = adv_diff_consdiff3d.f
