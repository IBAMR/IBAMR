## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran adv_diff advect navier_stokes utilities .

## Dimension-independent library
lib_LTLIBRARIES = libIBAMR.la
libIBAMR_la_LIBADD = $(FLIBS)			\
  $(top_builddir)/stools/src/libSTOOLS.la	\
  adv_diff/libIBAMR_adv_diff.la			\
  advect/libIBAMR_advect.la			\
  utilities/libIBAMR_utilities.la
libIBAMR_la_LDFLAGS = -version 0:0:0
libIBAMR_la_SOURCES =

## Dimension-dependent libraries
lib_LTLIBRARIES += libIBAMR2d.la
libIBAMR2d_la_LIBADD = $(FLIBS)			\
  $(top_builddir)/stools/src/libSTOOLS2d.la	\
  adv_diff/libIBAMR_adv_diff2d.la		\
  advect/libIBAMR_advect2d.la			\
  utilities/libIBAMR_utilities2d.la
libIBAMR2d_la_LDFLAGS = -version 0:0:0
libIBAMR2d_la_SOURCES =

lib_LTLIBRARIES += libIBAMR3d.la
libIBAMR3d_la_LIBADD = $(FLIBS)			\
  $(top_builddir)/stools/src/libSTOOLS3d.la	\
  adv_diff/libIBAMR_adv_diff3d.la		\
  advect/libIBAMR_advect3d.la			\
  utilities/libIBAMR_utilities3d.la
libIBAMR3d_la_LDFLAGS = -version 0:0:0
libIBAMR3d_la_SOURCES =
