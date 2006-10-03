## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

## Dimension-independent library
pkg_include_HEADERS = # intentionally blank

noinst_LTLIBRARIES = #libIBAMR_IB.la
#libIBAMR_IB_la_CXXFLAGS = $(AM_CXXFLAGS)
#libIBAMR_IB_la_LIBADD = fortran/libIBAMR_IB_fortran.la
#libIBAMR_IB_la_SOURCES =

## Dimension-dependent libraries
SOURCES =              IBEulerianForceSetter.C IBEulerianSourceSetter.C IBLagrangianForceStrategy.C IBLagrangianSourceStrategy.C IBLogicalCartMeshInitializer.C SpringForceSpec.C SpringForceSpecFactory.C
pkg_include_HEADERS += IBEulerianForceSetter.h IBEulerianSourceSetter.h IBLagrangianForceStrategy.h IBLagrangianSourceStrategy.h IBLogicalCartMeshInitializer.h SpringForceSpec.h SpringForceSpecFactory.h

noinst_LTLIBRARIES += libIBAMR_IB2d.la
libIBAMR_IB2d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_IB2d_la_LIBADD =
libIBAMR_IB2d_la_SOURCES = $(SOURCES)

noinst_LTLIBRARIES += libIBAMR_IB3d.la
libIBAMR_IB3d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_IB3d_la_LIBADD =
libIBAMR_IB3d_la_SOURCES = $(SOURCES)
