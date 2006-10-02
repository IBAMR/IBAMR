## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

## Dimension-independent library
pkg_include_HEADERS = # intentionally blank

noinst_LTLIBRARIES = libIBAMR_lagrangian.la
libIBAMR_lagrangian_la_CXXFLAGS = $(AM_CXXFLAGS)
libIBAMR_lagrangian_la_LIBADD = fortran/libIBAMR_lagrangian_fortran.la
libIBAMR_lagrangian_la_SOURCES =

## Dimension-dependent libraries
SOURCES =              LDataManager.C LEInteractor.C LNodeIndex.C LNodeIndexData.C LNodeIndexDataFactory.C LNodeIndexSet.C LNodeIndexVariable.C LNodeJacobianInitStrategy.C LNodeLevelData.C LNodePosnInitStrategy.C LagSiloDataWriter.C
pkg_include_HEADERS += LDataManager.h LEInteractor.h LNodeIndex.h LNodeIndexData.h LNodeIndexDataFactory.h LNodeIndexSet.h LNodeIndexVariable.h LNodeJacobianInitStrategy.h LNodeLevelData.h LNodePosnInitStrategy.h LagSiloDataWriter.h
pkg_include_HEADERS +=                               LNodeIndex.I LNodeIndexData.I                         LNodeIndexSet.I

noinst_LTLIBRARIES += libIBAMR_lagrangian2d.la
libIBAMR_lagrangian2d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_lagrangian2d_la_LIBADD = fortran/libIBAMR_lagrangian_fortran2d.la
libIBAMR_lagrangian2d_la_SOURCES = $(SOURCES)

noinst_LTLIBRARIES += libIBAMR_lagrangian3d.la
libIBAMR_lagrangian3d_la_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_lagrangian3d_la_LIBADD = fortran/libIBAMR_lagrangian_fortran3d.la
libIBAMR_lagrangian3d_la_SOURCES = $(SOURCES)
