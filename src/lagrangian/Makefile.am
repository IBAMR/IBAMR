## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

noinst_LIBRARIES = libIBAMR_lagrangian2d.a libIBAMR_lagrangian3d.a

## Dimension-dependent libraries
SOURCES              = LDataManager.C LEInteractor.C LNodeIndex.C LNodeIndexData.C LNodeIndexDataFactory.C LNodeIndexSet.C LNodeIndexVariable.C LNodeJacobianInitStrategy.C LNodeLevelData.C LNodePosnInitStrategy.C LagSiloDataWriter.C
pkg_include_HEADERS  = LDataManager.h LEInteractor.h LNodeIndex.h LNodeIndexData.h LNodeIndexDataFactory.h LNodeIndexSet.h LNodeIndexVariable.h LNodeJacobianInitStrategy.h LNodeLevelData.h LNodePosnInitStrategy.h LagSiloDataWriter.h
pkg_include_HEADERS +=                               LNodeIndex.I LNodeIndexData.I                         LNodeIndexSet.I                                                  LNodeLevelData.I

libIBAMR_lagrangian2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_lagrangian2d_a_SOURCES = $(SOURCES)

libIBAMR_lagrangian3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_lagrangian3d_a_SOURCES = $(SOURCES)
