## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

noinst_LIBRARIES =
if SAMRAI2D_ENABLED
noinst_LIBRARIES += libIBAMR_lagrangian2d.a
endif
if SAMRAI3D_ENABLED
noinst_LIBRARIES += libIBAMR_lagrangian3d.a
endif

## Dimension-dependent libraries
SOURCES = \
LDataManager.C \
LEInteractor.C \
LNodeIndex.C \
LNodeIndexData.C \
LNodeIndexDataFactory.C \
LNodeIndexSet.C \
LNodeIndexVariable.C \
LNodeLevelData.C \
LNodeInitStrategy.C \
LagSiloDataWriter.C

pkg_include_HEADERS = \
LDataManager.h \
LEInteractor.h \
LNodeIndex.h \
LNodeIndex.I \
LNodeIndexData.h \
LNodeIndexData.I \
LNodeIndexDataFactory.h \
LNodeIndexSet.h \
LNodeIndexSet.I \
LNodeIndexVariable.h \
LNodeLevelData.h \
LNodeLevelData.I \
LNodeInitStrategy.h \
LagSiloDataWriter.h

libIBAMR_lagrangian2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_lagrangian2d_a_SOURCES = $(SOURCES)

libIBAMR_lagrangian3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_lagrangian3d_a_SOURCES = $(SOURCES)
