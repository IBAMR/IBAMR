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
LNodeIndexData2.C \
LNodeIndexDataFactory2.C \
LNodeIndexSet.C \
LNodeIndexVariable2.C \
LNodeLevelData.C \
LNodeInitStrategy.C \
LagSiloDataWriter.C

pkg_include_HEADERS = \
ArrayData_specialized_LNodeIndexSet.I \
LDataManager.h \
LEInteractor.h \
LNodeIndex.h \
LNodeIndex.I \
LNodeIndexData2.h \
LNodeIndexData2.I \
LNodeIndexDataFactory2.h \
LNodeIndexSet.h \
LNodeIndexSet.I \
LNodeIndexVariable2.h \
LNodeLevelData.h \
LNodeLevelData.I \
LNodeInitStrategy.h \
LagM3dDataWriter.h \
LagSiloDataWriter.h

libIBAMR_lagrangian2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_lagrangian2d_a_SOURCES = $(SOURCES)

libIBAMR_lagrangian3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_lagrangian3d_a_SOURCES = $(SOURCES) LagM3dDataWriter.C
