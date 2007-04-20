## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

noinst_LIBRARIES = libIBAMR_IB2d.a libIBAMR_IB3d.a

## Dimension-dependent libraries
SOURCES = \
IBBeamForceGen.C \
IBBeamForceSpec.C \
IBBeamForceSpecFactory.C \
IBEulerianForceSetter.C \
IBEulerianSourceSetter.C \
IBHierarchyIntegrator.C \
IBLagrangianForceStrategy.C \
IBLagrangianForceStrategySet.C \
IBLagrangianSourceStrategy.C \
IBSpringForceGen.C \
IBSpringForceSpec.C \
IBSpringForceSpecFactory.C \
IBStandardForceGen.C \
IBStandardInitializer.C \
IBTargetPointForceGen.C \
IBTargetPointForceSpec.C \
IBTargetPointForceSpecFactory.C

pkg_include_HEADERS = \
IBBeamForceGen.h \
IBBeamForceSpec.h \
IBBeamForceSpecFactory.h \
IBEulerianForceSetter.h \
IBEulerianSourceSetter.h \
IBHierarchyIntegrator.h \
IBLagrangianForceStrategy.h \
IBLagrangianForceStrategySet.h \
IBLagrangianSourceStrategy.h \
IBSpringForceGen.h \
IBSpringForceSpec.h \
IBSpringForceSpecFactory.h \
IBStandardForceGen.h \
IBStandardInitializer.h \
IBTargetPointForceGen.h \
IBTargetPointForceSpec.h \
IBTargetPointForceSpecFactory.h

libIBAMR_IB2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
if SAMRAI2D_ENABLED
libIBAMR_IB2d_a_SOURCES = $(SOURCES)
else
libIBAMR_IB2d_a_SOURCES = # intentionally blank
endif

libIBAMR_IB3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
if SAMRAI3D_ENABLED
libIBAMR_IB3d_a_SOURCES = $(SOURCES)
else
libIBAMR_IB3d_a_SOURCES = # intentionally blank
endif
