## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

noinst_LIBRARIES =
if SAMRAI2D_ENABLED
noinst_LIBRARIES += libIBAMR_IB2d.a
endif
if SAMRAI3D_ENABLED
noinst_LIBRARIES += libIBAMR_IB3d.a
endif

## Dimension-dependent libraries
SOURCES = \
IBBeamForceGen.C \
IBBeamForceSpec.C \
IBBeamForceSpecFactory.C \
IBEulerianForceSetter.C \
IBEulerianSourceSetter.C \
IBHDF5Initializer.C \
IBHierarchyIntegrator.C \
IBInstrumentPanel.C \
IBInstrumentationSpec.C \
IBInstrumentationSpecFactory.C \
IBLagrangianForceStrategy.C \
IBLagrangianForceStrategySet.C \
IBLagrangianSourceStrategy.C \
IBSpringForceGen.C \
IBSpringForceSpec.C \
IBSpringForceSpecFactory.C \
IBStaggeredHierarchyIntegrator.C \
IBStandardForceGen.C \
IBStandardInitializer.C \
IBTargetPointForceGen.C \
IBTargetPointForceSpec.C \
IBTargetPointForceSpecFactory.C

pkg_include_HEADERS = \
IBBeamForceGen.h \
IBBeamForceSpec.h \
IBBeamForceSpec.I \
IBBeamForceSpecFactory.h \
IBEulerianForceSetter.h \
IBEulerianSourceSetter.h \
IBHDF5Initializer.h \
IBHierarchyIntegrator.h \
IBInstrumentPanel.h \
IBInstrumentationSpec.h \
IBInstrumentationSpec.I \
IBInstrumentationSpecFactory.h \
IBLagrangianForceStrategy.h \
IBLagrangianForceStrategySet.h \
IBLagrangianSourceStrategy.h \
IBSpringForceGen.h \
IBSpringForceSpec.h \
IBSpringForceSpec.I \
IBSpringForceSpecFactory.h \
IBStaggeredHierarchyIntegrator.h \
IBStandardForceGen.h \
IBStandardInitializer.h \
IBTargetPointForceGen.h \
IBTargetPointForceSpec.h \
IBTargetPointForceSpec.I \
IBTargetPointForceSpecFactory.h

libIBAMR_IB2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_IB2d_a_SOURCES = $(SOURCES)

libIBAMR_IB3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_IB3d_a_SOURCES = $(SOURCES)
