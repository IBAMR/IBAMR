## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = .

noinst_LIBRARIES = libIBAMR_IB2d.a libIBAMR_IB3d.a

## Dimension-dependent libraries
SOURCES              = IBEulerianForceSetter.C IBEulerianSourceSetter.C IBHierarchyIntegrator.C IBLagrangianForceStrategy.C IBLagrangianSourceStrategy.C IBStandardInitializer.C SpringForceGen.C SpringForceSpec.C SpringForceSpecFactory.C TargetPointForceGen.C TargetPointForceSpec.C TargetPointForceSpecFactory.C TargetPointInitializer.C pIBStandardForceGen.C pIBStandardForceSpec.C pIBStandardForceSpecFactory.C pIBStandardInitializer.C
pkg_include_HEADERS  = IBEulerianForceSetter.h IBEulerianSourceSetter.h IBHierarchyIntegrator.h IBLagrangianForceStrategy.h IBLagrangianSourceStrategy.h IBStandardInitializer.h SpringForceGen.h SpringForceSpec.h SpringForceSpecFactory.h TargetPointForceGen.h TargetPointForceSpec.h TargetPointForceSpecFactory.h TargetPointInitializer.h pIBStandardForceGen.h pIBStandardForceSpec.h pIBStandardForceSpecFactory.h pIBStandardInitializer.h

libIBAMR_IB2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_IB2d_a_SOURCES = $(SOURCES)

libIBAMR_IB3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_IB3d_a_SOURCES = $(SOURCES)
