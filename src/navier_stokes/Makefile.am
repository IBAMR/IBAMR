## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran .

noinst_LIBRARIES =
if SAMRAI2D_ENABLED
noinst_LIBRARIES += libIBAMR_navier_stokes2d.a
endif
if SAMRAI3D_ENABLED
noinst_LIBRARIES += libIBAMR_navier_stokes3d.a
endif

## Dimension-dependent libraries
SOURCES = \
HierarchyProjector.C \
INSHierarchyIntegrator.C \
INSIntermediateVelocityBcCoef.C \
INSProjectionBcCoef.C \
INSStaggeredBlockFactorizationPreconditioner.C \
INSStaggeredCenteredConvectiveOperator.C \
INSStaggeredHierarchyIntegrator.C \
INSStaggeredIntermediateVelocityBcCoef.C \
INSStaggeredPPMConvectiveOperator.C \
INSStaggeredPhysicalBoundaryHelper.C \
INSStaggeredPressureBcCoef.C \
INSStaggeredProjectionBcCoef.C \
INSStaggeredProjectionPreconditioner.C \
INSStaggeredStokesOperator.C \
INSStaggeredVelocityBcCoef.C

pkg_include_HEADERS = \
HierarchyProjector.h \
INSHierarchyIntegrator.h \
INSIntermediateVelocityBcCoef.h \
INSProjectionBcCoef.h \
INSStaggeredBlockFactorizationPreconditioner.h \
INSStaggeredCenteredConvectiveOperator.h \
INSStaggeredHierarchyIntegrator.h \
INSStaggeredIntermediateVelocityBcCoef.h \
INSStaggeredPPMConvectiveOperator.h \
INSStaggeredPhysicalBoundaryHelper.h \
INSStaggeredPressureBcCoef.h \
INSStaggeredProjectionBcCoef.h \
INSStaggeredProjectionPreconditioner.h \
INSStaggeredStokesOperator.h \
INSStaggeredVelocityBcCoef.h

libIBAMR_navier_stokes2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_navier_stokes2d_a_SOURCES = $(SOURCES)

libIBAMR_navier_stokes3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_navier_stokes3d_a_SOURCES = $(SOURCES)
