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
INSStaggeredHierarchyIntegrator.C \
INSStaggeredLinearOperator.C \
INSStaggeredNonlinearOperator.C \
INSStaggeredProjectionPreconditioner.C \
NormalVelocityBcCoefAdapter.C

pkg_include_HEADERS = \
HierarchyProjector.h \
INSHierarchyIntegrator.h \
INSIntermediateVelocityBcCoef.h \
INSProjectionBcCoef.h \
INSStaggeredHierarchyIntegrator.h \
INSStaggeredLinearOperator.h \
INSStaggeredNonlinearOperator.h \
INSStaggeredProjectionPreconditioner.h \
NormalVelocityBcCoefAdapter.h

libIBAMR_navier_stokes2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR_navier_stokes2d_a_SOURCES = $(SOURCES)

libIBAMR_navier_stokes3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR_navier_stokes3d_a_SOURCES = $(SOURCES)
