## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = IB fortran adv_diff advect navier_stokes utilities . tools

lib_LIBRARIES = libIBAMR.a
if SAMRAI2D_ENABLED
lib_LIBRARIES += libIBAMR2d.a
endif
if SAMRAI3D_ENABLED
lib_LIBRARIES += libIBAMR3d.a
endif

## Dimension-independent library
libIBAMR_a_SOURCES = \
get_version.C \
advect/fortran/godunov_helpers.f \
adv_diff/TGACoefs.C

pkg_include_HEADERS = \
get_version.h \
adv_diff/TGACoefs.h

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
IB/IBAnchorPointSpec.C \
IB/IBAnchorPointSpecFactory.C \
IB/IBBeamForceGen.C \
IB/IBBeamForceSpec.C \
IB/IBBeamForceSpecFactory.C \
IB/IBDataPostProcessor.C \
IB/IBEulerianForceSetter.C \
IB/IBEulerianSourceSetter.C \
IB/IBHDF5Initializer.C \
IB/IBHierarchyIntegrator.C \
IB/IBImplicitHierarchyIntegrator.C \
IB/IBInstrumentPanel.C \
IB/IBInstrumentationSpec.C \
IB/IBInstrumentationSpecFactory.C \
IB/IBLagrangianForceStrategy.C \
IB/IBLagrangianForceStrategySet.C \
IB/IBLagrangianSourceStrategy.C \
IB/IBMovingTargetPointForceGen.C \
IB/IBMovingTargetPointForceSpec.C \
IB/IBMovingTargetPointForceSpecFactory.C \
IB/IBSpringForceGen.C \
IB/IBSpringForceSpec.C \
IB/IBSpringForceSpecFactory.C \
IB/IBStaggeredHierarchyIntegrator.C \
IB/IBStandardForceGen.C \
IB/IBStandardInitializer.C \
IB/IBTargetPointForceGen.C \
IB/IBTargetPointForceSpec.C \
IB/IBTargetPointForceSpecFactory.C \
adv_diff/AdvDiffHierarchyIntegrator.C \
adv_diff/AdvDiffHypPatchOps.C \
advect/AdvectHypPatchOps.C \
advect/GodunovAdvector.C \
navier_stokes/HierarchyProjector.C \
navier_stokes/INSCoefs.C \
navier_stokes/INSHierarchyIntegrator.C \
navier_stokes/INSIntermediateVelocityBcCoef.C \
navier_stokes/INSProjectionBcCoef.C \
navier_stokes/INSStaggeredBlockFactorizationPreconditioner.C \
navier_stokes/INSStaggeredCenteredConvectiveOperator.C \
navier_stokes/INSStaggeredHierarchyIntegrator.C \
navier_stokes/INSStaggeredIntermediateVelocityBcCoef.C \
navier_stokes/INSStaggeredPPMConvectiveOperator.C \
navier_stokes/INSStaggeredPhysicalBoundaryHelper.C \
navier_stokes/INSStaggeredPressureBcCoef.C \
navier_stokes/INSStaggeredProjectionBcCoef.C \
navier_stokes/INSStaggeredProjectionPreconditioner.C \
navier_stokes/INSStaggeredStokesOperator.C \
navier_stokes/INSStaggeredVelocityBcCoef.C

pkg_include_HEADERS += \
IB/IBAnchorPointSpec.I \
IB/IBAnchorPointSpec.h \
IB/IBAnchorPointSpecFactory.h \
IB/IBBeamForceGen.h \
IB/IBBeamForceSpec.I \
IB/IBBeamForceSpec.h \
IB/IBBeamForceSpecFactory.h \
IB/IBDataPostProcessor.h \
IB/IBEulerianForceSetter.h \
IB/IBEulerianSourceSetter.h \
IB/IBHDF5Initializer.h \
IB/IBHierarchyIntegrator.h \
IB/IBImplicitHierarchyIntegrator.h \
IB/IBInstrumentPanel.h \
IB/IBInstrumentationSpec.I \
IB/IBInstrumentationSpec.h \
IB/IBInstrumentationSpecFactory.h \
IB/IBLagrangianForceStrategy.h \
IB/IBLagrangianForceStrategySet.h \
IB/IBLagrangianSourceStrategy.h \
IB/IBMovingTargetPointForceGen.h \
IB/IBMovingTargetPointForceSpec.I \
IB/IBMovingTargetPointForceSpec.h \
IB/IBMovingTargetPointForceSpecFactory.h \
IB/IBSpringForceGen.h \
IB/IBSpringForceSpec.I \
IB/IBSpringForceSpec.h \
IB/IBSpringForceSpecFactory.h \
IB/IBStaggeredHierarchyIntegrator.h \
IB/IBStandardForceGen.h \
IB/IBStandardInitializer.h \
IB/IBTargetPointForceGen.h \
IB/IBTargetPointForceSpec.I \
IB/IBTargetPointForceSpec.h \
IB/IBTargetPointForceSpecFactory.h \
adv_diff/AdvDiffHierarchyIntegrator.h \
adv_diff/AdvDiffHypPatchOps.h \
advect/AdvectHypPatchOps.h \
advect/GodunovAdvector.h \
navier_stokes/HierarchyProjector.h \
navier_stokes/INSCoefs.h \
navier_stokes/INSHierarchyIntegrator.h \
navier_stokes/INSIntermediateVelocityBcCoef.h \
navier_stokes/INSProjectionBcCoef.h \
navier_stokes/INSStaggeredBlockFactorizationPreconditioner.h \
navier_stokes/INSStaggeredCenteredConvectiveOperator.h \
navier_stokes/INSStaggeredHierarchyIntegrator.h \
navier_stokes/INSStaggeredIntermediateVelocityBcCoef.h \
navier_stokes/INSStaggeredPPMConvectiveOperator.h \
navier_stokes/INSStaggeredPhysicalBoundaryHelper.h \
navier_stokes/INSStaggeredPressureBcCoef.h \
navier_stokes/INSStaggeredProjectionBcCoef.h \
navier_stokes/INSStaggeredProjectionPreconditioner.h \
navier_stokes/INSStaggeredStokesOperator.h \
navier_stokes/INSStaggeredVelocityBcCoef.h

libIBAMR2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
adv_diff/fortran/adv_diff_consdiff2d.f \
advect/fortran/advect_detect2d.f \
advect/fortran/advect_diff2d.f \
advect/fortran/advect_stable2d.f \
advect/fortran/godunov_predictors2d.f \
advect/fortran/weno_fluxes2d.f \
navier_stokes/fortran/navier_stokes_bdryop2d.f \
navier_stokes/fortran/navier_stokes_divsource2d.f \
navier_stokes/fortran/navier_stokes_godunov_predictors2d.f \
navier_stokes/fortran/navier_stokes_stabledt2d.f \
navier_stokes/fortran/navier_stokes_staggered_derivatives2d.f \
navier_stokes/fortran/navier_stokes_staggered_helpers2d.f

libIBAMR3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
adv_diff/fortran/adv_diff_consdiff3d.f \
advect/fortran/advect_detect3d.f \
advect/fortran/advect_diff3d.f \
advect/fortran/advect_stable3d.f \
advect/fortran/godunov_predictors3d.f \
navier_stokes/fortran/navier_stokes_bdryop3d.f \
navier_stokes/fortran/navier_stokes_divsource3d.f \
navier_stokes/fortran/navier_stokes_stabledt3d.f \
navier_stokes/fortran/navier_stokes_staggered_derivatives3d.f \
navier_stokes/fortran/navier_stokes_staggered_helpers3d.f
