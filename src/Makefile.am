## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = IB fortran adv_diff advect navier_stokes utilities . tests tools

lib_LIBRARIES = libIBAMR.a
if SAMRAI2D_ENABLED
lib_LIBRARIES += libIBAMR2d.a
endif
if SAMRAI3D_ENABLED
lib_LIBRARIES += libIBAMR3d.a
endif

## Dimension-independent library
libIBAMR_a_LIBADD = \
advect/fortran/godunov_helpers.o \
adv_diff/libIBAMR_adv_diff_a-TGACoefs.o
libIBAMR_a_SOURCES = # intentionally blank

## Dimension-dependent libraries
libIBAMR2d_a_LIBADD = \
IB/libIBAMR_IB2d_a-IBAnchorPointSpec.o \
IB/libIBAMR_IB2d_a-IBAnchorPointSpecFactory.o \
IB/libIBAMR_IB2d_a-IBBeamForceGen.o \
IB/libIBAMR_IB2d_a-IBBeamForceSpec.o \
IB/libIBAMR_IB2d_a-IBBeamForceSpecFactory.o \
IB/libIBAMR_IB2d_a-IBEulerianForceSetter.o \
IB/libIBAMR_IB2d_a-IBEulerianSourceSetter.o \
IB/libIBAMR_IB2d_a-IBHDF5Initializer.o \
IB/libIBAMR_IB2d_a-IBHierarchyIntegrator.o \
IB/libIBAMR_IB2d_a-IBInstrumentPanel.o \
IB/libIBAMR_IB2d_a-IBInstrumentationSpec.o \
IB/libIBAMR_IB2d_a-IBInstrumentationSpecFactory.o \
IB/libIBAMR_IB2d_a-IBLagrangianForceStrategy.o \
IB/libIBAMR_IB2d_a-IBLagrangianForceStrategySet.o \
IB/libIBAMR_IB2d_a-IBLagrangianSourceStrategy.o \
IB/libIBAMR_IB2d_a-IBMovingTargetPointForceGen.o \
IB/libIBAMR_IB2d_a-IBMovingTargetPointForceSpecFactory.o \
IB/libIBAMR_IB2d_a-IBMovingTargetPointForceSpec.o \
IB/libIBAMR_IB2d_a-IBSpringForceGen.o \
IB/libIBAMR_IB2d_a-IBSpringForceSpecFactory.o \
IB/libIBAMR_IB2d_a-IBSpringForceSpec.o \
IB/libIBAMR_IB2d_a-IBStaggeredHierarchyIntegrator.o \
IB/libIBAMR_IB2d_a-IBStandardForceGen.o \
IB/libIBAMR_IB2d_a-IBStandardInitializer.o \
IB/libIBAMR_IB2d_a-IBTargetPointForceGen.o \
IB/libIBAMR_IB2d_a-IBTargetPointForceSpecFactory.o \
IB/libIBAMR_IB2d_a-IBTargetPointForceSpec.o \
adv_diff/fortran/adv_diff_consdiff2d.o \
adv_diff/libIBAMR_adv_diff2d_a-AdvDiffHierarchyIntegrator.o \
adv_diff/libIBAMR_adv_diff2d_a-AdvDiffHypPatchOps.o \
advect/fortran/advect_diff2d.o \
advect/fortran/advect_stable2d.o \
advect/fortran/advect_detect2d.o \
advect/fortran/godunov_predictors2d.o \
advect/fortran/weno_fluxes2d.o \
advect/libIBAMR_advect2d_a-AdvectHypPatchOps.o \
advect/libIBAMR_advect2d_a-GodunovAdvector.o \
navier_stokes/fortran/navier_stokes_bdryop2d.o \
navier_stokes/fortran/navier_stokes_divsource2d.o \
navier_stokes/fortran/navier_stokes_stabledt2d.o \
navier_stokes/fortran/navier_stokes_staggered_helpers2d.o \
navier_stokes/libIBAMR_navier_stokes2d_a-HierarchyProjector.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSHierarchyIntegrator.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSIntermediateVelocityBcCoef.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSProjectionBcCoef.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredHierarchyIntegrator.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredIntermediateVelocityBcCoef.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredPPMConvectiveOperator.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredPhysicalBoundaryHelper.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredPressureBcCoef.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredProjectionBcCoef.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredProjectionPreconditioner.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredStokesOperator.o \
navier_stokes/libIBAMR_navier_stokes2d_a-INSStaggeredVelocityBcCoef.o \
utilities/libIBAMR_utilities2d_a-get_version.o
libIBAMR2d_a_SOURCES = # intentionally blank

libIBAMR3d_a_LIBADD = \
IB/libIBAMR_IB3d_a-IBAnchorPointSpec.o \
IB/libIBAMR_IB3d_a-IBAnchorPointSpecFactory.o \
IB/libIBAMR_IB3d_a-IBBeamForceGen.o \
IB/libIBAMR_IB3d_a-IBBeamForceSpec.o \
IB/libIBAMR_IB3d_a-IBBeamForceSpecFactory.o \
IB/libIBAMR_IB3d_a-IBEulerianForceSetter.o \
IB/libIBAMR_IB3d_a-IBEulerianSourceSetter.o \
IB/libIBAMR_IB3d_a-IBHDF5Initializer.o \
IB/libIBAMR_IB3d_a-IBHierarchyIntegrator.o \
IB/libIBAMR_IB3d_a-IBInstrumentPanel.o \
IB/libIBAMR_IB3d_a-IBInstrumentationSpec.o \
IB/libIBAMR_IB3d_a-IBInstrumentationSpecFactory.o \
IB/libIBAMR_IB3d_a-IBLagrangianForceStrategy.o \
IB/libIBAMR_IB3d_a-IBLagrangianForceStrategySet.o \
IB/libIBAMR_IB3d_a-IBLagrangianSourceStrategy.o \
IB/libIBAMR_IB3d_a-IBMovingTargetPointForceGen.o \
IB/libIBAMR_IB3d_a-IBMovingTargetPointForceSpecFactory.o \
IB/libIBAMR_IB3d_a-IBMovingTargetPointForceSpec.o \
IB/libIBAMR_IB3d_a-IBSpringForceGen.o \
IB/libIBAMR_IB3d_a-IBSpringForceSpecFactory.o \
IB/libIBAMR_IB3d_a-IBSpringForceSpec.o \
IB/libIBAMR_IB3d_a-IBStaggeredHierarchyIntegrator.o \
IB/libIBAMR_IB3d_a-IBStandardForceGen.o \
IB/libIBAMR_IB3d_a-IBStandardInitializer.o \
IB/libIBAMR_IB3d_a-IBTargetPointForceGen.o \
IB/libIBAMR_IB3d_a-IBTargetPointForceSpecFactory.o \
IB/libIBAMR_IB3d_a-IBTargetPointForceSpec.o \
adv_diff/fortran/adv_diff_consdiff3d.o \
adv_diff/libIBAMR_adv_diff3d_a-AdvDiffHierarchyIntegrator.o \
adv_diff/libIBAMR_adv_diff3d_a-AdvDiffHypPatchOps.o \
advect/fortran/advect_diff3d.o \
advect/fortran/advect_stable3d.o \
advect/fortran/advect_detect3d.o \
advect/fortran/godunov_predictors3d.o \
advect/libIBAMR_advect3d_a-AdvectHypPatchOps.o \
advect/libIBAMR_advect3d_a-GodunovAdvector.o \
navier_stokes/fortran/navier_stokes_bdryop3d.o \
navier_stokes/fortran/navier_stokes_divsource3d.o \
navier_stokes/fortran/navier_stokes_stabledt3d.o \
navier_stokes/fortran/navier_stokes_staggered_helpers3d.o \
navier_stokes/libIBAMR_navier_stokes3d_a-HierarchyProjector.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSHierarchyIntegrator.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSIntermediateVelocityBcCoef.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSProjectionBcCoef.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredHierarchyIntegrator.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredIntermediateVelocityBcCoef.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredPPMConvectiveOperator.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredPhysicalBoundaryHelper.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredPressureBcCoef.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredProjectionBcCoef.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredProjectionPreconditioner.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredStokesOperator.o \
navier_stokes/libIBAMR_navier_stokes3d_a-INSStaggeredVelocityBcCoef.o \
utilities/libIBAMR_utilities3d_a-get_version.o
libIBAMR3d_a_SOURCES = # intentionally blank
