## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

lib_LIBRARIES = libIBAMR.a
if SAMRAI2D_ENABLED
  lib_LIBRARIES += libIBAMR2d.a
endif
if SAMRAI3D_ENABLED
  lib_LIBRARIES += libIBAMR3d.a
endif

all-local: $(lib_LIBRARIES)
	for file in "lib*.a" ; do \
	  ln -fs $(PWD)/$$file $(PWD)/../src; \
	done

clean-local:
	rm -f $(PWD)/../src/libIBAMR.a
	rm -f $(PWD)/../src/libIBAMR2d.a
	rm -f $(PWD)/../src/libIBAMR3d.a

## Dimension-independent library
libIBAMR_a_SOURCES = \
$(top_builddir)/src/advect/fortran/advect_helpers.f \
../src/navier_stokes/StokesSpecifications.cpp \
../src/utilities/RNG.cpp

if LIBMESH_ENABLED
  libIBAMR_a_SOURCES +=
endif

pkg_include_HEADERS = \
../include/ibamr/RNG.h \
../include/ibamr/StokesSpecifications.h \
../include/ibamr/app_namespaces.h \
../include/ibamr/ibamr_enums.h \
../include/ibamr/ibamr_utilities.h \
../include/ibamr/namespaces.h

if LIBMESH_ENABLED
  pkg_include_HEADERS +=
endif

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
../src/IB/CIBMethod.cpp \
../src/IB/CIBMobilitySolver.cpp \
../src/IB/CIBSaddlePointSolver.cpp \
../src/IB/CIBStaggeredStokesOperator.cpp \
../src/IB/CIBStaggeredStokesSolver.cpp \
../src/IB/CIBStrategy.cpp \
../src/IB/ConstraintIBKinematics.cpp \
../src/IB/ConstraintIBMethod.cpp \
../src/IB/DirectMobilitySolver.cpp \
../src/IB/GeneralizedIBMethod.cpp \
../src/IB/IBAnchorPointSpec.cpp \
../src/IB/IBAnchorPointSpecFactory.cpp \
../src/IB/IBBeamForceSpec.cpp \
../src/IB/IBBeamForceSpecFactory.cpp \
../src/IB/IBEulerianForceFunction.cpp \
../src/IB/IBEulerianSourceFunction.cpp \
../src/IB/IBExplicitHierarchyIntegrator.cpp \
../src/IB/IBHierarchyIntegrator.cpp \
../src/IB/IBHydrodynamicForceEvaluator.cpp \
../src/IB/IBImplicitStaggeredHierarchyIntegrator.cpp \
../src/IB/IBImplicitStrategy.cpp \
../src/IB/IBInstrumentPanel.cpp \
../src/IB/IBInstrumentationSpec.cpp \
../src/IB/IBInstrumentationSpecFactory.cpp \
../src/IB/IBKirchhoffRodForceGen.cpp \
../src/IB/IBLagrangianForceStrategy.cpp \
../src/IB/IBLagrangianForceStrategySet.cpp \
../src/IB/IBLagrangianSourceStrategy.cpp \
../src/IB/IBMethod.cpp \
../src/IB/IBMethodPostProcessStrategy.cpp \
../src/IB/IBRodForceSpec.cpp \
../src/IB/IBRodForceSpecFactory.cpp \
../src/IB/IBSourceSpec.cpp \
../src/IB/IBSourceSpecFactory.cpp \
../src/IB/IBSpringForceSpec.cpp \
../src/IB/IBSpringForceSpecFactory.cpp \
../src/IB/IBStandardForceGen.cpp \
../src/IB/IBStandardSourceGen.cpp \
../src/IB/IBStandardInitializer.cpp \
../src/IB/IBStrategy.cpp \
../src/IB/IBStrategySet.cpp \
../src/IB/IBTargetPointForceSpec.cpp \
../src/IB/IBTargetPointForceSpecFactory.cpp \
../src/IB/KrylovFreeBodyMobilitySolver.cpp \
../src/IB/KrylovMobilitySolver.cpp \
../src/IB/MobilityFunctions.cpp \
../src/IB/PenaltyIBMethod.cpp \
../src/IB/StaggeredStokesIBLevelRelaxationFACOperator.cpp \
../src/IB/NonbondedForceEvaluator.cpp \
../src/IB/Wall.cpp \
../src/IB/WallForceEvaluator.cpp \
../src/adv_diff/AdvDiffCenteredConvectiveOperator.cpp \
../src/adv_diff/AdvDiffConvectiveOperatorManager.cpp \
../src/adv_diff/AdvDiffPredictorCorrectorHierarchyIntegrator.cpp \
../src/adv_diff/AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp \
../src/adv_diff/AdvDiffHierarchyIntegrator.cpp \
../src/adv_diff/AdvDiffPPMConvectiveOperator.cpp \
../src/adv_diff/AdvDiffPhysicalBoundaryUtilities.cpp \
../src/adv_diff/AdvDiffSemiImplicitHierarchyIntegrator.cpp \
../src/adv_diff/AdvDiffStochasticForcing.cpp \
../src/advect/AdvectorExplicitPredictorPatchOps.cpp \
../src/advect/AdvectorPredictorCorrectorHyperbolicPatchOps.cpp \
../src/navier_stokes/INSCollocatedCenteredConvectiveOperator.cpp \
../src/navier_stokes/INSCollocatedConvectiveOperatorManager.cpp \
../src/navier_stokes/INSCollocatedHierarchyIntegrator.cpp \
../src/navier_stokes/INSCollocatedPPMConvectiveOperator.cpp \
../src/navier_stokes/INSCollocatedVelocityBcCoef.cpp \
../src/navier_stokes/INSHierarchyIntegrator.cpp \
../src/navier_stokes/INSIntermediateVelocityBcCoef.cpp \
../src/navier_stokes/INSProjectionBcCoef.cpp \
../src/navier_stokes/INSStaggeredCenteredConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredConvectiveOperatorManager.cpp \
../src/navier_stokes/INSStaggeredHierarchyIntegrator.cpp \
../src/navier_stokes/INSStaggeredPPMConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredPressureBcCoef.cpp \
../src/navier_stokes/INSStaggeredStabilizedPPMConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredStochasticForcing.cpp \
../src/navier_stokes/INSStaggeredUpwindConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredVelocityBcCoef.cpp \
../src/navier_stokes/KrylovLinearSolverStaggeredStokesSolverInterface.cpp \
../src/navier_stokes/PETScKrylovStaggeredStokesSolver.cpp \
../src/navier_stokes/SpongeLayerForceFunction.cpp \
../src/navier_stokes/StaggeredStokesBlockFactorizationPreconditioner.cpp \
../src/navier_stokes/StaggeredStokesBlockPreconditioner.cpp \
../src/navier_stokes/StaggeredStokesFACPreconditioner.cpp \
../src/navier_stokes/StaggeredStokesFACPreconditionerStrategy.cpp \
../src/navier_stokes/StaggeredStokesLevelRelaxationFACOperator.cpp \
../src/navier_stokes/StaggeredStokesOpenBoundaryStabilizer.cpp \
../src/navier_stokes/StaggeredStokesOperator.cpp \
../src/navier_stokes/StaggeredStokesPETScLevelSolver.cpp \
../src/navier_stokes/StaggeredStokesPETScMatUtilities.cpp \
../src/navier_stokes/StaggeredStokesPETScVecUtilities.cpp \
../src/navier_stokes/StaggeredStokesPhysicalBoundaryHelper.cpp \
../src/navier_stokes/StaggeredStokesProjectionPreconditioner.cpp \
../src/navier_stokes/StaggeredStokesSolver.cpp \
../src/navier_stokes/StaggeredStokesSolverManager.cpp \
../src/navier_stokes/StokesBcCoefStrategy.cpp \
../src/utilities/ConvectiveOperator.cpp

if LIBMESH_ENABLED
  DIM_INDEPENDENT_SOURCES += \
  ../src/IB/CIBFEMethod.cpp \
  ../src/IB/IBFECentroidPostProcessor.cpp \
  ../src/IB/IBFEMethod.cpp \
  ../src/IB/IBFEPostProcessor.cpp \
  ../src/IB/IMPMethod.cpp \
  ../src/IB/IMPInitializer.cpp \
  ../src/IB/MaterialPointSpec.cpp \
  ../src/IB/MaterialPointSpecFactory.cpp
endif

pkg_include_HEADERS += \
../include/ibamr/AdvDiffCenteredConvectiveOperator.h \
../include/ibamr/AdvDiffConvectiveOperatorManager.h \
../include/ibamr/AdvDiffHierarchyIntegrator.h \
../include/ibamr/AdvDiffPPMConvectiveOperator.h \
../include/ibamr/AdvDiffPhysicalBoundaryUtilities.h \
../include/ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h \
../include/ibamr/AdvDiffPredictorCorrectorHyperbolicPatchOps.h \
../include/ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h \
../include/ibamr/AdvDiffStochasticForcing.h \
../include/ibamr/AdvectorExplicitPredictorPatchOps.h \
../include/ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h \
../include/ibamr/CIBMethod.h \
../include/ibamr/CIBMobilitySolver.h \
../include/ibamr/CIBSaddlePointSolver.h \
../include/ibamr/CIBStaggeredStokesOperator.h \
../include/ibamr/CIBStaggeredStokesSolver.h \
../include/ibamr/CIBStrategy.h \
../include/ibamr/ConstraintIBKinematics.h \
../include/ibamr/ConstraintIBMethod.h \
../include/ibamr/ConvectiveOperator.h \
../include/ibamr/GeneralizedIBMethod.h \
../include/ibamr/DirectMobilitySolver.h \
../include/ibamr/IBAnchorPointSpec.h \
../include/ibamr/IBBeamForceSpec.h \
../include/ibamr/IBExplicitHierarchyIntegrator.h \
../include/ibamr/IBHierarchyIntegrator.h \
../include/ibamr/IBHydrodynamicForceEvaluator.h \
../include/ibamr/IBImplicitStaggeredHierarchyIntegrator.h \
../include/ibamr/IBImplicitStrategy.h \
../include/ibamr/IBInstrumentPanel.h \
../include/ibamr/IBInstrumentationSpec.h \
../include/ibamr/IBKirchhoffRodForceGen.h \
../include/ibamr/IBLagrangianForceStrategy.h \
../include/ibamr/IBLagrangianForceStrategySet.h \
../include/ibamr/IBLagrangianSourceStrategy.h \
../include/ibamr/IBMethod.h \
../include/ibamr/IBMethodPostProcessStrategy.h \
../include/ibamr/IBRodForceSpec.h \
../include/ibamr/IBSourceSpec.h \
../include/ibamr/IBSpringForceFunctions.h \
../include/ibamr/IBSpringForceSpec.h \
../include/ibamr/IBStandardForceGen.h \
../include/ibamr/IBStandardInitializer.h \
../include/ibamr/IBStandardSourceGen.h \
../include/ibamr/IBStrategy.h \
../include/ibamr/IBStrategySet.h \
../include/ibamr/IBTargetPointForceSpec.h \
../include/ibamr/INSCollocatedCenteredConvectiveOperator.h \
../include/ibamr/INSCollocatedConvectiveOperatorManager.h \
../include/ibamr/INSCollocatedHierarchyIntegrator.h \
../include/ibamr/INSCollocatedPPMConvectiveOperator.h \
../include/ibamr/INSCollocatedVelocityBcCoef.h \
../include/ibamr/INSHierarchyIntegrator.h \
../include/ibamr/INSIntermediateVelocityBcCoef.h \
../include/ibamr/INSProjectionBcCoef.h \
../include/ibamr/INSStaggeredCenteredConvectiveOperator.h \
../include/ibamr/INSStaggeredConvectiveOperatorManager.h \
../include/ibamr/INSStaggeredHierarchyIntegrator.h \
../include/ibamr/INSStaggeredPPMConvectiveOperator.h \
../include/ibamr/INSStaggeredPressureBcCoef.h \
../include/ibamr/INSStaggeredStabilizedPPMConvectiveOperator.h \
../include/ibamr/INSStaggeredStochasticForcing.h \
../include/ibamr/INSStaggeredUpwindConvectiveOperator.h \
../include/ibamr/INSStaggeredVelocityBcCoef.h \
../include/ibamr/KrylovFreeBodyMobilitySolver.h \
../include/ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h \
../include/ibamr/KrylovMobilitySolver.h \
../include/ibamr/MobilityFunctions.h \
../include/ibamr/NonbondedForceEvaluator.h \
../include/ibamr/PETScKrylovStaggeredStokesSolver.h \
../include/ibamr/PenaltyIBMethod.h \
../include/ibamr/SpongeLayerForceFunction.h \
../include/ibamr/StaggeredStokesBlockFactorizationPreconditioner.h \
../include/ibamr/StaggeredStokesBlockPreconditioner.h \
../include/ibamr/StaggeredStokesFACPreconditioner.h \
../include/ibamr/StaggeredStokesFACPreconditionerStrategy.h \
../include/ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h \
../include/ibamr/StaggeredStokesLevelRelaxationFACOperator.h \
../include/ibamr/StaggeredStokesOpenBoundaryStabilizer.h \
../include/ibamr/StaggeredStokesOperator.h \
../include/ibamr/StaggeredStokesPETScLevelSolver.h \
../include/ibamr/StaggeredStokesPETScMatUtilities.h \
../include/ibamr/StaggeredStokesPETScVecUtilities.h \
../include/ibamr/StaggeredStokesPhysicalBoundaryHelper.h \
../include/ibamr/StaggeredStokesProjectionPreconditioner.h \
../include/ibamr/StaggeredStokesSolver.h \
../include/ibamr/StaggeredStokesSolverManager.h \
../include/ibamr/StokesBcCoefStrategy.h \
../include/ibamr/StokesSpecifications.h \
../include/ibamr/Wall.h \
../include/ibamr/WallForceEvaluator.h

if LIBMESH_ENABLED
  pkg_include_HEADERS += \
  ../include/ibamr/CIBFEMethod.h \
  ../include/ibamr/IBFECentroidPostProcessor.h \
  ../include/ibamr/IBFEMethod.h \
  ../include/ibamr/IBFEPatchRecoveryPostProcessor.h \
  ../include/ibamr/IBFEPostProcessor.h \
  ../include/ibamr/IMPInitializer.h \
  ../include/ibamr/IMPMethod.h \
  ../include/ibamr/MaterialPointSpec.h
endif

libIBAMR2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR2d_a_CFLAGS = $(AM_CFLAGS) -DNDIM=2
libIBAMR2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff2d.f \
$(top_builddir)/src/advect/fortran/advect_centered_derivatives2d.f \
$(top_builddir)/src/advect/fortran/advect_detect2d.f \
$(top_builddir)/src/advect/fortran/advect_diff2d.f \
$(top_builddir)/src/advect/fortran/advect_stable2d.f \
$(top_builddir)/src/advect/fortran/advect_predictors2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stochastic_forcing2d.f

libIBAMR3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR3d_a_CFLAGS = $(AM_CFLAGS) -DNDIM=3
libIBAMR3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff3d.f \
$(top_builddir)/src/advect/fortran/advect_centered_derivatives3d.f \
$(top_builddir)/src/advect/fortran/advect_detect3d.f \
$(top_builddir)/src/advect/fortran/advect_diff3d.f \
$(top_builddir)/src/advect/fortran/advect_stable3d.f \
$(top_builddir)/src/advect/fortran/advect_predictors3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stochastic_forcing3d.f
