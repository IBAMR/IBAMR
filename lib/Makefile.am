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
../src/navier_stokes/StokesSpecifications.h \
../src/utilities/RNG.h \
../src/utilities/ibamr_enums.h \
../src/utilities/ibamr_utilities.h \
../src/utilities/app_namespaces.h \
../src/utilities/namespaces.h

if LIBMESH_ENABLED
  pkg_include_HEADERS +=
endif

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
../src/IB/GeneralizedIBMethod.cpp \
../src/IB/IBAnchorPointSpec.cpp \
../src/IB/IBAnchorPointSpecFactory.cpp \
../src/IB/IBBeamForceSpec.cpp \
../src/IB/IBBeamForceSpecFactory.cpp \
../src/IB/IBEulerianForceFunction.cpp \
../src/IB/IBEulerianSourceFunction.cpp \
../src/IB/IBExplicitHierarchyIntegrator.cpp \
../src/IB/IBHierarchyIntegrator.cpp \
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
../src/IB/PenaltyIBMethod.cpp \
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
../src/navier_stokes/StaggeredStokesBoxRelaxationFACOperator.cpp \
../src/navier_stokes/StaggeredStokesFACPreconditioner.cpp \
../src/navier_stokes/StaggeredStokesFACPreconditionerStrategy.cpp \
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
  ../src/IB/IBFECentroidPostProcessor.cpp \
  ../src/IB/IBFEMethod.cpp \
  ../src/IB/IBFEPostProcessor.cpp \
  ../src/IB/IMPMethod.cpp \
  ../src/IB/IMPInitializer.cpp \
  ../src/IB/MaterialPointSpec.cpp \
  ../src/IB/MaterialPointSpecFactory.cpp
endif

pkg_include_HEADERS += \
../src/IB/GeneralizedIBMethod.h \
../src/IB/IBAnchorPointSpec-inl.h \
../src/IB/IBAnchorPointSpec.h \
../src/IB/IBBeamForceSpec-inl.h \
../src/IB/IBBeamForceSpec.h \
../src/IB/IBExplicitHierarchyIntegrator.h \
../src/IB/IBHierarchyIntegrator.h \
../src/IB/IBImplicitStaggeredHierarchyIntegrator.h \
../src/IB/IBImplicitStrategy.h \
../src/IB/IBInstrumentPanel.h \
../src/IB/IBInstrumentationSpec-inl.h \
../src/IB/IBInstrumentationSpec.h \
../src/IB/IBKirchhoffRodForceGen.h \
../src/IB/IBLagrangianForceStrategy.h \
../src/IB/IBLagrangianForceStrategySet.h \
../src/IB/IBLagrangianSourceStrategy.h \
../src/IB/IBMethod.h \
../src/IB/IBMethodPostProcessStrategy.h \
../src/IB/IBRodForceSpec-inl.h \
../src/IB/IBRodForceSpec.h \
../src/IB/IBSourceSpec-inl.h \
../src/IB/IBSourceSpec.h \
../src/IB/IBSpringForceFunctions.h \
../src/IB/IBSpringForceSpec-inl.h \
../src/IB/IBSpringForceSpec.h \
../src/IB/IBStandardForceGen.h \
../src/IB/IBStandardSourceGen.h \
../src/IB/IBStandardInitializer.h \
../src/IB/IBStrategy.h \
../src/IB/IBStrategySet.h \
../src/IB/IBTargetPointForceSpec-inl.h \
../src/IB/IBTargetPointForceSpec.h \
../src/IB/PenaltyIBMethod.h \
../src/adv_diff/AdvDiffCenteredConvectiveOperator.h \
../src/adv_diff/AdvDiffConvectiveOperatorManager.h \
../src/adv_diff/AdvDiffPredictorCorrectorHierarchyIntegrator.h \
../src/adv_diff/AdvDiffPredictorCorrectorHyperbolicPatchOps.h \
../src/adv_diff/AdvDiffHierarchyIntegrator.h \
../src/adv_diff/AdvDiffPPMConvectiveOperator.h \
../src/adv_diff/AdvDiffPhysicalBoundaryUtilities.h \
../src/adv_diff/AdvDiffSemiImplicitHierarchyIntegrator.h \
../src/adv_diff/AdvDiffStochasticForcing.h \
../src/advect/AdvectorExplicitPredictorPatchOps.h \
../src/advect/AdvectorPredictorCorrectorHyperbolicPatchOps.h \
../src/navier_stokes/INSCollocatedCenteredConvectiveOperator.h \
../src/navier_stokes/INSCollocatedConvectiveOperatorManager.h \
../src/navier_stokes/INSCollocatedHierarchyIntegrator.h \
../src/navier_stokes/INSCollocatedPPMConvectiveOperator.h \
../src/navier_stokes/INSCollocatedVelocityBcCoef.h \
../src/navier_stokes/INSHierarchyIntegrator.h \
../src/navier_stokes/INSIntermediateVelocityBcCoef.h \
../src/navier_stokes/INSProjectionBcCoef.h \
../src/navier_stokes/INSStaggeredCenteredConvectiveOperator.h \
../src/navier_stokes/INSStaggeredConvectiveOperatorManager.h \
../src/navier_stokes/INSStaggeredHierarchyIntegrator.h \
../src/navier_stokes/INSStaggeredPPMConvectiveOperator.h \
../src/navier_stokes/INSStaggeredPressureBcCoef.h \
../src/navier_stokes/INSStaggeredStabilizedPPMConvectiveOperator.h \
../src/navier_stokes/INSStaggeredStochasticForcing.h \
../src/navier_stokes/INSStaggeredUpwindConvectiveOperator.h \
../src/navier_stokes/INSStaggeredVelocityBcCoef.h \
../src/navier_stokes/KrylovLinearSolverStaggeredStokesSolverInterface.h \
../src/navier_stokes/PETScKrylovStaggeredStokesSolver.h \
../src/navier_stokes/SpongeLayerForceFunction.h \
../src/navier_stokes/StaggeredStokesBlockFactorizationPreconditioner.h \
../src/navier_stokes/StaggeredStokesBlockPreconditioner.h \
../src/navier_stokes/StaggeredStokesBoxRelaxationFACOperator.h \
../src/navier_stokes/StaggeredStokesFACPreconditioner.h \
../src/navier_stokes/StaggeredStokesFACPreconditionerStrategy.h \
../src/navier_stokes/StaggeredStokesOpenBoundaryStabilizer.h \
../src/navier_stokes/StaggeredStokesOperator.h \
../src/navier_stokes/StaggeredStokesPETScLevelSolver.h \
../src/navier_stokes/StaggeredStokesPETScMatUtilities.h \
../src/navier_stokes/StaggeredStokesPETScVecUtilities.h \
../src/navier_stokes/StaggeredStokesPhysicalBoundaryHelper.h \
../src/navier_stokes/StaggeredStokesProjectionPreconditioner.h \
../src/navier_stokes/StaggeredStokesSolver.h \
../src/navier_stokes/StaggeredStokesSolverManager.h \
../src/navier_stokes/StokesBcCoefStrategy.h \
../src/utilities/ConvectiveOperator.h

if LIBMESH_ENABLED
  pkg_include_HEADERS += \
  ../src/IB/IBFECentroidPostProcessor.h \
  ../src/IB/IBFEMethod.h \
  ../src/IB/IBFEPostProcessor.h \
  ../src/IB/IMPMethod.h \
  ../src/IB/IMPInitializer.h \
  ../src/IB/MaterialPointSpec.h
endif

libIBAMR2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
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
