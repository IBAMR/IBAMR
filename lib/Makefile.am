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
$(top_builddir)/src/fortran/sqrtm.f \
$(top_srcdir)/src/navier_stokes/StokesSpecifications.cpp \
$(top_srcdir)/src/utilities/RNG.cpp

if LIBMESH_ENABLED
libIBAMR_a_SOURCES += \
$(top_srcdir)/src/solid_mechanics/ExplicitFEMechanicsSolver.cpp
endif

pkg_include_HEADERS = \
$(top_srcdir)/src/navier_stokes/StokesSpecifications.h \
$(top_srcdir)/src/utilities/RNG.h \
$(top_srcdir)/src/utilities/ibamr_enums.h \
$(top_srcdir)/src/utilities/ibamr_utilities.h \
$(top_srcdir)/src/utilities/app_namespaces.h \
$(top_srcdir)/src/utilities/namespaces.h

if LIBMESH_ENABLED
pkg_include_HEADERS += \
$(top_srcdir)/src/solid_mechanics/ExplicitFEMechanicsSolver.h
endif

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
$(top_srcdir)/src/IB/GeneralizedIBMethod.cpp \
$(top_srcdir)/src/IB/IBAnchorPointSpec.cpp \
$(top_srcdir)/src/IB/IBAnchorPointSpecFactory.cpp \
$(top_srcdir)/src/IB/IBBeamForceSpec.cpp \
$(top_srcdir)/src/IB/IBBeamForceSpecFactory.cpp \
$(top_srcdir)/src/IB/IBEulerianForceFunction.cpp \
$(top_srcdir)/src/IB/IBEulerianSourceFunction.cpp \
$(top_srcdir)/src/IB/IBExplicitHierarchyIntegrator.cpp \
$(top_srcdir)/src/IB/IBHierarchyIntegrator.cpp \
$(top_srcdir)/src/IB/IBImplicitStaggeredHierarchyIntegrator.cpp \
$(top_srcdir)/src/IB/IBImplicitStaggeredPETScLevelSolver.cpp \
$(top_srcdir)/src/IB/IBInstrumentPanel.cpp \
$(top_srcdir)/src/IB/IBInstrumentationSpec.cpp \
$(top_srcdir)/src/IB/IBInstrumentationSpecFactory.cpp \
$(top_srcdir)/src/IB/IBKirchhoffRodForceGen.cpp \
$(top_srcdir)/src/IB/IBLagrangianForceStrategy.cpp \
$(top_srcdir)/src/IB/IBLagrangianForceStrategySet.cpp \
$(top_srcdir)/src/IB/IBLagrangianSourceStrategy.cpp \
$(top_srcdir)/src/IB/IBMethod.cpp \
$(top_srcdir)/src/IB/IBMethodPostProcessStrategy.cpp \
$(top_srcdir)/src/IB/IBRodForceSpec.cpp \
$(top_srcdir)/src/IB/IBRodForceSpecFactory.cpp \
$(top_srcdir)/src/IB/IBSourceSpec.cpp \
$(top_srcdir)/src/IB/IBSourceSpecFactory.cpp \
$(top_srcdir)/src/IB/IBSpringForceSpec.cpp \
$(top_srcdir)/src/IB/IBSpringForceSpecFactory.cpp \
$(top_srcdir)/src/IB/IBStandardForceGen.cpp \
$(top_srcdir)/src/IB/IBStandardSourceGen.cpp \
$(top_srcdir)/src/IB/IBStandardInitializer.cpp \
$(top_srcdir)/src/IB/IBStrategy.cpp \
$(top_srcdir)/src/IB/IBStrategySet.cpp \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.cpp \
$(top_srcdir)/src/IB/IBTargetPointForceSpecFactory.cpp \
$(top_srcdir)/src/IB/PenaltyIBMethod.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffCenteredConvectiveOperator.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffConvectiveOperatorManager.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffPredictorCorrectorHierarchyIntegrator.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffHierarchyIntegrator.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffPPMConvectiveOperator.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffPhysicalBoundaryUtilities.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffSemiImplicitHierarchyIntegrator.cpp \
$(top_srcdir)/src/adv_diff/AdvDiffStochasticForcing.cpp \
$(top_srcdir)/src/advect/AdvectorPredictorCorrectorHyperbolicPatchOps.cpp \
$(top_srcdir)/src/advect/AdvectorExplicitPredictorStrategy.cpp \
$(top_srcdir)/src/navier_stokes/INSCollocatedCenteredConvectiveOperator.cpp \
$(top_srcdir)/src/navier_stokes/INSCollocatedConvectiveOperatorManager.cpp \
$(top_srcdir)/src/navier_stokes/INSCollocatedHierarchyIntegrator.cpp \
$(top_srcdir)/src/navier_stokes/INSCollocatedPPMConvectiveOperator.cpp \
$(top_srcdir)/src/navier_stokes/INSHierarchyIntegrator.cpp \
$(top_srcdir)/src/navier_stokes/INSIntermediateVelocityBcCoef.cpp \
$(top_srcdir)/src/navier_stokes/INSProjectionBcCoef.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredCenteredConvectiveOperator.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredConvectiveOperatorManager.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredHierarchyIntegrator.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredPPMConvectiveOperator.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredPressureBcCoef.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredStabilizedPPMConvectiveOperator.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredStochasticForcing.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredUpwindConvectiveOperator.cpp \
$(top_srcdir)/src/navier_stokes/INSStaggeredVelocityBcCoef.cpp \
$(top_srcdir)/src/navier_stokes/KrylovLinearSolverStaggeredStokesSolverInterface.cpp \
$(top_srcdir)/src/navier_stokes/PETScKrylovStaggeredStokesSolver.cpp \
$(top_srcdir)/src/navier_stokes/SpongeLayerForceFunction.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesBlockFactorizationPreconditioner.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesBlockPreconditioner.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesBoxRelaxationFACOperator.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesFACPreconditioner.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesFACPreconditionerStrategy.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesOpenBoundaryStabilizer.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesOperator.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPETScLevelSolver.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPETScMatUtilities.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPETScVecUtilities.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPhysicalBoundaryHelper.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesProjectionPreconditioner.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesSolver.cpp \
$(top_srcdir)/src/navier_stokes/StaggeredStokesSolverManager.cpp \
$(top_srcdir)/src/navier_stokes/StokesBcCoefStrategy.cpp \
$(top_srcdir)/src/utilities/ConvectiveOperator.cpp

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
$(top_srcdir)/src/IB/IBFEMethod.cpp \
$(top_srcdir)/src/IB/IMPMethod.cpp \
$(top_srcdir)/src/IB/IMPInitializer.cpp \
$(top_srcdir)/src/IB/MaterialPointSpec.cpp \
$(top_srcdir)/src/IB/MaterialPointSpecFactory.cpp \
$(top_srcdir)/src/IB/SimplifiedIBFEMethod.cpp
endif

pkg_include_HEADERS += \
$(top_srcdir)/src/IB/GeneralizedIBMethod.h \
$(top_srcdir)/src/IB/IBAnchorPointSpec-inl.h \
$(top_srcdir)/src/IB/IBAnchorPointSpec.h \
$(top_srcdir)/src/IB/IBBeamForceSpec-inl.h \
$(top_srcdir)/src/IB/IBBeamForceSpec.h \
$(top_srcdir)/src/IB/IBExplicitHierarchyIntegrator.h \
$(top_srcdir)/src/IB/IBHierarchyIntegrator.h \
$(top_srcdir)/src/IB/IBImplicitStaggeredHierarchyIntegrator.h \
$(top_srcdir)/src/IB/IBImplicitStaggeredPETScLevelSolver.h \
$(top_srcdir)/src/IB/IBInstrumentPanel.h \
$(top_srcdir)/src/IB/IBInstrumentationSpec-inl.h \
$(top_srcdir)/src/IB/IBInstrumentationSpec.h \
$(top_srcdir)/src/IB/IBKirchhoffRodForceGen.h \
$(top_srcdir)/src/IB/IBLagrangianForceStrategy.h \
$(top_srcdir)/src/IB/IBLagrangianForceStrategySet.h \
$(top_srcdir)/src/IB/IBLagrangianSourceStrategy.h \
$(top_srcdir)/src/IB/IBMethod.h \
$(top_srcdir)/src/IB/IBMethodPostProcessStrategy.h \
$(top_srcdir)/src/IB/IBRodForceSpec-inl.h \
$(top_srcdir)/src/IB/IBRodForceSpec.h \
$(top_srcdir)/src/IB/IBSourceSpec-inl.h \
$(top_srcdir)/src/IB/IBSourceSpec.h \
$(top_srcdir)/src/IB/IBSpringForceFunctions.h \
$(top_srcdir)/src/IB/IBSpringForceSpec-inl.h \
$(top_srcdir)/src/IB/IBSpringForceSpec.h \
$(top_srcdir)/src/IB/IBStandardForceGen.h \
$(top_srcdir)/src/IB/IBStandardSourceGen.h \
$(top_srcdir)/src/IB/IBStandardInitializer.h \
$(top_srcdir)/src/IB/IBStrategy.h \
$(top_srcdir)/src/IB/IBStrategySet.h \
$(top_srcdir)/src/IB/IBTargetPointForceSpec-inl.h \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.h \
$(top_srcdir)/src/IB/PenaltyIBMethod.h \
$(top_srcdir)/src/adv_diff/AdvDiffCenteredConvectiveOperator.h \
$(top_srcdir)/src/adv_diff/AdvDiffConvectiveOperatorManager.h \
$(top_srcdir)/src/adv_diff/AdvDiffPredictorCorrectorHierarchyIntegrator.h \
$(top_srcdir)/src/adv_diff/AdvDiffPredictorCorrectorHyperbolicPatchOps.h \
$(top_srcdir)/src/adv_diff/AdvDiffHierarchyIntegrator.h \
$(top_srcdir)/src/adv_diff/AdvDiffPPMConvectiveOperator.h \
$(top_srcdir)/src/adv_diff/AdvDiffPhysicalBoundaryUtilities.h \
$(top_srcdir)/src/adv_diff/AdvDiffSemiImplicitHierarchyIntegrator.h \
$(top_srcdir)/src/adv_diff/AdvDiffStochasticForcing.h \
$(top_srcdir)/src/advect/AdvectorPredictorCorrectorHyperbolicPatchOps.h \
$(top_srcdir)/src/advect/AdvectorExplicitPredictorStrategy.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedCenteredConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedConvectiveOperatorManager.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedPPMConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSIntermediateVelocityBcCoef.h \
$(top_srcdir)/src/navier_stokes/INSProjectionBcCoef.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredCenteredConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredConvectiveOperatorManager.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredPPMConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredPressureBcCoef.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredStabilizedPPMConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredStochasticForcing.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredUpwindConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredVelocityBcCoef.h \
$(top_srcdir)/src/navier_stokes/KrylovLinearSolverStaggeredStokesSolverInterface.h \
$(top_srcdir)/src/navier_stokes/PETScKrylovStaggeredStokesSolver.h \
$(top_srcdir)/src/navier_stokes/SpongeLayerForceFunction.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesBlockFactorizationPreconditioner.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesBlockPreconditioner.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesBoxRelaxationFACOperator.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesFACPreconditioner.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesFACPreconditionerStrategy.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesOpenBoundaryStabilizer.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesOperator.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPETScLevelSolver.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPETScMatUtilities.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPETScVecUtilities.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesPhysicalBoundaryHelper.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesProjectionPreconditioner.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesSolver.h \
$(top_srcdir)/src/navier_stokes/StaggeredStokesSolverManager.h \
$(top_srcdir)/src/navier_stokes/StokesBcCoefStrategy.h \
$(top_srcdir)/src/utilities/ConvectiveOperator.h

if LIBMESH_ENABLED
pkg_include_HEADERS += \
$(top_srcdir)/src/IB/IBFEMethod.h \
$(top_srcdir)/src/IB/IMPMethod.h \
$(top_srcdir)/src/IB/IMPInitializer.h \
$(top_srcdir)/src/IB/MaterialPointSpec.h \
$(top_srcdir)/src/IB/SimplifiedIBFEMethod.cpp
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
