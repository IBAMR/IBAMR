## ---------------------------------------------------------------------
##
## Copyright (c) 2009 - 2024 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

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
libIBAMR_a_SOURCES = ../src/dummy.cpp

pkg_include_HEADERS = \
../include/ibamr/app_namespaces.h \
../include/ibamr/ibamr_enums.h \
../include/ibamr/ibamr_utilities.h \
../include/ibamr/namespaces.h

if LIBMESH_ENABLED
  pkg_include_HEADERS +=
endif

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
../src/IB/BrinkmanAdvDiffBcHelper.cpp \
../src/IB/BrinkmanPenalizationRigidBodyDynamics.cpp \
../src/IB/BrinkmanPenalizationStrategy.cpp \
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
../src/IB/IBInstrumentPanel.cpp \
../src/IB/IBInstrumentationSpec.cpp \
../src/IB/IBInstrumentationSpecFactory.cpp \
../src/IB/IBInterpolantHierarchyIntegrator.cpp \
../src/IB/IBInterpolantMethod.cpp \
../src/IB/IBKirchhoffRodForceGen.cpp \
../src/IB/IBLagrangianForceStrategy.cpp \
../src/IB/IBLagrangianForceStrategySet.cpp \
../src/IB/IBLagrangianSourceStrategy.cpp \
../src/IB/IBLevelSetMethod.cpp \
../src/IB/IBMethod.cpp \
../src/IB/IBRedundantInitializer.cpp \
../src/IB/IBRodForceSpec.cpp \
../src/IB/IBRodForceSpecFactory.cpp \
../src/IB/IBHydrodynamicSurfaceForceEvaluator.cpp \
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
../src/adv_diff/AdvDiffCUIConvectiveOperator.cpp \
../src/adv_diff/AdvDiffConservativeMassScalarTransportRKIntegrator.cpp \
../src/adv_diff/AdvDiffConvectiveOperatorManager.cpp \
../src/adv_diff/AdvDiffPredictorCorrectorHierarchyIntegrator.cpp \
../src/adv_diff/AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp \
../src/adv_diff/AdvDiffHierarchyIntegrator.cpp \
../src/adv_diff/AdvDiffPPMConvectiveOperator.cpp \
../src/adv_diff/AdvDiffPhysicalBoundaryUtilities.cpp \
../src/adv_diff/AdvDiffSemiImplicitHierarchyIntegrator.cpp \
../src/adv_diff/AdvDiffStochasticForcing.cpp \
../src/adv_diff/AdvDiffWavePropConvectiveOperator.cpp \
../src/adv_diff/BrinkmanAdvDiffSemiImplicitHierarchyIntegrator.cpp \
../src/advect/AdvectorExplicitPredictorPatchOps.cpp \
../src/advect/AdvectorPredictorCorrectorHyperbolicPatchOps.cpp \
../src/complex_fluids/CFUpperConvectiveOperator.cpp \
../src/complex_fluids/CFStrategy.cpp \
../src/complex_fluids/CFGiesekusStrategy.cpp \
../src/complex_fluids/CFOldroydBStrategy.cpp \
../src/complex_fluids/CFRoliePolyStrategy.cpp \
../src/complex_fluids/CFINSForcing.cpp \
../src/level_set/FastSweepingLSMethod.cpp \
../src/level_set/LevelSetUtilities.cpp \
../src/level_set/LSInitStrategy.cpp \
../src/level_set/RelaxationLSBcCoefs.cpp \
../src/level_set/RelaxationLSMethod.cpp \
../src/navier_stokes/INSCollocatedCenteredConvectiveOperator.cpp \
../src/navier_stokes/INSCollocatedConvectiveOperatorManager.cpp \
../src/navier_stokes/INSCollocatedHierarchyIntegrator.cpp \
../src/navier_stokes/INSCollocatedPPMConvectiveOperator.cpp \
../src/navier_stokes/INSCollocatedWavePropConvectiveOperator.cpp \
../src/navier_stokes/INSCollocatedVelocityBcCoef.cpp \
../src/navier_stokes/INSHierarchyIntegrator.cpp \
../src/navier_stokes/INSIntermediateVelocityBcCoef.cpp \
../src/navier_stokes/INSProjectionBcCoef.cpp \
../src/navier_stokes/INSStaggeredCenteredConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredConvectiveOperatorManager.cpp \
../src/navier_stokes/INSStaggeredCUIConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredHierarchyIntegrator.cpp \
../src/navier_stokes/INSStaggeredPPMConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredPressureBcCoef.cpp \
../src/navier_stokes/INSStaggeredStabilizedPPMConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredStochasticForcing.cpp \
../src/navier_stokes/INSStaggeredUpwindConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredWavePropConvectiveOperator.cpp \
../src/navier_stokes/INSStaggeredVelocityBcCoef.cpp \
../src/navier_stokes/INSVCStaggeredConservativeHierarchyIntegrator.cpp \
../src/navier_stokes/INSVCStaggeredConservativeMassMomentumRKIntegrator.cpp \
../src/navier_stokes/INSVCStaggeredConservativeMassMomentumSSPRKIntegrator.cpp \
../src/navier_stokes/INSVCStaggeredHierarchyIntegrator.cpp \
../src/navier_stokes/INSVCStaggeredNonConservativeHierarchyIntegrator.cpp \
../src/navier_stokes/INSVCStaggeredPressureBcCoef.cpp \
../src/navier_stokes/INSVCStaggeredVelocityBcCoef.cpp \
../src/navier_stokes/MarangoniSurfaceTensionForceFunction.cpp \
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
../src/navier_stokes/StokesSpecifications.cpp \
../src/navier_stokes/SurfaceTensionForceFunction.cpp \
../src/navier_stokes/VCStaggeredStokesOperator.cpp \
../src/navier_stokes/VCStaggeredStokesProjectionPreconditioner.cpp \
../src/phase_change/AllenCahnHierarchyIntegrator.cpp \
../src/phase_change/CarmanKozenyDragForce.cpp \
../src/phase_change/EnthalpyHierarchyIntegrator.cpp \
../src/phase_change/HeavisideForcingFunction.cpp \
../src/phase_change/PhaseChangeDivUSourceFunction.cpp \
../src/phase_change/PhaseChangeHierarchyIntegrator.cpp \
../src/utilities/CellConvectiveOperator.cpp \
../src/utilities/ConvectiveOperator.cpp \
../src/wave_generation/FifthOrderStokesWaveGenerator.cpp \
../src/wave_generation/FirstOrderStokesWaveGenerator.cpp \
../src/wave_generation/IrregularWaveBcCoef.cpp \
../src/wave_generation/IrregularWaveGenerator.cpp \
../src/wave_generation/StokesFifthOrderWaveBcCoef.cpp \
../src/wave_generation/StokesFirstOrderWaveBcCoef.cpp \
../src/wave_generation/StokesSecondOrderWaveBcCoef.cpp \
../src/wave_generation/StokesWaveGeneratorStrategy.cpp \
../src/wave_generation/WaveDampingFunctions.cpp \
../src/wave_generation/WaveGenerationFunctions.cpp \
../src/utilities/STSMassFluxIntegrator.cpp \
../src/utilities/RNG.cpp \
../src/utilities/ins_utilities.cpp \
../src/utilities/vc_ins_utilities.cpp

if LIBMESH_ENABLED
  DIM_INDEPENDENT_SOURCES += \
  ../src/IB/FEMechanicsBase.cpp \
  ../src/IB/FEMechanicsExplicitIntegrator.cpp \
  ../src/IB/IBFECentroidPostProcessor.cpp \
  ../src/IB/IBFEDirectForcingKinematics.cpp \
  ../src/IB/IBFEInstrumentPanel.cpp \
  ../src/IB/IBFEMethod.cpp \
  ../src/IB/IBFEPostProcessor.cpp \
  ../src/IB/IBFESurfaceMethod.cpp \
  ../src/IB/IIMethod.cpp \
  ../src/IB/IMPInitializer.cpp \
  ../src/IB/IMPMethod.cpp \
  ../src/IB/MaterialPointSpec.cpp \
  ../src/IB/MaterialPointSpecFactory.cpp \
  ../src/level_set/FESurfaceDistanceEvaluator.cpp
endif

pkg_include_HEADERS += \
../include/ibamr/AdvDiffCenteredConvectiveOperator.h \
../include/ibamr/AdvDiffCUIConvectiveOperator.h \
../include/ibamr/AdvDiffConservativeMassScalarTransportRKIntegrator.h \
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
../include/ibamr/AllenCahnHierarchyIntegrator.h \
../include/ibamr/BrinkmanAdvDiffSemiImplicitHierarchyIntegrator.h \
../include/ibamr/BrinkmanAdvDiffBcHelper.h \
../include/ibamr/BrinkmanPenalizationStrategy.h \
../include/ibamr/BrinkmanPenalizationRigidBodyDynamics.h \
../include/ibamr/CIBMethod.h \
../include/ibamr/CIBMobilitySolver.h \
../include/ibamr/CIBSaddlePointSolver.h \
../include/ibamr/CIBStaggeredStokesOperator.h \
../include/ibamr/CIBStaggeredStokesSolver.h \
../include/ibamr/CIBStrategy.h \
../include/ibamr/CarmanKozenyDragForce.h \
../include/ibamr/CellConvectiveOperator.h \
../include/ibamr/ConstraintIBKinematics.h \
../include/ibamr/ConstraintIBMethod.h \
../include/ibamr/ConvectiveOperator.h \
../include/ibamr/DirectMobilitySolver.h \
../include/ibamr/EnthalpyHierarchyIntegrator.h \
../include/ibamr/FastSweepingLSMethod.h \
../include/ibamr/FifthOrderStokesWaveGenerator.h \
../include/ibamr/FirstOrderStokesWaveGenerator.h \
../include/ibamr/GeneralizedIBMethod.h \
../include/ibamr/HeavisideForcingFunction.h \
../include/ibamr/IBAnchorPointSpec.h \
../include/ibamr/IBBeamForceSpec.h \
../include/ibamr/IBExplicitHierarchyIntegrator.h \
../include/ibamr/IBHierarchyIntegrator.h \
../include/ibamr/IBHydrodynamicForceEvaluator.h \
../include/ibamr/IBImplicitStaggeredHierarchyIntegrator.h \
../include/ibamr/IBImplicitStrategy.h \
../include/ibamr/IBInstrumentPanel.h \
../include/ibamr/IBInstrumentationSpec.h \
../include/ibamr/IBInterpolantHierarchyIntegrator.h \
../include/ibamr/IBInterpolantMethod.h \
../include/ibamr/IBKirchhoffRodForceGen.h \
../include/ibamr/IBLagrangianForceStrategy.h \
../include/ibamr/IBLagrangianForceStrategySet.h \
../include/ibamr/IBLagrangianSourceStrategy.h \
../include/ibamr/IBLevelSetMethod.h \
../include/ibamr/IBMethod.h \
../include/ibamr/IBMethodPostProcessStrategy.h \
../include/ibamr/IBRedundantInitializer.h \
../include/ibamr/IBRodForceSpec.h \
../include/ibamr/IBHydrodynamicSurfaceForceEvaluator.h \
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
../include/ibamr/INSStaggeredCUIConvectiveOperator.h \
../include/ibamr/INSStaggeredHierarchyIntegrator.h \
../include/ibamr/INSStaggeredPPMConvectiveOperator.h \
../include/ibamr/INSStaggeredPressureBcCoef.h \
../include/ibamr/INSStaggeredStabilizedPPMConvectiveOperator.h \
../include/ibamr/INSStaggeredStochasticForcing.h \
../include/ibamr/INSStaggeredUpwindConvectiveOperator.h \
../include/ibamr/INSStaggeredVelocityBcCoef.h \
../include/ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h \
../include/ibamr/INSVCStaggeredConservativeMassMomentumRKIntegrator.h \
../include/ibamr/INSVCStaggeredConservativeMassMomentumSSPRKIntegrator.h \
../include/ibamr/INSVCStaggeredHierarchyIntegrator.h \
../include/ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h \
../include/ibamr/INSVCStaggeredPressureBcCoef.h \
../include/ibamr/INSVCStaggeredVelocityBcCoef.h \
../include/ibamr/IrregularWaveBcCoef.h \
../include/ibamr/IrregularWaveGenerator.h \
../include/ibamr/ins_utilities.h \
../include/ibamr/KrylovFreeBodyMobilitySolver.h \
../include/ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h \
../include/ibamr/KrylovMobilitySolver.h \
../include/ibamr/LevelSetUtilities.h \
../include/ibamr/LSInitStrategy.h \
../include/ibamr/MarangoniSurfaceTensionForceFunction.h \
../include/ibamr/MobilityFunctions.h \
../include/ibamr/NonbondedForceEvaluator.h \
../include/ibamr/PETScKrylovStaggeredStokesSolver.h \
../include/ibamr/PenaltyIBMethod.h \
../include/ibamr/PhaseChangeDivUSourceFunction.h \
../include/ibamr/PhaseChangeHierarchyIntegrator.h \
../include/ibamr/RelaxationLSBcCoefs.h \
../include/ibamr/RelaxationLSMethod.h \
../include/ibamr/SpongeLayerForceFunction.h \
../include/ibamr/STSMassFluxIntegrator.h \
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
../include/ibamr/StokesFifthOrderWaveBcCoef.h \
../include/ibamr/StokesFirstOrderWaveBcCoef.h \
../include/ibamr/StokesSecondOrderWaveBcCoef.h \
../include/ibamr/StokesSpecifications.h \
../include/ibamr/StokesSpecifications.h \
../include/ibamr/StokesWaveGeneratorStrategy.h \
../include/ibamr/SurfaceTensionForceFunction.h \
../include/ibamr/VCStaggeredStokesOperator.h \
../include/ibamr/VCStaggeredStokesProjectionPreconditioner.h \
../include/ibamr/vc_ins_utilities.h \
../include/ibamr/Wall.h \
../include/ibamr/WallForceEvaluator.h \
../include/ibamr/WaveDampingFunctions.h \
../include/ibamr/WaveGenerationFunctions.h \
../include/ibamr/RNG.h

if LIBMESH_ENABLED
  pkg_include_HEADERS += \
  ../include/ibamr/FEMechanicsBase.h \
  ../include/ibamr/FEMechanicsExplicitIntegrator.h \
  ../include/ibamr/FESurfaceDistanceEvaluator.h \
  ../include/ibamr/IBFECentroidPostProcessor.h \
  ../include/ibamr/IBFEDirectForcingKinematics.h \
  ../include/ibamr/IBFEInstrumentPanel.h \
  ../include/ibamr/IBFEMethod.h \
  ../include/ibamr/IBFEPostProcessor.h \
  ../include/ibamr/IBFESurfaceMethod.h \
  ../include/ibamr/IIMethod.h \
  ../include/ibamr/IMPInitializer.h \
  ../include/ibamr/IMPMethod.h \
  ../include/ibamr/MaterialPointSpec.h \
  ../include/ibamr/private/MaterialPointSpec-inl.h
endif

libIBAMR2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR2d_a_CFLAGS = $(AM_CFLAGS) -DNDIM=2
libIBAMR2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff2d.f \
$(top_builddir)/src/adv_diff/fortran/adv_diff_wp_convective_op2d.f \
$(top_builddir)/src/advect/fortran/advect_centered_derivatives2d.f \
$(top_builddir)/src/advect/fortran/advect_detect2d.f \
$(top_builddir)/src/advect/fortran/advect_diff2d.f \
$(top_builddir)/src/advect/fortran/advect_helpers2d.f \
$(top_builddir)/src/advect/fortran/advect_stable2d.f \
$(top_builddir)/src/advect/fortran/advect_predictors2d.f \
$(top_builddir)/src/complex_fluids/fortran/div_tensor2d.f \
$(top_builddir)/src/complex_fluids/fortran/log_upper_convective_op2d.f \
$(top_builddir)/src/complex_fluids/fortran/sqrt_upper_convective_op2d.f \
$(top_builddir)/src/complex_fluids/fortran/upper_convective_op2d.f \
$(top_builddir)/src/level_set/fortran/levelsetops2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stochastic_forcing2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_surface_tension_forcing2d.f

libIBAMR3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR3d_a_CFLAGS = $(AM_CFLAGS) -DNDIM=3
libIBAMR3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff3d.f \
$(top_builddir)/src/adv_diff/fortran/adv_diff_wp_convective_op3d.f \
$(top_builddir)/src/advect/fortran/advect_centered_derivatives3d.f \
$(top_builddir)/src/advect/fortran/advect_detect3d.f \
$(top_builddir)/src/advect/fortran/advect_diff3d.f \
$(top_builddir)/src/advect/fortran/advect_helpers3d.f \
$(top_builddir)/src/advect/fortran/advect_stable3d.f \
$(top_builddir)/src/advect/fortran/advect_predictors3d.f \
$(top_builddir)/src/complex_fluids/fortran/div_tensor3d.f \
$(top_builddir)/src/complex_fluids/fortran/log_upper_convective_op3d.f \
$(top_builddir)/src/complex_fluids/fortran/sqrt_upper_convective_op3d.f \
$(top_builddir)/src/complex_fluids/fortran/upper_convective_op3d.f \
$(top_builddir)/src/level_set/fortran/levelsetops3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stochastic_forcing3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_surface_tension_forcing3d.f
