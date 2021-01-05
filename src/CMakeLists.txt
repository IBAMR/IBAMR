## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

ADD_LIBRARY(IBAMR2d)
ADD_LIBRARY(IBAMR3d)

# 'base' names (without, 2d.f.m4 or 3d.f.m4) for all the fortran that is
# explicitly compiled
SET(FORTRAN_SRC_BASE
  adv_diff/fortran/adv_diff_consdiff
  adv_diff/fortran/adv_diff_wp_convective_op
  advect/fortran/advect_centered_derivatives
  advect/fortran/advect_detect
  advect/fortran/advect_diff
  advect/fortran/advect_helpers
  advect/fortran/advect_predictors
  advect/fortran/advect_stable
  complex_fluids/fortran/div_tensor
  complex_fluids/fortran/log_upper_convective_op
  complex_fluids/fortran/sqrt_upper_convective_op
  complex_fluids/fortran/upper_convective_op
  level_set/fortran/levelsetops
  navier_stokes/fortran/navier_stokes_bdryop
  navier_stokes/fortran/navier_stokes_divsource
  navier_stokes/fortran/navier_stokes_stabledt
  navier_stokes/fortran/navier_stokes_staggered_derivatives
  navier_stokes/fortran/navier_stokes_staggered_helpers
  navier_stokes/fortran/navier_stokes_stochastic_forcing
  navier_stokes/fortran/navier_stokes_surface_tension_forcing
  )

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    IBAMR_PROCESS_M4(NAME ${_fortran_root_name}${_d}d.f)
  ENDFOREACH()

  SET("FORTRAN_GENERATED_SRC${_d}D" "")
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    LIST(APPEND "FORTRAN_GENERATED_SRC${_d}D" "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_root_name}${_d}d.f")
  ENDFOREACH()
ENDFOREACH()

#
# set up C++:
#
SET(CXX_SRC
  # adv diff
  adv_diff/AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp
  adv_diff/AdvDiffConvectiveOperatorManager.cpp
  adv_diff/AdvDiffCUIConvectiveOperator.cpp
  adv_diff/AdvDiffCenteredConvectiveOperator.cpp
  adv_diff/AdvDiffPredictorCorrectorHierarchyIntegrator.cpp
  adv_diff/AdvDiffSemiImplicitHierarchyIntegrator.cpp
  adv_diff/AdvDiffHierarchyIntegrator.cpp
  adv_diff/AdvDiffPhysicalBoundaryUtilities.cpp
  adv_diff/AdvDiffWavePropConvectiveOperator.cpp
  adv_diff/AdvDiffStochasticForcing.cpp
  adv_diff/AdvDiffPPMConvectiveOperator.cpp
  adv_diff/BrinkmanAdvDiffSemiImplicitHierarchyIntegrator.cpp

  # wave generation
  wave_generation/StokesFirstOrderWaveBcCoef.cpp
  wave_generation/StokesSecondOrderWaveBcCoef.cpp
  wave_generation/StokesFifthOrderWaveBcCoef.cpp
  wave_generation/StokesWaveGeneratorStrategy.cpp
  wave_generation/WaveGenerationFunctions.cpp
  wave_generation/WaveDampingFunctions.cpp
  wave_generation/IrregularWaveGenerator.cpp
  wave_generation/FifthOrderStokesWaveGenerator.cpp
  wave_generation/FirstOrderStokesWaveGenerator.cpp
  wave_generation/IrregularWaveBcCoef.cpp

  # navier stokes
  navier_stokes/INSIntermediateVelocityBcCoef.cpp
  navier_stokes/INSStaggeredCenteredConvectiveOperator.cpp
  navier_stokes/StaggeredStokesPETScMatUtilities.cpp
  navier_stokes/PETScKrylovStaggeredStokesSolver.cpp
  navier_stokes/StokesSpecifications.cpp
  navier_stokes/INSStaggeredStochasticForcing.cpp
  navier_stokes/VCStaggeredStokesProjectionPreconditioner.cpp
  navier_stokes/StaggeredStokesSolverManager.cpp
  navier_stokes/INSStaggeredUpwindConvectiveOperator.cpp
  navier_stokes/INSCollocatedHierarchyIntegrator.cpp
  navier_stokes/VCStaggeredStokesOperator.cpp
  navier_stokes/StaggeredStokesPETScVecUtilities.cpp
  navier_stokes/INSCollocatedVelocityBcCoef.cpp
  navier_stokes/INSVCStaggeredConservativeHierarchyIntegrator.cpp
  navier_stokes/StaggeredStokesPhysicalBoundaryHelper.cpp
  navier_stokes/INSStaggeredVelocityBcCoef.cpp
  navier_stokes/INSStaggeredCUIConvectiveOperator.cpp
  navier_stokes/INSCollocatedPPMConvectiveOperator.cpp
  navier_stokes/INSStaggeredWavePropConvectiveOperator.cpp
  navier_stokes/INSCollocatedConvectiveOperatorManager.cpp
  navier_stokes/INSStaggeredHierarchyIntegrator.cpp
  navier_stokes/INSProjectionBcCoef.cpp
  navier_stokes/INSHierarchyIntegrator.cpp
  navier_stokes/StaggeredStokesPETScLevelSolver.cpp
  navier_stokes/StaggeredStokesBlockFactorizationPreconditioner.cpp
  navier_stokes/StaggeredStokesBoxRelaxationFACOperator.cpp
  navier_stokes/INSStaggeredConvectiveOperatorManager.cpp
  navier_stokes/INSStaggeredStabilizedPPMConvectiveOperator.cpp
  navier_stokes/INSCollocatedWavePropConvectiveOperator.cpp
  navier_stokes/SpongeLayerForceFunction.cpp
  navier_stokes/INSVCStaggeredVelocityBcCoef.cpp
  navier_stokes/KrylovLinearSolverStaggeredStokesSolverInterface.cpp
  navier_stokes/StaggeredStokesFACPreconditioner.cpp
  navier_stokes/INSStaggeredPPMConvectiveOperator.cpp
  navier_stokes/INSStaggeredPressureBcCoef.cpp
  navier_stokes/StokesBcCoefStrategy.cpp
  navier_stokes/INSVCStaggeredNonConservativeHierarchyIntegrator.cpp
  navier_stokes/SurfaceTensionForceFunction.cpp
  navier_stokes/StaggeredStokesLevelRelaxationFACOperator.cpp
  navier_stokes/INSVCStaggeredPressureBcCoef.cpp
  navier_stokes/StaggeredStokesProjectionPreconditioner.cpp
  navier_stokes/StaggeredStokesFACPreconditionerStrategy.cpp
  navier_stokes/INSVCStaggeredConservativeMassMomentumIntegrator.cpp
  navier_stokes/StaggeredStokesOpenBoundaryStabilizer.cpp
  navier_stokes/StaggeredStokesBlockPreconditioner.cpp
  navier_stokes/StaggeredStokesOperator.cpp
  navier_stokes/INSCollocatedCenteredConvectiveOperator.cpp
  navier_stokes/StaggeredStokesSolver.cpp
  navier_stokes/INSVCStaggeredHierarchyIntegrator.cpp

  # IB
  IB/IBInstrumentPanel.cpp
  IB/BrinkmanPenalizationStrategy.cpp
  IB/BrinkmanAdvDiffBcHelper.cpp
  IB/PenaltyIBMethod.cpp
  IB/GeneralizedIBMethod.cpp
  IB/IBLagrangianForceStrategy.cpp
  IB/KrylovMobilitySolver.cpp
  IB/IBHydrodynamicForceEvaluator.cpp
  IB/DirectMobilitySolver.cpp
  IB/IBHydrodynamicSurfaceForceEvaluator.cpp
  IB/IBRodForceSpecFactory.cpp
  IB/IBExplicitHierarchyIntegrator.cpp
  IB/CIBStaggeredStokesSolver.cpp
  IB/IBKirchhoffRodForceGen.cpp
  IB/IBAnchorPointSpec.cpp
  IB/IBImplicitStaggeredHierarchyIntegrator.cpp
  IB/IBFECentroidPostProcessor.cpp
  IB/IBRedundantInitializer.cpp
  IB/IBAnchorPointSpecFactory.cpp
  IB/CIBStrategy.cpp
  IB/IBInterpolantHierarchyIntegrator.cpp
  IB/IBFESurfaceMethod.cpp
  IB/IBInstrumentationSpecFactory.cpp
  # This file is not yet finished
  # IB/IBFEPatchRecoveryPostProcessor.cpp
  IB/IBFEPostProcessor.cpp
  IB/BrinkmanPenalizationRigidBodyDynamics.cpp
  IB/CIBMobilitySolver.cpp
  IB/IBStrategy.cpp
  IB/IBTargetPointForceSpec.cpp
  IB/IBStrategySet.cpp
  IB/IBHierarchyIntegrator.cpp
  IB/MaterialPointSpecFactory.cpp
  IB/IBInstrumentationSpec.cpp
  IB/IBFEMethod.cpp
  IB/ConstraintIBMethod.cpp
  IB/IBEulerianSourceFunction.cpp
  IB/IBBeamForceSpecFactory.cpp
  IB/MobilityFunctions.cpp
  IB/IMPMethod.cpp
  IB/IMPInitializer.cpp
  IB/IBStandardForceGen.cpp
  IB/CIBStaggeredStokesOperator.cpp
  IB/Wall.cpp
  IB/IBStandardInitializer.cpp
  IB/IBInterpolantMethod.cpp
  IB/IBEulerianForceFunction.cpp
  IB/NonbondedForceEvaluator.cpp
  IB/IBFEInstrumentPanel.cpp
  IB/StaggeredStokesIBLevelRelaxationFACOperator.cpp
  IB/MaterialPointSpec.cpp
  IB/KrylovFreeBodyMobilitySolver.cpp
  IB/IBSpringForceSpecFactory.cpp
  IB/ConstraintIBKinematics.cpp
  IB/IBLevelSetMethod.cpp
  IB/IBLagrangianForceStrategySet.cpp
  IB/IBRodForceSpec.cpp
  IB/FEMechanicsBase.cpp
  IB/FEMechanicsExplicitIntegrator.cpp
  IB/IBSourceSpec.cpp
  IB/IBSourceSpecFactory.cpp
  IB/CIBSaddlePointSolver.cpp
  IB/WallForceEvaluator.cpp
  IB/IBTargetPointForceSpecFactory.cpp
  IB/IBMethod.cpp
  IB/IBFEDirectForcingKinematics.cpp
  IB/IBBeamForceSpec.cpp
  IB/IBSpringForceSpec.cpp
  IB/CIBMethod.cpp
  IB/IBLagrangianSourceStrategy.cpp
  IB/IBStandardSourceGen.cpp

  # complex fluids
  complex_fluids/CFRoliePolyRelaxation.cpp
  complex_fluids/CFGiesekusRelaxation.cpp
  complex_fluids/CFRelaxationOperator.cpp
  complex_fluids/CFINSForcing.cpp
  complex_fluids/CFOldroydBRelaxation.cpp
  complex_fluids/CFUpperConvectiveOperator.cpp

  # advect
  advect/AdvectorExplicitPredictorPatchOps.cpp
  advect/AdvectorPredictorCorrectorHyperbolicPatchOps.cpp

  # level set
  level_set/FESurfaceDistanceEvaluator.cpp
  level_set/LSInitStrategy.cpp
  level_set/FastSweepingLSMethod.cpp
  level_set/RelaxationLSBcCoefs.cpp
  level_set/RelaxationLSMethod.cpp
  utilities/RNG.cpp
  utilities/ConvectiveOperator.cpp
  )

TARGET_SOURCES(IBAMR2d PRIVATE ${FORTRAN_GENERATED_SRC2D} ${CXX_SRC})
TARGET_SOURCES(IBAMR3d PRIVATE ${FORTRAN_GENERATED_SRC3D} ${CXX_SRC})

TARGET_LINK_LIBRARIES(IBAMR2d PUBLIC IBAMRHeaders IBTKHeaders)
TARGET_LINK_LIBRARIES(IBAMR3d PUBLIC IBAMRHeaders IBTKHeaders)

# Since libIBAMR and libIBTK have the same dependencies (and libIBAMR always
# depends on libIBTK), satisfy the dependencies (and NDIM definition) for
# libIBAMR by linking against libIBTK.
#
# The order we resolve libIBAMR's dependencies in is important since SAMRAI is
# usually staticly linked - if we do SAMRAI then IBTK we will end up with
# duplicate symbols, so it has to be the other way around (so that we pick up
# any SAMRAI symbols in libIBTK rather than adding them to libIBAMR).
TARGET_LINK_LIBRARIES(IBAMR2d PUBLIC IBTK2d)
TARGET_LINK_LIBRARIES(IBAMR3d PUBLIC IBTK3d)

INSTALL(TARGETS IBAMR2d EXPORT IBAMRTargets COMPONENT library)
INSTALL(TARGETS IBAMR3d EXPORT IBAMRTargets COMPONENT library)
