srcs-core.cpp += $(call thisdir, \
	CIBMethod.cpp \
	CIBMobilitySolver.cpp \
	CIBSaddlePointSolver.cpp \
	CIBStaggeredStokesOperator.cpp \
	CIBStaggeredStokesSolver.cpp \
	CIBStrategy.cpp \
	ConstraintIBKinematics.cpp \
	ConstraintIBMethod.cpp \
	DirectMobilitySolver.cpp \
	GeneralizedIBMethod.cpp \
	IBAnchorPointSpec.cpp \
	IBAnchorPointSpecFactory.cpp \
	IBBeamForceSpec.cpp \
	IBBeamForceSpecFactory.cpp \
	IBEulerianForceFunction.cpp \
	IBEulerianSourceFunction.cpp \
	IBExplicitHierarchyIntegrator.cpp \
	IBHierarchyIntegrator.cpp \
	IBHydrodynamicForceEvaluator.cpp \
	IBHydrodynamicSurfaceForceEvaluator.cpp \
	IBImplicitStaggeredHierarchyIntegrator.cpp \
	IBImplicitStrategy.cpp \
	IBInstrumentPanel.cpp \
	IBInstrumentationSpec.cpp \
	IBInstrumentationSpecFactory.cpp \
	IBKirchhoffRodForceGen.cpp \
	IBLagrangianForceStrategy.cpp \
	IBLagrangianForceStrategySet.cpp \
	IBLagrangianSourceStrategy.cpp \
	IBMethod.cpp \
	IBMethodPostProcessStrategy.cpp \
	IBRedundantInitializer.cpp \
	IBRodForceSpec.cpp \
	IBRodForceSpecFactory.cpp \
	IBSourceSpec.cpp \
	IBSourceSpecFactory.cpp \
	IBSpringForceSpec.cpp \
	IBSpringForceSpecFactory.cpp \
	IBStandardForceGen.cpp \
	IBStandardInitializer.cpp \
	IBStandardSourceGen.cpp \
	IBStrategy.cpp \
	IBStrategySet.cpp \
	IBTargetPointForceSpec.cpp \
	IBTargetPointForceSpecFactory.cpp \
	KrylovFreeBodyMobilitySolver.cpp \
	KrylovMobilitySolver.cpp \
	MobilityFunctions.cpp \
	NonbondedForceEvaluator.cpp \
	PenaltyIBMethod.cpp \
	StaggeredStokesIBLevelRelaxationFACOperator.cpp \
	Wall.cpp \
	WallForceEvaluator.cpp \
	)

ifneq ($(IBAMR_LIBMESH_LIB),)

srcs-core.cpp += $(call thisdir, \
	CIBFEMethod.cpp \
	IBFECentroidPostProcessor.cpp \
	IBFEInstrumentPanel.cpp \
	IBFEMethod.cpp \
	IBFEPatchRecoveryPostProcessor.cpp \
	IBFEPostProcessor.cpp \
	IBFESurfaceMethod.cpp \
	IMPInitializer.cpp \
	IMPMethod.cpp \
	MaterialPointSpec.cpp \
	MaterialPointSpecFactory.cpp \
	)

endif
