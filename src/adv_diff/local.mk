srcs-core.cpp += $(call thisdir, \
	AdvDiffCUIConvectiveOperator.cpp \
	AdvDiffCenteredConvectiveOperator.cpp \
	AdvDiffConvectiveOperatorManager.cpp \
	AdvDiffHierarchyIntegrator.cpp \
	AdvDiffPPMConvectiveOperator.cpp \
	AdvDiffPhysicalBoundaryUtilities.cpp \
	AdvDiffPredictorCorrectorHierarchyIntegrator.cpp \
	AdvDiffPredictorCorrectorHyperbolicPatchOps.cpp \
	AdvDiffSemiImplicitHierarchyIntegrator.cpp \
	AdvDiffStochasticForcing.cpp \
	AdvDiffWavePropConvectiveOperator.cpp \
	)

include $(call incsubdirs,fortran)
