srcs-core.cpp += $(call thisdir, \
	AdvectorExplicitPredictorPatchOps.cpp \
	AdvectorPredictorCorrectorHyperbolicPatchOps.cpp \
	)

include $(call incsubdirs,fortran)
