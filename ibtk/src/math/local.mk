srcs-ibtk.cpp += $(call thisdir, \
	HierarchyMathOps.cpp \
	PatchMathOps.cpp \
	PETScMatUtilities.cpp \
	PETScVecUtilities.cpp \
	PoissonUtilities.cpp \
	)

include $(call incsubdirs,fortran)
