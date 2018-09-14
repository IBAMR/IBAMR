srcs-ibtk.cpp += $(call thisdir, \
	LData.cpp \
	LDataManager.cpp \
	LEInteractor.cpp \
	LIndexSetData.cpp \
	LIndexSetDataFactory.cpp \
	LIndexSetVariable.cpp \
	LInitStrategy.cpp \
	LMarker.cpp \
	LMesh.cpp \
	LNode.cpp \
	LNodeIndex.cpp \
	LSet.cpp \
	LSetData.cpp \
	LSetDataFactory.cpp \
	LSetDataIterator.cpp \
	LSetVariable.cpp \
	LSiloDataWriter.cpp \
	LTransaction.cpp \
	)

ifneq ($(LIBMESH_LIB),)
srcs-ibtk.cpp += $(call thisdir, \
	FEDataInterpolation.cpp \
	FEDataManager.cpp \
	)
endif

include $(call incsubdirs,fortran)
