ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	patchsmoothers2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	patchsmoothers3d.F \
	)
endif
