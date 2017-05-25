ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	cubiccoarsen2d.F \
	rt0coarsen2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	cubiccoarsen3d.F \
	rt0coarsen3d.F \
	)
endif
