ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	cubiccoarsen2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	cubiccoarsen3d.F \
	)
endif
