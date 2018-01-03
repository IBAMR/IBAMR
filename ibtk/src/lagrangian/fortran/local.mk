ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	lagrangian_delta.F \
	lagrangian_interaction2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	lagrangian_delta.F \
	lagrangian_interaction3d.F \
	)
endif
