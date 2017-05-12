ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	linearcfinterpolation2d.F \
	quadcfinterpolation2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	linearcfinterpolation3d.F \
    quadcfinterpolation3d.F \
	)
endif