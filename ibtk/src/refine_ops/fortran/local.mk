ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	cart_side_refine2d.F \
	divpreservingrefine2d.F \
	)
else
srcs-core.F += $(call thisdir, \
    cart_side_refine3d.F \
	divpreservingrefine3d.F \
	)
endif
