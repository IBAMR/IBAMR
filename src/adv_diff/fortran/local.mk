ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
    adv_diff_consdiff2d.F \
	)
else
srcs-core.F += $(call thisdir, \
    adv_diff_consdiff3d.F \
	)
endif