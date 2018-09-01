ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
    adv_diff_consdiff2d.F \
    adv_diff_wp_convective_op2d.F \
	)
else
srcs-core.F += $(call thisdir, \
    adv_diff_consdiff3d.F \
    adv_diff_wp_convective_op3d.F \
	)
endif
