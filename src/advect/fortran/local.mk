ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
    advect_centered_derivatives2d.F \
    advect_detect2d.F \
    advect_diff2d.F \
    advect_helpers.F \
    advect_predictors2d.F \
    advect_stable2d.F \
	)
else
srcs-core.F += $(call thisdir, \
    advect_centered_derivatives3d.F \
    advect_detect3d.F \
    advect_diff3d.F \
    advect_helpers.F \
    advect_predictors3d.F \
    advect_stable3d.F \
	)
endif
