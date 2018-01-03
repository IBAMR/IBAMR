ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	curl2d.F \
    div2d.F \
    flux2d.F \
    grad2d.F \
    graddetect2d.F \
    interp2d.F \
    laplace2d.F \
    miscmath2d.F \
    rot2d.F \
    strain2d.F \
    vclaplace2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	curl3d.F \
    div3d.F \
    flux3d.F \
    grad3d.F \
    graddetect3d.F \
    interp3d.F \
    laplace3d.F \
    miscmath3d.F \
    rot3d.F \
    strain3d.F \
	)
endif
