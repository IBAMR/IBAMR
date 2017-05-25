ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
	cartphysbdryop2d.F \
	)
else
srcs-core.F += $(call thisdir, \
	cartphysbdryop3d.F \
	)
endif
