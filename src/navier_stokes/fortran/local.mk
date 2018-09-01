ifeq ($(DIM),2)
srcs-core.F += $(call thisdir, \
    navier_stokes_bdryop2d.F \
    navier_stokes_divsource2d.F \
    navier_stokes_stabledt2d.F \
    navier_stokes_staggered_derivatives2d.F \
    navier_stokes_staggered_helpers2d.F \
    navier_stokes_stochastic_forcing2d.F \
    navier_stokes_surface_tension_forcing2d.F \
	)
else
srcs-core.F += $(call thisdir, \
    navier_stokes_bdryop3d.F \
    navier_stokes_divsource3d.F \
    navier_stokes_stabledt3d.F \
    navier_stokes_staggered_derivatives3d.F \
    navier_stokes_staggered_helpers3d.F \
    navier_stokes_stochastic_forcing3d.F \
    navier_stokes_surface_tension_forcing3d.F \
	)
endif
