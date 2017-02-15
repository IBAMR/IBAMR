## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = eel2d eel3d

## Standard make targets.
examples:
	@(cd eel2d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd eel3d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd falling_sphere && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd flow_past_cylinder && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd knifefish && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd impulsively_started_cylinder && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd moving_plate && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd oscillating_rigid_cylinder && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd stokes_first_problem && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
