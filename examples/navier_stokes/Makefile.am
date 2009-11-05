## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = cell_centered staggered_grid

## Standard make targets.
examples:
	@(cd cell_centered  && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd staggered_grid && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;