## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = operator_test solver_test

## Standard make targets.
examples:
	@(cd operator_test && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd solver_test   && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
