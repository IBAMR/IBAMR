## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ex0 ex1 ex3

## Standard make targets.
examples:
	@(cd ex0 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex1 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex3 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;


