## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = test0 test1 test2

## Standard make targets.
tests:
	@(cd test0 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd test1 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd test2 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
