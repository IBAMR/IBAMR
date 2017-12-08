## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = explicit

## Standard make targets.
examples:
	@(cd explicit && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

gtest:
	@(cd explicit && make gtest) || exit 1;

gtest-long:
	@(cd explicit && make gtest-long) || exit 1;

gtest-all:
	@(cd explicit && make gtest-all) || exit 1;
