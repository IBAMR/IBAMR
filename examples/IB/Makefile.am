## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = explicit implicit

## Standard make targets.
examples:
	@(cd explicit && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd implicit && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

gtest: 
	@(cd explicit && make gtest) || exit 1;
	@(cd implicit && make gtest) || exit 1;

gtest-long:
	@(cd explicit && make gtest) || exit 1;
	@(cd implicit && make gtest-long) || exit 1;

gtest-all:
	@(cd explicit && make gtest-long) || exit 1;
	@(cd implicit && make gtest-long) || exit 1;
