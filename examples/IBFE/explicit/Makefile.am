## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ex0 ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 

## Standard make targets.
examples:
	@(cd ex0 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex1 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex2 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex3 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex4 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex5 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex6 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex7 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ex8 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
gtest:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/runTests.sh $(PWD) ; \
	fi ;
	./runTests.sh
