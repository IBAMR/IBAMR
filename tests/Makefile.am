## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = Stokes-IB

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) tests ; \
	fi ;

## Standard make targets.
tests:
	@(cd Stokes-IB        && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
