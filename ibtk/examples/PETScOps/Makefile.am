## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ProlongationMat 

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;


## Standard make targets.
examples:
	@(cd ProlongationMat && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
