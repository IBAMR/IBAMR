## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = explicit

## Standard make targets.
examples:
	@(cd explicit && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
