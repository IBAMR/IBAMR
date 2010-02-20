## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = approx_projection exact_projection

## Standard make targets.
examples:
	@(cd approx_projection && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd exact_projection && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
