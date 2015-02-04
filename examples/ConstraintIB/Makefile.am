## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = eel2d eel3d

## Standard make targets.
examples:
	@(cd eel2d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd eel3d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

