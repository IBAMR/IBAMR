## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = explicit

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

gtest:
	@(cd explicit && make gtest) || exit 1;

gtest-long:
	@(cd explicit && make gtest) || exit 1;

gtest-all:
	@(cd explicit && make gtest-long) || exit 1;
