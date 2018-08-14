## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ex0 ex1

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

gtest:

# tests are disabled for CIBFE
#	@(cd ex0 && make gtest) || exit 1;
#	@(cd ex1 && make gtest) || exit 1;

gtest-long:

#	@(cd ex0 && make gtest) || exit 1;
#	@(cd ex1 && make gtest-long) || exit 1;

gtest-all:

#	@(cd ex0 && make gtest-long) || exit 1;
#	@(cd ex1 && make gtest-long) || exit 1;
