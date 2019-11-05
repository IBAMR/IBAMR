## ---------------------------------------------------------------------
##
## Copyright (c) 2017 - 2019 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ex0 ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 ex10 ex11 ex12 ex13

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

gtest:
	@(cd ex0 && make gtest) || exit 1;
	@(cd ex1 && make gtest) || exit 1;
	@(cd ex2 && make gtest) || exit 1;
	@(cd ex3 && make gtest) || exit 1;
	@(cd ex4 && make gtest) || exit 1;
	@(cd ex5 && make gtest) || exit 1;
	@(cd ex6 && make gtest) || exit 1;
	@(cd ex7 && make gtest) || exit 1;
	@(cd ex8 && make gtest) || exit 1;
	@(cd ex9 && make gtest) || exit 1;
	@(cd ex10 && make gtest) || exit 1;
	@(cd ex11 && make gtest) || exit 1;
	@(cd ex12 && make gtest) || exit 1;
	@(cd ex13 && make gtest) || exit 1;

gtest-long:
	@(cd ex0 && make gtest-long) || exit 1;
	@(cd ex1 && make gtest-long) || exit 1;
	@(cd ex2 && make gtest-long) || exit 1;
	@(cd ex3 && make gtest-long) || exit 1;
	@(cd ex4 && make gtest-long) || exit 1;
	@(cd ex5 && make gtest-long) || exit 1;
	@(cd ex6 && make gtest-long) || exit 1;
	@(cd ex7 && make gtest-long) || exit 1;
	@(cd ex8 && make gtest-long) || exit 1;
	@(cd ex9 && make gtest-long) || exit 1;
	@(cd ex10 && make gtest-long) || exit 1;
	@(cd ex11 && make gtest-long) || exit 1;
	@(cd ex12 && make gtest-long) || exit 1;
	@(cd ex13 && make gtest-long) || exit 1;

gtest-all:
	@(cd ex0 && make gtest-long) || exit 1;
	@(cd ex1 && make gtest-long) || exit 1;
	@(cd ex2 && make gtest-long) || exit 1;
	@(cd ex3 && make gtest-long) || exit 1;
	@(cd ex4 && make gtest-long) || exit 1;
	@(cd ex5 && make gtest-long) || exit 1;
	@(cd ex6 && make gtest-long) || exit 1;
	@(cd ex7 && make gtest-long) || exit 1;
	@(cd ex8 && make gtest-long) || exit 1;
	@(cd ex9 && make gtest-long) || exit 1;
	@(cd ex10 && make gtest-long) || exit 1;
	@(cd ex11 && make gtest-long) || exit 1;
	@(cd ex12 && make gtest-long) || exit 1;
	@(cd ex13 && make gtest-long) || exit 1;
