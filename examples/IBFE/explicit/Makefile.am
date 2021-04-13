## ---------------------------------------------------------------------
##
## Copyright (c) 2010 - 2019 by the IBAMR developers
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
SUBDIRS = ex0 ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 ex10 ex11 ex12

COMPLETE_EXAMPLES = $(SUBDIRS)
INCOMPLETE_EXAMPLES =
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples



