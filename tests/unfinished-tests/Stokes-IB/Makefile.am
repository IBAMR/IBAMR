## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2019 by the IBAMR developers
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
SUBDIRS = test0 test1 test2

## Standard make targets.
tests:
	@(cd test0 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd test1 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd test2 && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
