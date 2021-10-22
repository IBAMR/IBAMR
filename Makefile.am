## ---------------------------------------------------------------------
##
## Copyright (c) 2006 - 2019 by the IBAMR developers
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
SUBDIRS = ibtk src lib
ACLOCAL_AMFLAGS = -I m4

## Standard make targets.
lib: all
examples: lib
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) $@

MPIEXEC = @MPIEXEC@
NUMDIFF = @NUMDIFF@
# attest uses parameters set in attest.conf; said parameters are overridden by
# command line arguments. Use configuration info to set up the path to mpirun:
attest.conf: $(abs_builddir)/tests/attest.conf
	ln -f -s ./tests/attest.conf $(abs_builddir)

tests: lib attest.conf
	ln -f -s $(top_srcdir)/attest $(abs_builddir)
	@cd $@ && $(MAKE) $(AM_MAKEFLAGS) $@
.PHONY: tests

install-exec-local:
	@echo "===================================================================="
	@echo "===================================================================="
	@echo "WARNING: make install is not officially supported by IBAMR!"
	@echo "         use at your own risk!"
	@echo "===================================================================="
	@echo "===================================================================="

## Add various "extra" files to the IBAMR distribution.
EXTRA_DIST  = $(top_srcdir)/COPYRIGHT
EXTRA_DIST += $(top_srcdir)/autogen.sh
EXTRA_DIST += $(top_srcdir)/doc/Doxyfile.in
EXTRA_DIST += $(top_srcdir)/doc/footer.html
EXTRA_DIST += $(top_srcdir)/doc/header.html
EXTRA_DIST += $(top_srcdir)/doc/libstdc++.tag
EXTRA_DIST += $(top_srcdir)/doc/mainpage.dox
EXTRA_DIST += $(top_srcdir)/doc/online_Doxyfile.in
EXTRA_DIST += $(top_srcdir)/doc/samrai.tag
EXTRA_DIST += $(top_srcdir)/doc/stylesheet.css

## Refresh emacs TAGS file.
refresh_tags:
	cd ibtk && $(MAKE) refresh_tags
	cd $(top_srcdir) && find . -name "*.[ChIf]" -print | etags -
