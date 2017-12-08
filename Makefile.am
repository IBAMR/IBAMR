## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ibtk src lib
ACLOCAL_AMFLAGS = -I m4

## Standard make targets.
lib: all
examples: lib
	@(cd examples && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

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



if GTEST_ENABLED
gtest:
	@(cd examples && make gtest) || exit 1;
	@(cd ibtk     && make gtest) || exit 1;

gtest-long:
	@(cd examples && make gtest-long) || exit 1;
	@(cd ibtk     && make gtest-long) || exit 1;

gtest-all:
	@(cd examples && make gtest-all) || exit 1;
	@(cd ibtk     && make gtest-all) || exit 1;
else
gtest:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo " To run tests, Google Test Framework must be enabled.                "
	@echo " Configure IBAMR with additional options:                            "
	@echo "                      --enable-gtest --with-gtest=path               "
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
gtest-long:
	make gtest
endif
