## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = ibtk src lib
ACLOCAL_AMFLAGS = -I m4

## Standard make targets.
lib: all
examples: lib
	@(cd examples && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

MPIEXEC = @MPIEXEC@
NUMDIFF = @NUMDIFF@
# attest uses parameters set in attest.conf; said parameters are overridden by
# command line arguments. Use configuration info to set up the path to mpirun:
attest.conf:
	@printf "# attest default configuration file\n" > attest.conf-tmp
	@printf "[attest]\n" >> attest.conf-tmp
	@printf "jobs = 1\n" >> attest.conf-tmp
	@printf "keep_work_directories = False\n"  >> attest.conf-tmp
	@printf "mpiexec = ${MPIEXEC}\n" >> attest.conf-tmp
	@printf "numdiff = ${NUMDIFF}\n" >> attest.conf-tmp
	@printf "show_only = False\n" >> attest.conf-tmp
	@printf "test_directory = /tests/\n" >> attest.conf-tmp
	@printf "test_regex = .*\n" >> attest.conf-tmp
	@mv attest.conf-tmp attest.conf

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
