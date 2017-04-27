## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = IB IBFE IMP adv_diff advect navier_stokes ConstraintIB CIB CIBFE

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

## Standard make targets.
examples:
	@(cd IB               && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd IBFE             && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd IMP              && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd adv_diff         && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd advect           && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd navier_stokes    && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd ConstraintIB     && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd CIB              && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd CIBFE            && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

if GTEST_ENABLED
gtest:
	@(cd IB               && make gtest) || exit 1;
	@(cd IBFE             && make gtest) || exit 1;
##      @(cd IMP              && make gtest) || exit 1;
	@(cd adv_diff         && make gtest) || exit 1;
	@(cd advect           && make gtest) || exit 1;
	@(cd navier_stokes    && make gtest) || exit 1;
	@(cd ConstraintIB     && make gtest) || exit 1;
	@(cd CIB              && make gtest) || exit 1;
	@(cd CIBFE            && make gtest) || exit 1;

gtest-long:
	@(cd IB               && make gtest-long) || exit 1;
	@(cd IBFE             && make gtest-long) || exit 1;
##      @(cd IMP              && make gtest-long) || exit 1;
	@(cd adv_diff         && make gtest-long) || exit 1;
	@(cd advect           && make gtest-long) || exit 1;
	@(cd navier_stokes    && make gtest-long) || exit 1;
	@(cd ConstraintIB     && make gtest-long) || exit 1;
	@(cd CIB              && make gtest-long) || exit 1;
	@(cd CIBFE            && make gtest-long) || exit 1;

gtest-all:
	@(cd IB               && make gtest-all) || exit 1;
	@(cd IBFE             && make gtest-all) || exit 1;
##      @(cd IMP              && make gtest-all) || exit 1;
	@(cd adv_diff         && make gtest-all) || exit 1;
	@(cd advect           && make gtest-all) || exit 1;
	@(cd navier_stokes    && make gtest-all) || exit 1;
	@(cd ConstraintIB     && make gtest-all) || exit 1;
	@(cd CIB              && make gtest-all) || exit 1;
	@(cd CIBFE            && make gtest-all) || exit 1;
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
