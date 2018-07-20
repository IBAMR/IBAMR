## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = IB IBFE IMP adv_diff advect navier_stokes ConstraintIB CIB CIBFE level_set multiphase_flow vc_navier_stokes

all:
	if test "$(CONFIGURATION_BUILD_DIR)" != ""; then \
	  $(MAKE) examples ; \
	fi ;

COMPLETE_EXAMPLES = IB IBFE adv_diff advect navier_stokes ConstraintIB CIB level_set multiphase_flow vc_navier_stokes
# Not all examples can be compiled at the present time
INCOMPLETE_EXAMPLES = IMP CIBFE
.PHONY: examples $(COMPLETE_EXAMPLES)
examples: $(COMPLETE_EXAMPLES)

$(COMPLETE_EXAMPLES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) examples

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
##	@(cd CIBFE            && make gtest) || exit 1;
	@(cd level_set        && make gtest) || exit 1;
	@(cd multiphase_flow  && make gtest) || exit 1;
	@(cd vc_navier_stokes && make gtest) || exit 1;

gtest-long:
	@(cd IB               && make gtest-long) || exit 1;
	@(cd IBFE             && make gtest-long) || exit 1;
##      @(cd IMP              && make gtest-long) || exit 1;
	@(cd adv_diff         && make gtest-long) || exit 1;
	@(cd advect           && make gtest-long) || exit 1;
	@(cd navier_stokes    && make gtest-long) || exit 1;
	@(cd ConstraintIB     && make gtest-long) || exit 1;
	@(cd CIB              && make gtest-long) || exit 1;
##	@(cd CIBFE            && make gtest-long) || exit 1;
	@(cd level_set        && make gtest-long) || exit 1;
	@(cd multiphase_flow  && make gtest-long) || exit 1;
	@(cd vc_navier_stokes && make gtest-long) || exit 1;

gtest-all:
	@(cd IB               && make gtest-all) || exit 1;
	@(cd IBFE             && make gtest-all) || exit 1;
##      @(cd IMP              && make gtest-all) || exit 1;
	@(cd adv_diff         && make gtest-all) || exit 1;
	@(cd advect           && make gtest-all) || exit 1;
	@(cd navier_stokes    && make gtest-all) || exit 1;
	@(cd ConstraintIB     && make gtest-all) || exit 1;
	@(cd CIB              && make gtest-all) || exit 1;
##	@(cd CIBFE            && make gtest-all) || exit 1;
	@(cd level_set        && make gtest-all) || exit 1;
	@(cd multiphase_flow  && make gtest-all) || exit 1;
	@(cd vc_navier_stokes && make gtest-all) || exit 1;

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
