## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidDensity.h SetFluidGasSolidViscosity.cpp SetFluidGasSolidViscosity.h LSLocateGasInterface.cpp LSLocateGasInterface.h  TagLSRefinementCells.cpp TagLSRefinementCells.h SetLSProperties.cpp SetLSProperties.h GravityForcing.cpp GravityForcing.h FlowGravityForcing.cpp FlowGravityForcing.h
GTEST_DRIVER   = test_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidDensity.h SetFluidGasSolidViscosity.cpp SetFluidGasSolidViscosity.h LSLocateGasInterface.cpp LSLocateGasInterface.h  TagLSRefinementCells.cpp TagLSRefinementCells.h SetLSProperties.cpp SetLSProperties.h GravityForcing.cpp GravityForcing.h FlowGravityForcing.cpp FlowGravityForcing.h
EXTRA_DIST     = input2d input2d.test CircleMesh.msh

GTESTS   =
EXAMPLES =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
if GTEST_ENABLED
GTESTS   += test2d
endif
EXAMPLES += main2d
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[2]d* $(PWD) ; \
	  cp -f $(srcdir)/*.msh $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[2]d* $(PWD) ; \
	  cp -f $(srcdir)/*.msh $(PWD) ; \
	fi ;
	./test2d input2d.test

gtest-long:
	make gtest
	./test2d input2d.test

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

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input[2]d ; \
	  rm -f $(builddir)/input[2]d.test ; \
	  rm -f $(builddir)/*.msh ; \
	fi ;