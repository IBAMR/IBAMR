## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp SetFluidProperties.cpp SetLSProperties.cpp LSLocateColumnInterface.cpp GravityForcing.cpp TagLSRefinementCells.cpp
GTEST_DRIVER   = test_main.cpp SetFluidProperties.cpp SetLSProperties.cpp LSLocateColumnInterface.cpp GravityForcing.cpp TagLSRefinementCells.cpp
EXTRA_DIST     = input2d input2d.test
EXTRA_DIST    += input2d_CaseA input2d_CaseB input2d_CaseC input2d_5thOrderStokes
EXTRA_DIST    += input3d input3d.test
EXTRA_DIST    += compareStokesFirstOrder.m compareStokesSecondOrder.m compareStokesFifthOrder.m

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
if SAMRAI3D_ENABLED
if GTEST_ENABLED
GTESTS   += test3d
endif
EXAMPLES += main3d
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d* $(PWD) ; \
	  cp -f $(srcdir)/*.m $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d* $(PWD) ; \
	fi ;
	./test2d input2d.test

gtest-long:
	make gtest
	./test3d input3d.test

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
	  rm -f $(builddir)/input[23]d               ; \
	  rm -f $(builddir)/input[23]d.test          ; \
	  rm -f $(builddir)/input[2]d_CaseA          ; \
	  rm -f $(builddir)/input[2]d_CaseB          ; \
	  rm -f $(builddir)/input[2]d_CaseC          ; \
	  rm -f $(builddir)/input[2]d_5thOrderStokes ; \
	fi ;