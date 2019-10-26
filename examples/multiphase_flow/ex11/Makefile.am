## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateStructureInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
GTEST_DRIVER   = test_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateStructureInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
EXTRA_DIST     = input2d input2d.test

EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
if GTEST_ENABLED
GTESTS   += test2d
endif
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES) bargeGen
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

bargeGen_CXXFLAGS = $(AM_CXXFLAGS)
bargeGen_LDADD =
bargeGen_SOURCES = bargeGen.cpp

barge2d.vertex: bargeGen
	./bargeGen

examples: $(EXAMPLES) barge2d.vertex
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[2]d* $(PWD) ; \
	  cp -f $(srcdir)/Barge2d.m $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.test $(PWD) ; \
	  cp -f $(srcdir)/barge2d.vertex $(PWD) ; \
	fi ;
	./test2d input2d.test

gtest-long:
	make gtest

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
	  rm -f $(builddir)/input2d.test ; \
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/barge2d.vertex ; \
	  rm -f $(builddir)/bargeGen ; \
	  rm -f $(builddir)/Barge2d.m ; \
	fi ;
