## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateCircularInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
GTEST_DRIVER   = test_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateCircularInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
EXTRA_DIST     = input2d.cyl input2d.test
EXTRA_DIST    += input3d


EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
if GTEST_ENABLED
GTESTS   += test2d
endif
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif
if SAMRAI3D_ENABLED
if GTEST_ENABLED
GTESTS   += test3d
endif
EXAMPLES += main3d
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES) sphereGen
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

sphereGen_CXXFLAGS = $(AM_CXXFLAGS)
sphereGen_LDADD =
sphereGen_SOURCES = sphereGen.cpp

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = $(GTEST_DRIVER)

sphere3d.vertex: sphereGen
	./sphereGen

examples: $(EXAMPLES) sphere3d.vertex
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d* $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d.vertex $(PWD) ; \
	  cp -f $(srcdir)/Cylinder2d_cubic.m $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.test $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d.vertex $(PWD) ; \
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
	  rm -f $(builddir)/input2d.cyl ; \
	  rm -f $(builddir)/input2d.test ; \
	  rm -f $(builddir)/input3d.sphere ; \
	  rm -f $(builddir)/cylinder2d.vertex ; \
	  rm -f $(builddir)/Cylinder2d_cubic.m ; \
	  rm -f $(builddir)/sphereGen.cpp ; \
	fi ;
