## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateStructureInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
GTEST_DRIVER   = test_main.cpp SetFluidGasSolidDensity.cpp SetFluidGasSolidViscosity.cpp LSLocateStructureInterface.cpp LSLocateGasInterface.cpp RigidBodyKinematics.cpp TagLSRefinementCells.cpp SetLSProperties.cpp GravityForcing.cpp
EXTRA_DIST     = input2d input2d.test input3d

EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS = wedgeGen
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

wedgeGen_CXXFLAGS = $(AM_CXXFLAGS)
wedgeGen_LDADD =
wedgeGen_SOURCES = wedgeGen.cpp

# wedgeGen creates both 2d and 3d input files
wedge2d.vertex: wedgeGen
	./wedgeGen

wedge3d.vertex: wedge2d.vertex

examples: $(EXAMPLES) wedge2d.vertex wedge3d.vertex
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d* $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.test $(PWD) ; \
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
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/wedge2d.vertex ; \
	  rm -f $(builddir)/Wedge2d.m ; \
	  rm -f $(builddir)/Wedge3d.m ; \
	fi ;
