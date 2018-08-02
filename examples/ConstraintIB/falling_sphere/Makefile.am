## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = CartGridBodyForce.cpp ForceProjector.cpp RigidBodyKinematics.cpp ex_main.cpp
GTEST_DRIVER = CartGridBodyForce.cpp ForceProjector.cpp RigidBodyKinematics.cpp test_main.cpp
EXTRA_DIST  = input3d 
EXTRA_DIST += sphereGen3d.cpp


EXTRA_PROGRAMS =
GTESTS =
EXAMPLES =
if SAMRAI3D_ENABLED
if GTEST_ENABLED
GTESTS += test3d
endif
EXAMPLES += main3d sphereGen3d
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

sphereGen3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
sphereGen3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
sphereGen3d_SOURCES = sphereGen3d.cpp

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = $(GTEST_DRIVER)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/sphere3d.vertex $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/sphere3d.vertex ;\
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d.test $(PWD) ; \
	  cp -f $(srcdir)/sphere3d.vertex $(PWD) ; \
	fi ;

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
