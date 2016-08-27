## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER    = test_main.cpp

EXTRA_DIST = input2d input3d input2d.test input3d.test

EXAMPLES       =
GTESTS         =
EXTRA_PROGRAMS =
if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXAMPLES      += main2d convergence_tester2d
GTESTS        += test2d
EXTRA_PROGRAMS += $(EXAMPLES) $(GTESTS)
endif
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

convergence_tester2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
convergence_tester2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
convergence_tester2d_SOURCES = convergence_tester.cpp

test2d_CXXFLAGS = $(main2d_CXXFLAGS)
test2d_LDADD = $(main2d_LDADD)
test2d_SOURCES = $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;
	./test2d input2d.test --gtest_filter=*.2d
clean-local:
	rm -f $(EXTRA_PROGRAMS) $(GTEST_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input* ; \
	fi ;
