## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

GTESTdir = $(top_builddir)/examples/IBFE/explicit/ex4/

## Dimension-dependent testers
EXAMPLE_DRIVER = ex4_main.cpp
TEST_DRIVER = test_main.cpp

EXTRA_DIST = input2d input3d input2d.test input3d.test

EXTRA_PROGRAMS =
GTEST_PROGRAMS =
if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d convergence_tester2d
GTEST_PROGRAMS += test2d
endif
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d convergence_tester3d
GTEST_PROGRAMS += test3d
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
test2d_SOURCES = $(TEST_DRIVER)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

convergence_tester3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
convergence_tester3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
convergence_tester3d_SOURCES = convergence_tester.cpp

test3d_CXXFLAGS = $(main3d_CXXFLAGS)
test3d_LDADD = $(main3d_LDADD)
test3d_SOURCES = $(TEST_DRIVER)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

gtest: $(GTEST_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS) $(GTEST_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input* ; \
	fi ;
