## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
	
## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp

EXTRA_DIST = input2d input3d input2d.test input3d.test

EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS =

if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXAMPLES += main2d
GTESTS   += test2d
EXTRA_PROGRAMS += $(EXAMPLES) $(GTESTS)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d 
GTESTS   += test3d
EXTRA_PROGRAMS += $(EXAMPLES) $(GTESTS)
endif
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
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;
	./test2d input2d.test --gtest_filter=*.2d
	./test3d input3d.test --gtest_filter=*.3d

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input* ; \
	fi ;
