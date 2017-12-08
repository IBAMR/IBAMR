## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp
EXTRA_DIST  = input2d input2d.test
EXTRA_DIST += input3d input3d.test
EXTRA_DIST += points2d.vertex
EXTRA_DIST += points3d.vertex


EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
if GTEST_ENABLED
GTESTS += test2d
endif
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif
if SAMRAI3D_ENABLED
EXAMPLES += main3d
if GTEST_ENABLED
GTESTS   += test3d
endif
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
	  cp -f $(srcdir)/input{2,3}d $(PWD) ; \
	  cp -f $(srcdir)/points{2,3}d.vertex $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input{2,3}d* $(PWD) ; \
	  cp -f $(srcdir)/points{2,3}d.vertex $(PWD) ; \
	fi ;

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
	  rm -f $(builddir)/input{2,3}d ; \
	  rm -f $(builddir)/points{2,3}d.vertex ; \
	fi ;
