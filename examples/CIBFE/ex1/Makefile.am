## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp
EXTRA_DIST     = input2d input3d README input2d.test input3d.test

EXTRA_PROGRAMS =
EXAMPLES =
GTESTS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
if GTEST_ENABLED
GTESTS   += test2d
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

test2d_CXXFLAGS = $(main2d_CXXFLAGS)
test2d_LDADD = $(main2d_LDADD)
test2d_SOURCES = $(GTEST_DRIVER)

test3d_CXXFLAGS = $(main3d_CXXFLAGS)
test3d_LDADD = $(main3d_LDADD)
test3d_SOURCES = $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/README $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.test $(PWD) ; \
	  cp -f $(srcdir)/input3d.test $(PWD) ; \
	  cp -f $(srcdir)/README $(PWD) ; \
	fi ;

gtest-long:
	make gtest
	./test2d input2d.test
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
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/README ; \
	fi ;
