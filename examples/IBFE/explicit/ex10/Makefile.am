## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
	
## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp

EXTRA_DIST = input2d  input2d.test 

EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS =

if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXAMPLES += main2d
if GTEST_ENABLED
GTESTS   += test2d
endif
EXTRA_PROGRAMS += $(EXAMPLES) $(GTESTS)
endif
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)


examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	fi ;

gtest-long:
	make gtest
	./test2d input2d.test
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
	  rm -f $(builddir)/input* ; \
	fi ;
