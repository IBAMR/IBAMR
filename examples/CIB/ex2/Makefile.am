## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER = test_main.cpp
EXTRA_DIST  = input2d.plate input2d.plate.amr plate2d.vertex plate2d.m petsc_options.dat
EXTRA_DIST += input2d.test

EXTRA_PROGRAMS =
EXAMPLES = 
GTESTS =
if SAMRAI2D_ENABLED
EXAMPLES += main2d
if GTEST_ENABLED
GTESTS += test2d
endif
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(EXAMPLE_DRIVER)

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d* $(PWD) ; \
	  cp -f $(srcdir)/*.vertex $(PWD) ; \
	  cp -f $(srcdir)/*.m $(PWD) ; \
	  cp -f $(srcdir)/*.dat $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d* $(PWD) ; \
	  cp -f $(srcdir)/*.vertex $(PWD) ; \
	  cp -f $(srcdir)/*.m $(PWD) ; \
	  cp -f $(srcdir)/*.dat $(PWD) ; \
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
	  rm -f $(builddir)/input2d.plate ; \
	  rm -f $(builddir)/input2d.plate.amr ; \
          rm -f $(builddir)/plate2d.vertex ; \
          rm -f $(builddir)/plate2d.m ; \
          rm -f $(builddir)/petsc_options.dat ; \
	fi ;
