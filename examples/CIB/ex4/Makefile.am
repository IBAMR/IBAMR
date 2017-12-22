## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp
EXTRA_DIST  = input3d.sphere sphere_3d.vertex input3d.test.long
EXTRA_DIST += input3d.test

EXTRA_PROGRAMS =
GTESTS = 
EXAMPLES =
if SAMRAI3D_ENABLED
EXAMPLES += main3d
if GTEST_ENABLED
GTESTS   += test3d
endif
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(EXAMPLE_DRIVER)

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d.sphere $(PWD) ; \
	  cp -f $(srcdir)/sphere_3d.vertex $(PWD) ; \
	  cp -f $(srcdir)/check_sol.m $(PWD) ; \
	  cp -f $(srcdir)/petsc_options.dat $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d* $(PWD) ; \
	  cp -f $(srcdir)/sphere_3d.vertex $(PWD) ; \
	  cp -f $(srcdir)/check_sol.m $(PWD) ; \
	  cp -f $(srcdir)/petsc_options.dat $(PWD) ; \
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

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d.sphere* ; \
          rm -f $(builddir)/sphere_3d.vertex ; \
	  rm -f $(builddir)/check_sol.m ; \
	  rm -f $(builddir)/petsc_options.dat ; \
	fi ;
