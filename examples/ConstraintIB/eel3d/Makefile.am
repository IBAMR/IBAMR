## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

if GSL_ENABLED

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER = test_main.cpp
EXTRA_DIST  = input3d input3d.test
EXTRA_DIST += eelgenerator3d.cpp


EXTRA_PROGRAMS =
GTESTS =
EXAMPLES =
if SAMRAI3D_ENABLED
EXAMPLES += main3d eelgenerator3d
if GTEST_ENABLED
GTESTS += test3d
endif
EXTRA_PROGRAMS += $(GTESTS) $(EXAMPLES)
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = IBEELKinematics3d.cpp $(EXAMPLE_DRIVER)

eelgenerator3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
eelgenerator3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
eelgenerator3d_SOURCES = eelgenerator3d.cpp

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = IBEELKinematics3d.cpp $(GTEST_DRIVER)

examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/eel3d.vertex $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/eel3d.vertex ;\
	fi ;

else
examples:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo " This example needs GNU Scientific Library.                          "
	@echo " Configure IBAMR with additional options --enable-gsl --with-gsl=path"
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d.test $(PWD) ; \
	  cp -f $(srcdir)/eel3d.vertex $(PWD) ; \
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
