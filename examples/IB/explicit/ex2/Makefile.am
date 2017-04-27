## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXAMPLE_DRIVER = ex_main.cpp
GTEST_DRIVER   = test_main.cpp
EXTRA_DIST  = input2d.cylinder input2d.cylinder_stabilized input3d.cylinder input3d.sphere
EXTRA_DIST += input2d.test input3d.cylinder.test
EXTRA_DIST += cylinder2d_128.vertex cylinder2d_256.vertex cylinder2d_512.vertex cylinder2d_1024.vertex
EXTRA_DIST += cylinder2d_128.spring cylinder2d_256.spring cylinder2d_512.spring cylinder2d_1024.spring
EXTRA_DIST += cylinder3d_128.vertex cylinder3d_256.vertex cylinder3d_512.vertex
EXTRA_DIST += sphere3d_32.vertex sphere3d_64.vertex sphere3d_128.vertex sphere3d_256.vertex
EXTRA_DIST += sphere3d_32.spring sphere3d_64.spring sphere3d_128.spring sphere3d_256.spring

EXAMPLES =
GTESTS   =
EXTRA_PROGRAMS =
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

test2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
test2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
test2d_SOURCES = $(GTEST_DRIVER)

test3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
test3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
test3d_SOURCES = $(GTEST_DRIVER)


examples: $(EXAMPLES)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input{2d.cylinder,2d.cylinder_stabilized,3d.cylinder,3d.sphere} $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d_{128,256,512,1024}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/cylinder3d_{128,256,512}.vertex $(PWD) ; \
	  cp -f $(srcdir)/sphere3d_{32,64,128,256}.{vertex,spring} $(PWD) ; \
	fi ;

if GTEST_ENABLED
gtest: $(GTESTS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input* $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d_{128,256,512,1024}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/cylinder3d_{128,256,512}.vertex $(PWD) ; \
	  cp -f $(srcdir)/sphere3d_{32,64,128,256}.{vertex,spring} $(PWD) ; \
	fi ;
	./test2d input2d.test

gtest-long:
	make gtest
	./test3d input3d.cylinder.test

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
	  rm -f $(builddir)/input{2d.cylinder,2d.cylinder_stabilized,3d.cylinder,3d.sphere} ; \
	  rm -f $(builddir)/cylinder2d_{128,256,512,1024}.{vertex,spring} ; \
	  rm -f $(builddir)/cylinder3d_{128,256,512}.vertex ; \
	  rm -f $(builddir)/sphere3d_{32,64,128,256}.{vertex,spring} ; \
	fi ;
