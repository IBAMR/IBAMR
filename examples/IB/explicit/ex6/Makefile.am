## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp IBSimpleHierarchyIntegrator.cpp IBSimpleHierarchyIntegrator.h
EXTRA_DIST  = input2d
EXTRA_DIST += generate_curve2d.m

EXTRA_DIST += curve2d_64.vertex curve2d_128.vertex curve2d_256.vertex curve2d_512.vertex
EXTRA_DIST += curve2d_64.spring curve2d_128.spring curve2d_256.spring curve2d_512.spring

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	  cp -f $(srcdir)/curve2d_{64,128,256,512}.{vertex,spring} $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/curve2d_{64,128,256,512}.{vertex,spring} ; \
	fi ;
