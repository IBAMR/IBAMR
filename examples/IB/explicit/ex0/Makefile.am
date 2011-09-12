## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C main_static.C
EXTRA_DIST  = input2d input2d.shell input2d.shell_circum_fibers
EXTRA_DIST += generate_curve2d.m generate_shell2d.m

EXTRA_DIST += curve2d_64.vertex curve2d_128.vertex curve2d_256.vertex curve2d_512.vertex
EXTRA_DIST += curve2d_64.spring curve2d_128.spring curve2d_256.spring curve2d_512.spring
EXTRA_DIST += shell2d_64.vertex shell2d_128.vertex shell2d_256.vertex ##shell2d_512.vertex
EXTRA_DIST += shell2d_64.spring shell2d_128.spring shell2d_256.spring ##shell2d_512.spring
EXTRA_DIST += shell2d_radial_64.vertex shell2d_radial_128.vertex shell2d_radial_256.vertex ##shell2d_radial_512.vertex
EXTRA_DIST += shell2d_radial_64.spring shell2d_radial_128.spring shell2d_radial_256.spring ##shell2d_radial_512.spring

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d main2d_static
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = main.C

main2d_static_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_static_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_static_SOURCES = main_static.C

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d{,.shell,.shell_circum_fibers} $(PWD) ; \
	  cp -f $(srcdir)/curve2d_{64,128,256,512}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{64,128,256}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_radial_{64,128,256}.{vertex,spring} $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d{,.shell,.shell_circum_fibers} ; \
	  rm -f $(builddir)/curve2d_{64,128,256,512}.{vertex,spring} ; \
	  rm -f $(builddir)/shell2d_{64,128,256}.{vertex,spring} ; \
	  rm -f $(builddir)/shell2d_radial_{64,128,256}.{vertex,spring} ; \
	fi ;
