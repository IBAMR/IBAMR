## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d input3d
EXTRA_DIST += generate_mesh2d.m
EXTRA_DIST += curve2d_64.vertex curve2d_64.spring
EXTRA_DIST += curve2d_128.vertex curve2d_128.spring
EXTRA_DIST += curve2d_256.vertex curve2d_256.spring
EXTRA_DIST += curve2d_512.vertex curve2d_512.spring
EXTRA_DIST += shell2d_64.vertex shell2d_64.spring
EXTRA_DIST += shell2d_128.vertex shell2d_128.spring
EXTRA_DIST += shell2d_256.vertex shell2d_256.spring
EXTRA_DIST += shell2d_512.vertex shell2d_512.spring
EXTRA_DIST += sphere3d.vertex sphere3d.spring

check_PROGRAMS = main2d main3d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	  cp -f $(srcdir)/curve2d_{64,128,256,512}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{64,128,256,512}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/sphere3d.{vertex,spring} $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS2d.a $(top_builddir)/stools/src/libSTOOLS.a
main2d_SOURCES = $(SOURCES)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(top_builddir)/src/libIBAMR3d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS3d.a $(top_builddir)/stools/src/libSTOOLS.a
main3d_SOURCES = $(SOURCES)
