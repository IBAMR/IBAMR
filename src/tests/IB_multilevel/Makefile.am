## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d
EXTRA_DIST += generate_mesh2d.m
EXTRA_DIST += shell2d_64.vertices shell2d_64.edges
EXTRA_DIST += shell2d_128.vertices shell2d_128.edges
EXTRA_DIST += shell2d_256.vertices shell2d_256.edges
EXTRA_DIST += shell2d_512.vertices shell2d_512.edges
EXTRA_DIST += shell2d_left_64.vertices shell2d_left_64.edges
EXTRA_DIST += shell2d_left_128.vertices shell2d_left_128.edges
EXTRA_DIST += shell2d_left_256.vertices shell2d_left_256.edges
EXTRA_DIST += shell2d_left_512.vertices shell2d_left_512.edges
EXTRA_DIST += shell2d_right_64.vertices shell2d_right_64.edges
EXTRA_DIST += shell2d_right_128.vertices shell2d_right_128.edges
EXTRA_DIST += shell2d_right_256.vertices shell2d_right_256.edges
EXTRA_DIST += shell2d_right_512.vertices shell2d_right_512.edges

check_PROGRAMS = main2d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{64,128,256,512}.{vertices,edges} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_left_{64,128,256,512}.{vertices,edges} $(PWD) ; \
	  cp -f $(srcdir)/shell2d_right_{64,128,256,512}.{vertices,edges} $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS2d.a $(top_builddir)/stools/src/libSTOOLS.a
main2d_SOURCES = $(SOURCES)
