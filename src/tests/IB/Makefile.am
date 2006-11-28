## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d input3d
EXTRA_DIST += generate_mesh2d.m
EXTRA_DIST += shell2d_64.vertices shell2d_64.edges
EXTRA_DIST += shell2d_128.vertices shell2d_128.edges
EXTRA_DIST += shell2d_256.vertices shell2d_256.edges
EXTRA_DIST += shell2d_512.vertices shell2d_512.edges
EXTRA_DIST += sphere3d.vertices sphere3d.edges

check_PROGRAMS = main2d main3d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{64,128,256,512}.{vertices,edges} $(PWD) ; \
	  cp -f $(srcdir)/sphere3d.{vertices,edges} $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = ../../libIBAMR2d.la ../../libIBAMR.la
main2d_SOURCES = $(SOURCES)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = ../../libIBAMR3d.la ../../libIBAMR.la
main3d_SOURCES = $(SOURCES)
