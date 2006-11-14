## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C UInit.C UInit.h
EXTRA_DIST = input2d input3d cylinder2d.256 cylinder2d.512 cylinder2d.1024 cylinder2d.2048

check_PROGRAMS = main2d main3d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d.{256,512,1024,2048} $(PWD) ; \
	  cp -f $(srcdir)/sphere3d_vertices.dat $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = ../../libIBAMR2d.la ../../libIBAMR.la
main2d_SOURCES = $(SOURCES)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = ../../libIBAMR3d.la ../../libIBAMR.la
main3d_SOURCES = $(SOURCES)
