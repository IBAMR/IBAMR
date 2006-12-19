## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C UInit.C UInit.h
EXTRA_DIST = input2d ##input3d

check_PROGRAMS = main2d ##main3d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
##	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS2d.a $(top_builddir)/stools/src/libSTOOLS.a
main2d_SOURCES = $(SOURCES)

##main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
##main3d_LDADD = $(top_builddir)/src/libIBAMR3d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS3d.a $(top_builddir)/stools/src/libSTOOLS.a
##main3d_SOURCES = $(SOURCES)
