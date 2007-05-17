## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d
EXTRA_DIST += generate_mesh2d.m

EXTRA_DIST += shell2d_crse_128.vertex shell2d_crse_128.spring
EXTRA_DIST += shell2d_fine_128.vertex shell2d_fine_128.spring
EXTRA_DIST += shell2d_crse_256.vertex shell2d_crse_256.spring
EXTRA_DIST += shell2d_fine_256.vertex shell2d_fine_256.spring

if DATA_FILES_ENABLED
EXTRA_DIST += shell2d_crse_512.vertex shell2d_crse_512.spring
EXTRA_DIST += shell2d_fine_512.vertex shell2d_fine_512.spring
endif

check_PROGRAMS =
if SAMRAI2D_ENABLED
check_PROGRAMS += main2d
endif

check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{crse,fine}_{128,256}.{vertex,spring} $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS2d.a $(top_builddir)/stools/src/libSTOOLS.a
main2d_SOURCES = $(SOURCES)
