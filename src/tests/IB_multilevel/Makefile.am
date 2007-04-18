## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d
EXTRA_DIST += generate_mesh2d.m
EXTRA_DIST += shell2d_crse_128.vertex shell2d_crse_128.spring shell2d_fine_128.vertex shell2d_fine_128.spring
EXTRA_DIST += shell2d_crse_256.vertex shell2d_crse_256.spring shell2d_fine_256.vertex shell2d_fine_256.spring
EXTRA_DIST += shell2d_crse_512.vertex shell2d_crse_512.spring shell2d_fine_512.vertex shell2d_fine_512.spring
EXTRA_DIST += shell2d_crse_1024.vertex shell2d_crse_1024.spring shell2d_fine_1024.vertex shell2d_fine_1024.spring

check_PROGRAMS = main2d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input[23]d $(PWD) ; \
	  cp -f $(srcdir)/shell2d_{crse,fine}_{128,256,512,1024}.{vertex,spring} $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS2d.a $(top_builddir)/stools/src/libSTOOLS.a
main2d_SOURCES = $(SOURCES)
