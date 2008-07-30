## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d.cylinder
EXTRA_DIST += generate_cylinder2d.m compute_St.m

EXTRA_DIST += cylinder2d_128.vertex cylinder2d_256.vertex cylinder2d_512.vertex cylinder2d_1024.vertex
EXTRA_DIST += cylinder2d_128.spring cylinder2d_256.spring cylinder2d_512.spring cylinder2d_1024.spring

check_PROGRAMS =
if SAMRAI2D_ENABLED
check_PROGRAMS += main2d
endif

check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.cylinder $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d_{128,256,512,1024}.{vertex,spring} $(PWD) ; \
          cp -f $(srcdir)/compute_St.m $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/ibtk/src/libIBTK2d.a $(top_builddir)/ibtk/src/libIBTK.a
main2d_SOURCES = $(SOURCES)
