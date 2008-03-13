## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C GravitationalBodyForce.C GravitationalBodyForce.h PressureBcCoefs.C PressureBcCoefs.h
EXTRA_DIST  = input2d
EXTRA_DIST +=  fila_256.vertex  fila_256.spring  fila_256.beam  fila_256.target
EXTRA_DIST +=  fila_512.vertex  fila_512.spring  fila_512.beam  fila_512.target
EXTRA_DIST += fila_1024.vertex fila_1024.spring fila_1024.beam fila_1024.target

check_PROGRAMS =
if SAMRAI2D_ENABLED
check_PROGRAMS += main2d
endif

check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	  cp -f $(srcdir)/fila_{256,512,1024}.{vertex,spring,beam,target} $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/ibtk/src/libIBTK2d.a $(top_builddir)/ibtk/src/libIBTK.a
main2d_SOURCES = $(SOURCES)
