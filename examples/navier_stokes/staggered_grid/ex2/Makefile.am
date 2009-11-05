## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST = input2d

check_PROGRAMS =
if SAMRAI2D_ENABLED
check_PROGRAMS += main2d
endif

check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/ibtk/src/libIBTK2d.a $(top_builddir)/ibtk/src/libIBTK.a $(top_builddir)/ibtk/third_party/qd-2.3.7/fortran/libqdmod.a $(top_builddir)/ibtk/third_party/qd-2.3.7/src/libqd.a
main2d_SOURCES = $(SOURCES)
