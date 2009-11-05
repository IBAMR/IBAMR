## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C UInit.C UInit.h
EXTRA_DIST = input2d

bin_PROGRAMS =
if SAMRAI2D_ENABLED
bin_PROGRAMS += main2d
endif

all-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/lib/libIBAMR2d.a $(top_builddir)/lib/libIBAMR.a $(top_builddir)/ibtk/lib/libIBTK2d.a $(top_builddir)/ibtk/lib/libIBTK.a $(top_builddir)/ibtk/third_party/qd-2.3.7/fortran/libqdmod.a $(top_builddir)/ibtk/third_party/qd-2.3.7/src/libqd.a
main2d_SOURCES = $(SOURCES)
