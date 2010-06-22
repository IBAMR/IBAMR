## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input3d
EXTRA_DIST += generate_mesh3d.m

EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(top_builddir)/lib/libIBAMR3d.a $(top_builddir)/lib/libIBAMR.a $(top_builddir)/ibtk/lib/libIBTK3d.a $(top_builddir)/ibtk/lib/libIBTK.a $(top_builddir)/ibtk/third_party/qd-2.3.7/fortran/libqdmod.a $(top_builddir)/ibtk/third_party/qd-2.3.7/src/libqd.a
main3d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	fi ;
