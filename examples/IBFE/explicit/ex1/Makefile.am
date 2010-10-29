## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
EXTRA_DIST = input2d input3d

EXTRA_PROGRAMS =
if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d convergence_tester2d
endif
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/lib/libIBAMR2d.a $(top_builddir)/lib/libIBAMR.a $(top_builddir)/ibtk/lib/libIBTK2d.a $(top_builddir)/ibtk/lib/libIBTK.a
main2d_SOURCES = main.C

convergence_tester2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
convergence_tester2d_LDADD = $(top_builddir)/lib/libIBAMR2d.a $(top_builddir)/lib/libIBAMR.a $(top_builddir)/ibtk/lib/libIBTK2d.a $(top_builddir)/ibtk/lib/libIBTK.a
convergence_tester2d_SOURCES = convergence_tester.C

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	fi ;
