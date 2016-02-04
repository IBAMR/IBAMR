## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp
EXTRA_DIST = input2d.shell README

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(SOURCES)

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.shell $(PWD) ; \
	  cp -f $(srcdir)/shell2d_64.vertex $(PWD) ; \
          cp -f $(srcdir)/shell2d_64.spring $(PWD) ; \
          cp -f $(srcdir)/generate_shell2d.m $(PWD) ; \
          cp -f $(srcdir)/README $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d.shell ; \
	  rm -f $(builddir)/shell2d_64.vertex ; \
          rm -f $(builddir)/shell2d_64.spring ; \
          rm -f $(builddir)/generate_shell2d.m ; \
          rm -f $(builddir)/README ; \
	fi ;
