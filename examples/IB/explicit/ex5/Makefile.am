## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp
EXTRA_DIST  = input2d
EXTRA_DIST += input3d
EXTRA_DIST += points2d.vertex
EXTRA_DIST += points3d.vertex

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d
endif
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = $(SOURCES)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input{2,3}d $(PWD) ; \
	  cp -f $(srcdir)/points{2,3}d.vertex $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input{2,3}d ; \
	  rm -f $(builddir)/points{2,3}d.vertex ; \
	fi ;
