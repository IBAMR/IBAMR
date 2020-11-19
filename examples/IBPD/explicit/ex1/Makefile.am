## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp generate_sheet3d.cpp
EXTRA_DIST  = input3d
EXTRA_DIST += generate_sheet3d.cpp 

EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d sheet3d_generator
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = main.cpp

sheet3d_generator_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
sheet3d_generator_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
sheet3d_generator_SOURCES = generate_sheet3d.cpp

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/generate_sheet3d.cpp $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/generate_sheet3d.cpp ; \
	fi ;
