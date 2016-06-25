## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp generate_membrane3d.cpp
EXTRA_DIST  = input3d read_file_ibamr3d.m
EXTRA_DIST += generate_membrane3d.cpp

EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d membrane3d_generator
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = main.cpp

membrane3d_generator_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
membrane3d_generator_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
membrane3d_generator_SOURCES = generate_membrane3d.cpp

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/generate_membrane3d.cpp $(PWD) ; \
	  cp -f $(srcdir)/read_file_ibamr3d.m $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/generate_membrane3d.cpp ; \
	  rm -f $(builddir)/read_file_ibamr3d.m ; \
	fi ;
