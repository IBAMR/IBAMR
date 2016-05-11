## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp generate_membrane2d.cpp
EXTRA_DIST  = input2d
EXTRA_DIST += generate_membrane2d.cpp 

EXTRA_DIST += membrane2d.vertex 
EXTRA_DIST += membrane2d.spring 

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d membrane2d_generator
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
main2d_SOURCES = main.cpp

membrane2d_generator_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
membrane2d_generator_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
membrane2d_generator_SOURCES = generate_membrane2d.cpp

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d $(PWD) ; \
	  cp -f $(srcdir)/membrane2d.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/generate_membrane2d.cpp $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/membrane2d.{vertex,spring} ; \
	  rm -f $(builddir)/generate_membrane2d.cpp ; \
	fi ;
