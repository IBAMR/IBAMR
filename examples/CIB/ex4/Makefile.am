## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp
EXTRA_DIST  = input3d.sphere sphere_3d.vertex 

EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d.sphere $(PWD) ; \
	  cp -f $(srcdir)/sphere_3d.vertex $(PWD) ; \
	  cp -f $(srcdir)/check_sol.m $(PWD) ; \
	  cp -f $(srcdir)/PetscOptions.dat $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d.sphere ; \
          rm -f $(builddir)/sphere_3d.vertex ; \
	  rm -f $(builddir)/check_sol.m ; \
	  rm -f $(builddir)/PetscOptions.dat ; \
	fi ;
