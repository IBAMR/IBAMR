## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.cpp
EXTRA_DIST  = input3d shell_3d_42.vertex shell_3d_162.vertex shell_3d_42.clones shell_3d_162.clones PetscOptions.dat shell_3d_42.constraint

EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/shell_3d_42.vertex $(PWD) ; \
	  cp -f $(srcdir)/shell_3d_162.vertex $(PWD) ; \
	  cp -f $(srcdir)/shell_3d_42.clones $(PWD) ; \
	  cp -f $(srcdir)/shell_3d_162.clones $(PWD) ; \
	  cp -f $(srcdir)/PetscOptions.dat $(PWD) ; \
	  cp -f $(srcdir)/shell_3d_42.constraint $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
          rm -f $(builddir)/shell_3d_42.vertex ; \
          rm -f $(builddir)/shell_3d_162.vertex ; \
          rm -f $(builddir)/shell_3d_42.clones ; \
          rm -f $(builddir)/shell_3d_162.clones ; \
          rm -f $(builddir)/PetscOptions.dat ; \
          rm -f $(builddir)/shell_3d_42.constraint ; \
	fi ;
