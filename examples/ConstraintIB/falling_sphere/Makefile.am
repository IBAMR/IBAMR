## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

if GSL_ENABLED

## Dimension-dependent testers
EXTRA_DIST  = input3d 
EXTRA_DIST += sphereGen3d.cpp


EXTRA_PROGRAMS =
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d sphereGen3d
endif

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
main3d_SOURCES = CartGridBodyForce.cpp ForceProjector.cpp RigidBodyKinematics.cpp main.cpp

sphereGen3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
sphereGen3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
sphereGen3d_SOURCES = sphereGen3d.cpp

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input3d $(PWD) ; \
	  cp -f $(srcdir)/sphere3d.vertex $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input3d ; \
	  rm -f $(builddir)/sphere3d.vertex ;\
	fi ;

else
examples:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo " This example needs GNU Scientific Library.                          "
	@echo " Configure IBAMR with additional options --enable-gsl --with-gsl=path"
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif
