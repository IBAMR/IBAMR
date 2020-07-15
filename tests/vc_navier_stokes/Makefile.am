include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = vc_navier_stokes_01_2d vc_navier_stokes_01_3d

vc_navier_stokes_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
vc_navier_stokes_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
vc_navier_stokes_01_2d_SOURCES = vc_navier_stokes_01.cpp

vc_navier_stokes_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3 -DSOURCE_DIR=\"$(abs_srcdir)\"
vc_navier_stokes_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
vc_navier_stokes_01_3d_SOURCES = vc_navier_stokes_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
