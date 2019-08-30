include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = extrapolation_01 physical_boundary_2d_01 physical_boundary_3d_01

extrapolation_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
extrapolation_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
extrapolation_01_SOURCES = extrapolation_01.cpp

physical_boundary_2d_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
physical_boundary_2d_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
physical_boundary_2d_01_SOURCES = physical_boundary_01.cpp

physical_boundary_3d_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
physical_boundary_3d_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
physical_boundary_3d_01_SOURCES = physical_boundary_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
