include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = extrapolation_01 01_2d 01_3d

extrapolation_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
extrapolation_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
extrapolation_01_SOURCES = extrapolation_01.cpp

01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
01_2d_SOURCES = physical_boundary_01.cpp

01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
01_3d_SOURCES = physical_boundary_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
