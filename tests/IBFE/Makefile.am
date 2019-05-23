include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = interpolate_velocity_01 interpolate_velocity_02 explicit_ex4_2d \
explicit_ex4_3d

interpolate_velocity_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
interpolate_velocity_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
interpolate_velocity_01_SOURCES = interpolate_velocity_01.cpp

interpolate_velocity_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
interpolate_velocity_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
interpolate_velocity_02_SOURCES = interpolate_velocity_02.cpp

explicit_ex4_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
explicit_ex4_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
explicit_ex4_2d_SOURCES = explicit_ex4.cpp

explicit_ex4_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
explicit_ex4_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
explicit_ex4_3d_SOURCES = explicit_ex4.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
