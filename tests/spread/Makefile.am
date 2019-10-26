include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = spread_01_2d spread_01_3d

spread_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
spread_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
spread_01_2d_SOURCES = spread_01.cpp

spread_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
spread_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
spread_01_3d_SOURCES = spread_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
