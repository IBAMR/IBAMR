include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = elem_hmax_01 elem_hmax_02

elem_hmax_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
elem_hmax_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_01_SOURCES = elem_hmax_01.cpp

elem_hmax_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
elem_hmax_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_02_SOURCES = elem_hmax_02.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -L -f -s $(srcdir)/*input $(PWD) ; \
	  ln -L -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
