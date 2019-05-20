include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = free_falling_cyl_cib rotating_barge_cib

# this test needs some extra input files, so make SOURCE_DIR available:
free_falling_cyl_cib_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
free_falling_cyl_cib_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
free_falling_cyl_cib_SOURCES = free_falling_cyl_cib.cpp

rotating_barge_cib_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
rotating_barge_cib_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
rotating_barge_cib_SOURCES = rotating_barge_cib.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
