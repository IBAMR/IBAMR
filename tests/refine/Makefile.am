include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = rt0_refine_01_2d rt0_refine_01_3d

rt0_refine_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
rt0_refine_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
rt0_refine_01_2d_SOURCES = rt0_refine_01.cpp

rt0_refine_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
rt0_refine_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
rt0_refine_01_3d_SOURCES = rt0_refine_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
