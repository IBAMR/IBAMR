include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = elem_hmax_01 elem_hmax_02 mpi_type_wrappers

elem_hmax_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
elem_hmax_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_01_SOURCES = elem_hmax_01.cpp

elem_hmax_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
elem_hmax_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_02_SOURCES = elem_hmax_02.cpp

mpi_type_wrappers_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
mpi_type_wrappers_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
mpi_type_wrappers_SOURCES = mpi_type_wrappers.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
