include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = elem_hmax_01 elem_hmax_02 mpi_type_wrappers poisson_01_2d \
poisson_01_3d samraidatacache_01_2d samraidatacache_01_3d laplace_01_2d \
laplace_01_3d ldata_01

elem_hmax_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
elem_hmax_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_01_SOURCES = elem_hmax_01.cpp

elem_hmax_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
elem_hmax_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_02_SOURCES = elem_hmax_02.cpp

mpi_type_wrappers_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
mpi_type_wrappers_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
mpi_type_wrappers_SOURCES = mpi_type_wrappers.cpp

laplace_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
laplace_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
laplace_01_2d_SOURCES = laplace_01.cpp

laplace_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
laplace_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
laplace_01_3d_SOURCES = laplace_01.cpp

poisson_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
poisson_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
poisson_01_2d_SOURCES = poisson_01.cpp

poisson_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
poisson_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
poisson_01_3d_SOURCES = poisson_01.cpp

samraidatacache_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
samraidatacache_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
samraidatacache_01_2d_SOURCES = samraidatacache_01.cpp

samraidatacache_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
samraidatacache_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
samraidatacache_01_3d_SOURCES = samraidatacache_01.cpp

ldata_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ldata_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ldata_01_SOURCES = ldata_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
