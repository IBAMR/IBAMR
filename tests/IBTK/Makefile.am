include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = elem_hmax_01 elem_hmax_02 mpi_type_wrappers poisson_01_2d \
poisson_01_3d samraidatacache_01_2d samraidatacache_01_3d laplace_01_2d \
laplace_01_3d laplace_02_2d laplace_02_3d laplace_03_2d laplace_03_3d ldata_01 \
prolongation_mat_2d prolongation_mat_3d phys_boundary_ops_2d phys_boundary_ops_3d \
vc_viscous_solver_2d vc_viscous_solver_3d jacobian_calc_01

elem_hmax_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
elem_hmax_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_01_SOURCES = elem_hmax_01.cpp

elem_hmax_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
elem_hmax_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
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

laplace_02_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
laplace_02_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
laplace_02_2d_SOURCES = laplace_02.cpp

laplace_02_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
laplace_02_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
laplace_02_3d_SOURCES = laplace_02.cpp

laplace_03_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
laplace_03_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
laplace_03_2d_SOURCES = laplace_03.cpp

laplace_03_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
laplace_03_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
laplace_03_3d_SOURCES = laplace_03.cpp

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

prolongation_mat_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
prolongation_mat_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
prolongation_mat_2d_SOURCES = prolongation_mat.cpp

prolongation_mat_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
prolongation_mat_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
prolongation_mat_3d_SOURCES = prolongation_mat.cpp

phys_boundary_ops_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
phys_boundary_ops_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
phys_boundary_ops_2d_SOURCES = phys_boundary_ops.cpp

phys_boundary_ops_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
phys_boundary_ops_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
phys_boundary_ops_3d_SOURCES = phys_boundary_ops.cpp

vc_viscous_solver_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
vc_viscous_solver_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
vc_viscous_solver_2d_SOURCES = vc_viscous_solver.cpp

vc_viscous_solver_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
vc_viscous_solver_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
vc_viscous_solver_3d_SOURCES = vc_viscous_solver.cpp

jacobian_calc_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
jacobian_calc_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
jacobian_calc_01_SOURCES = jacobian_calc_01.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
