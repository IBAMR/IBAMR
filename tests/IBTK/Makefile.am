## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = mpi_type_wrappers poisson_01_2d \
poisson_01_3d samraidatacache_01_2d samraidatacache_01_3d laplace_01_2d \
laplace_01_3d laplace_02_2d laplace_02_3d laplace_03_2d laplace_03_3d ldata_01 \
prolongation_mat_2d prolongation_mat_3d phys_boundary_ops_2d phys_boundary_ops_3d \
vc_viscous_solver_2d vc_viscous_solver_3d box_utilities_01_2d box_utilities_01_3d \
ghost_accumulation_01_2d ghost_accumulation_01_3d ghost_indices_01_2d \
ghost_indices_01_3d ibtk_init hierarchy_callbacks ibtk_mpi equal_eps helmholtz_2d \
helmholtz_3d secondary_hierarchy_01_2d child_integrators_2d version_macros

if LIBMESH_ENABLED
EXTRA_PROGRAMS += elem_hmax_01 elem_hmax_02 jacobian_calc_01 bounding_boxes_01_2d \
bounding_boxes_01_3d mapping_01 fe_values_01 fe_values_02 \
multilevel_fe_01_2d multilevel_fe_01_3d subdomain_level_translation_01 \
fischer_guess_01
endif

if LIBMESH_ENABLED
elem_hmax_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
elem_hmax_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
elem_hmax_01_SOURCES = elem_hmax_01.cpp
endif

if LIBMESH_ENABLED
elem_hmax_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
elem_hmax_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
elem_hmax_02_SOURCES = elem_hmax_02.cpp
endif

equal_eps_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
equal_eps_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
equal_eps_SOURCES = equal_eps.cpp

ibtk_init_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ibtk_init_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ibtk_init_SOURCES = ibtk_init.cpp

ibtk_mpi_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ibtk_mpi_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ibtk_mpi_SOURCES = ibtk_mpi.cpp

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

if LIBMESH_ENABLED
jacobian_calc_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
jacobian_calc_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
jacobian_calc_01_SOURCES = jacobian_calc_01.cpp
endif

if LIBMESH_ENABLED
bounding_boxes_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
bounding_boxes_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
bounding_boxes_01_2d_SOURCES = bounding_boxes_01.cpp
endif

if LIBMESH_ENABLED
bounding_boxes_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
bounding_boxes_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
bounding_boxes_01_3d_SOURCES = bounding_boxes_01.cpp
endif

if LIBMESH_ENABLED
mapping_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
mapping_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
mapping_01_SOURCES = mapping_01.cpp
endif

if LIBMESH_ENABLED
# this test reads meshes stored in the present directory so make SOURCE_DIR available:
fe_values_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
fe_values_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
fe_values_01_SOURCES = fe_values_01.cpp
endif

if LIBMESH_ENABLED
# this test reads meshes stored in the present directory so make SOURCE_DIR available:
fe_values_02_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
fe_values_02_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
fe_values_02_SOURCES = fe_values_02.cpp
endif

box_utilities_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
box_utilities_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
box_utilities_01_2d_SOURCES = box_utilities_01.cpp

box_utilities_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
box_utilities_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
box_utilities_01_3d_SOURCES = box_utilities_01.cpp

ghost_indices_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ghost_indices_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ghost_indices_01_2d_SOURCES = ghost_indices_01.cpp

ghost_indices_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
ghost_indices_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
ghost_indices_01_3d_SOURCES = ghost_indices_01.cpp

ghost_accumulation_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
ghost_accumulation_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
ghost_accumulation_01_2d_SOURCES = ghost_accumulation_01.cpp

ghost_accumulation_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
ghost_accumulation_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
ghost_accumulation_01_3d_SOURCES = ghost_accumulation_01.cpp

hierarchy_callbacks_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
hierarchy_callbacks_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
hierarchy_callbacks_SOURCES = hierarchy_callbacks.cpp

if LIBMESH_ENABLED
subdomain_level_translation_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
subdomain_level_translation_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
subdomain_level_translation_01_SOURCES = subdomain_level_translation_01.cpp
endif

if LIBMESH_ENABLED
multilevel_fe_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
multilevel_fe_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
multilevel_fe_01_2d_SOURCES = multilevel_fe_01.cpp
multilevel_fe_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3 -DSOURCE_DIR=\"$(abs_srcdir)\"
multilevel_fe_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
multilevel_fe_01_3d_SOURCES = multilevel_fe_01.cpp
endif

if LIBMESH_ENABLED
fischer_guess_01_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2 -DSOURCE_DIR=\"$(abs_srcdir)\"
fischer_guess_01_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
fischer_guess_01_SOURCES = fischer_guess_01.cpp
endif

helmholtz_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
helmholtz_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
helmholtz_2d_SOURCES = helmholtz.cpp

helmholtz_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
helmholtz_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
helmholtz_3d_SOURCES = helmholtz.cpp

secondary_hierarchy_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
secondary_hierarchy_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
secondary_hierarchy_01_2d_SOURCES = secondary_hierarchy_01.cpp

child_integrators_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
child_integrators_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
child_integrators_2d_SOURCES = child_integrators.cpp

version_macros_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
version_macros_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
version_macros_SOURCES = version_macros.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
