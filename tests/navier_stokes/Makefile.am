## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2024 by the IBAMR developers
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

EXTRA_PROGRAMS = navier_stokes_01_2d navier_stokes_01_3d poiseuille_flow_2d stokes_operator_2d stokes_operator_3d \
                 stokes_petsc_mat_utilities_ops_2d stokes_petsc_mat_utilities_prolongation_2d \
                 stokes_ib_cav_sparse_map_closure_2d stokes_ib_cav_patch_level_maps_2d \
                 stokes_ib_cav_matlab_relaxed_parity_2d stokes_ib_cav_matlab_strict_parity_2d \
                 stokes_petsc_level_solver_asm_modes_2d

navier_stokes_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
navier_stokes_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
navier_stokes_01_2d_SOURCES = navier_stokes_01.cpp

navier_stokes_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
navier_stokes_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
navier_stokes_01_3d_SOURCES = navier_stokes_01.cpp

poiseuille_flow_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
poiseuille_flow_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
poiseuille_flow_2d_SOURCES = poiseuille_flow_2d.cpp

stokes_operator_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_operator_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_operator_2d_SOURCES = stokes_operator.cpp

stokes_operator_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
stokes_operator_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
stokes_operator_3d_SOURCES = stokes_operator.cpp

stokes_petsc_mat_utilities_ops_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_petsc_mat_utilities_ops_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_petsc_mat_utilities_ops_2d_SOURCES = stokes_petsc_mat_utilities_ops.cpp

stokes_petsc_mat_utilities_prolongation_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_petsc_mat_utilities_prolongation_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_petsc_mat_utilities_prolongation_2d_SOURCES = stokes_petsc_mat_utilities_prolongation.cpp

stokes_ib_cav_sparse_map_closure_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_ib_cav_sparse_map_closure_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_ib_cav_sparse_map_closure_2d_SOURCES = stokes_ib_cav_sparse_map_closure.cpp

stokes_ib_cav_patch_level_maps_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_ib_cav_patch_level_maps_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_ib_cav_patch_level_maps_2d_SOURCES = stokes_ib_cav_patch_level_maps.cpp

stokes_ib_cav_matlab_relaxed_parity_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_ib_cav_matlab_relaxed_parity_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_ib_cav_matlab_relaxed_parity_2d_SOURCES = stokes_ib_cav_matlab_relaxed_parity.cpp

stokes_ib_cav_matlab_strict_parity_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_ib_cav_matlab_strict_parity_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_ib_cav_matlab_strict_parity_2d_SOURCES = stokes_ib_cav_matlab_strict_parity.cpp

stokes_petsc_level_solver_asm_modes_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
stokes_petsc_level_solver_asm_modes_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
stokes_petsc_level_solver_asm_modes_2d_SOURCES = stokes_petsc_level_solver_asm_modes.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
