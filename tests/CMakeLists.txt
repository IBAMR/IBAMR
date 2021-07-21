## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# adv_diff:
ibamr_add_2d_tests(adv_diff_02
  DIRECTORY adv_diff)
ibamr_add_2d_tests(adv_diff_03
  DIRECTORY adv_diff)
ibamr_add_2d_tests(adv_diff_convec_opers
  DIRECTORY adv_diff)
ibamr_add_2d_tests(bp_adv_diff_01
  DIRECTORY adv_diff)
ibamr_add_2d_tests(bp_adv_diff_02
  DIRECTORY adv_diff)

ibamr_add_3d_tests(adv_diff_01
  DIRECTORY adv_diff)
ibamr_add_3d_tests(adv_diff_02
  DIRECTORY adv_diff)
ibamr_add_3d_tests(adv_diff_convec_opers
  DIRECTORY adv_diff )

# advect:
ibamr_add_2d_tests(advect_01
  DIRECTORY advect)
ibamr_add_3d_tests(advect_01
  DIRECTORY advect)

#coarsen:
ibamr_add_2d_tests(rt0_coarsen_01 DIRECTORY coarsen)
ibamr_add_3d_tests(rt0_coarsen_01 DIRECTORY coarsen)

# complex_fluids:
ibamr_add_tests(cf_four_roll_mill DIRECTORY complex_fluids LIBRARIES IBAMR2d)
ibamr_add_2d_tests(cf_relaxation_op_01 DIRECTORY complex_fluids)
ibamr_add_2d_tests(cf_forcing_op_01 DIRECTORY complex_fluids)

ibamr_add_3d_tests(cf_relaxation_op_01 DIRECTORY complex_fluids)
ibamr_add_3d_tests(cf_forcing_op_01 DIRECTORY complex_fluids)

# fe_mechanics:
if(IBAMR_HAVE_LIBMESH)
  ibamr_add_2d_tests(fe_mechanics_ex0 DIRECTORY fe_mechanics)
  ibamr_add_3d_tests(fe_mechanics_ex0 DIRECTORY fe_mechanics)
endif()

# interpolate:
ibamr_add_2d_tests(interpolate_01 DIRECTORY interpolate)
ibamr_add_3d_tests(interpolate_01 DIRECTORY interpolate)

# level_set:
ibamr_add_2d_tests(fe_surface_distance DIRECTORY level_set)
ibamr_add_3d_tests(fe_surface_distance DIRECTORY level_set)

# multiphase_flow:
ibamr_add_tests(free_falling_cyl_cib DIRECTORY multiphase_flow LIBRARIES IBAMR2d)
ibamr_add_tests(rotating_barge_cib DIRECTORY multiphase_flow LIBRARIES IBAMR2d)
ibamr_add_tests(check_hydro_force DIRECTORY multiphase_flow LIBRARIES IBAMR2d)
ibamr_add_tests(high_density_droplet DIRECTORY multiphase_flow LIBRARIES IBAMR2d)

# navier_stokes:
ibamr_add_2d_tests(navier_stokes_01 DIRECTORY navier_stokes)
ibamr_add_3d_tests(navier_stokes_01 DIRECTORY navier_stokes)

# physical_boundary:
ibamr_add_tests(extrapolation_01 DIRECTORY physical_boundary LIBRARIES IBAMR2d)
ibamr_add_2d_tests(01 DIRECTORY physical_boundary)
ibamr_add_3d_tests(01 DIRECTORY physical_boundary)

# refine:
ibamr_add_2d_tests(rt0_refine_01 DIRECTORY refine)
ibamr_add_3d_tests(rt0_refine_01 DIRECTORY refine)

# spread:
if(IBAMR_HAVE_LIBMESH)
  ibamr_add_2d_tests(spread_01 DIRECTORY spread)
  ibamr_add_3d_tests(spread_01 DIRECTORY spread)
  ibamr_add_2d_tests(spread_02 DIRECTORY spread)
  ibamr_add_3d_tests(spread_02 DIRECTORY spread)
endif()

# vc_navier_stokes:
ibamr_add_2d_tests(vc_navier_stokes_01 DIRECTORY vc_navier_stokes)
ibamr_add_3d_tests(vc_navier_stokes_01 DIRECTORY vc_navier_stokes)

# wave_tank:
IF (IBAMR_HAVE_LIBMESH)
ibamr_add_tests(nwt_cylinder DIRECTORY wave_tank LIBRARIES IBAMR2d)
ibamr_add_tests(nwt DIRECTORY wave_tank LIBRARIES IBAMR2d)
ENDIF ()

# CIB:
ibamr_add_tests(cib_plate DIRECTORY CIB LIBRARIES IBAMR2d)
ibamr_add_tests(cib_double_shell DIRECTORY CIB LIBRARIES IBAMR3d)

# ConstraintIB:
ibamr_add_tests(oscillating_rigid_cylinder DIRECTORY ConstraintIB LIBRARIES IBAMR2d)

# IB:
ibamr_add_tests(explicit_ex0 DIRECTORY IB LIBRARIES IBAMR2d)
ibamr_add_tests(explicit_ex1 DIRECTORY IB LIBRARIES IBAMR2d)
ibamr_add_tests(ib_body_force DIRECTORY IB LIBRARIES IBAMR2d)
ibamr_add_tests(ib_body_force_kirchhoff DIRECTORY IB LIBRARIES IBAMR3d)

# IBFE:
if(IBAMR_HAVE_LIBMESH)
  ibamr_add_tests(interpolate_velocity_02 DIRECTORY IBFE LIBRARIES IBAMR2d)

  ibamr_add_2d_tests(explicit_ex0 DIRECTORY IBFE)
  ibamr_add_2d_tests(explicit_ex1 DIRECTORY IBFE)
  ibamr_add_2d_tests(explicit_ex4 DIRECTORY IBFE)
  ibamr_add_2d_tests(explicit_ex5 DIRECTORY IBFE)
  ibamr_add_2d_tests(explicit_ex8 DIRECTORY IBFE)
  ibamr_add_2d_tests(interpolate_velocity_01 DIRECTORY IBFE)
  ibamr_add_2d_tests(ib_partitioning_01 DIRECTORY IBFE)
  ibamr_add_2d_tests(ib_partitioning_02 DIRECTORY IBFE)
  ibamr_add_2d_tests(zero_exterior_values DIRECTORY IBFE)

  ibamr_add_3d_tests(explicit_ex2 DIRECTORY IBFE)
  ibamr_add_3d_tests(explicit_ex4 DIRECTORY IBFE)
  ibamr_add_3d_tests(explicit_ex5 DIRECTORY IBFE)
  ibamr_add_3d_tests(ib_partitioning_01 DIRECTORY IBFE)
  ibamr_add_3d_tests(ib_partitioning_02 DIRECTORY IBFE)
  ibamr_add_3d_tests(interpolate_velocity_01 DIRECTORY IBFE)
  ibamr_add_3d_tests(zero_exterior_values DIRECTORY IBFE)
endif()

# IBTK:
ibamr_add_tests(equal_eps DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(hierarchy_callbacks DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(ibtk_init DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(ibtk_mpi DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(ldata_01 DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(mpi_type_wrappers DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(child_integrators DIRECTORY IBTK LIBRARIES IBAMR2d)
ibamr_add_tests(version_macros DIRECTORY IBTK  LIBRARIES IBAMR2d)

if (IBAMR_HAVE_LIBMESH)
  ibamr_add_tests(elem_hmax_01 DIRECTORY IBTK LIBRARIES IBAMR2d)
  ibamr_add_tests(elem_hmax_02 DIRECTORY IBTK LIBRARIES IBAMR3d)
  ibamr_add_tests(fe_values_01 DIRECTORY IBTK LIBRARIES IBAMR2d)
  ibamr_add_tests(fe_values_02 DIRECTORY IBTK LIBRARIES IBAMR2d)
  ibamr_add_tests(fischer_guess_01 DIRECTORY IBTK LIBRARIES IBAMR2d)
  ibamr_add_tests(jacobian_calc_01 DIRECTORY IBTK LIBRARIES IBAMR2d)
  ibamr_add_tests(mapping_01 DIRECTORY IBTK LIBRARIES IBAMR2d)
  ibamr_add_tests(subdomain_level_translation_01 DIRECTORY IBTK LIBRARIES IBAMR2d)

  ibamr_add_2d_tests(bounding_boxes_01 DIRECTORY IBTK)
  ibamr_add_2d_tests(multilevel_fe_01 DIRECTORY IBTK)
endif ()
ibamr_add_2d_tests(box_utilities_01 DIRECTORY IBTK)
ibamr_add_2d_tests(ghost_accumulation_01 DIRECTORY IBTK)
ibamr_add_2d_tests(ghost_indices_01 DIRECTORY IBTK)
ibamr_add_2d_tests(laplace_01 DIRECTORY IBTK)
ibamr_add_2d_tests(laplace_02 DIRECTORY IBTK)
ibamr_add_2d_tests(laplace_03 DIRECTORY IBTK)
ibamr_add_2d_tests(phys_boundary_ops DIRECTORY IBTK)
ibamr_add_2d_tests(poisson_01 DIRECTORY IBTK)
ibamr_add_2d_tests(prolongation_mat DIRECTORY IBTK)
ibamr_add_2d_tests(samraidatacache_01 DIRECTORY IBTK)
ibamr_add_2d_tests(secondary_hierarchy_01 DIRECTORY IBTK)
ibamr_add_2d_tests(vc_viscous_solver DIRECTORY IBTK)
ibamr_add_2d_tests(helmholtz DIRECTORY IBTK)

if(IBAMR_HAVE_LIBMESH)
  ibamr_add_3d_tests(bounding_boxes_01 DIRECTORY IBTK)
  ibamr_add_3d_tests(multilevel_fe_01 DIRECTORY IBTK)
endif()
ibamr_add_3d_tests(box_utilities_01 DIRECTORY IBTK)
ibamr_add_3d_tests(ghost_accumulation_01 DIRECTORY IBTK)
ibamr_add_3d_tests(ghost_indices_01 DIRECTORY IBTK)
ibamr_add_3d_tests(laplace_01 DIRECTORY IBTK)
ibamr_add_3d_tests(laplace_02 DIRECTORY IBTK)
ibamr_add_3d_tests(laplace_03 DIRECTORY IBTK)
ibamr_add_3d_tests(phys_boundary_ops DIRECTORY IBTK)
ibamr_add_3d_tests(poisson_01 DIRECTORY IBTK)
ibamr_add_3d_tests(prolongation_mat DIRECTORY IBTK)
ibamr_add_3d_tests(samraidatacache_01 DIRECTORY IBTK)
ibamr_add_3d_tests(vc_viscous_solver DIRECTORY IBTK)
ibamr_add_3d_tests(helmholtz DIRECTORY IBTK)
