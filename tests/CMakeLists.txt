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

ADD_CUSTOM_TARGET(tests)

# We use the list of test directories in two ways:
# 1. each test is added to a target tests-dir so that, e.g., 'make tests-IBFE'
#    only compiles the IBFE tests.
# 2. all input and output files in these source directories are symlinked into
#    their corresponding build directories.
SET(TEST_DIRECTORIES CIB IB IBFE IBTK adv_diff advect complex_fluids interpolate
  level_set multiphase_flow navier_stokes physical_boundary refine spread
  vc_navier_stokes wave_tank)

FOREACH(_dir ${TEST_DIRECTORIES})
  ADD_CUSTOM_TARGET("tests-${_dir}")
  ADD_DEPENDENCIES(tests "tests-${_dir}")
ENDFOREACH()

# TODO - we should clean up these macros so we don't have 3 macros that do
# basically the same thing

# Convenience macro that sets up an executable target which links against
# _target_link. For example, if the inputs are Foo, bar.cpp, and quux then we
# create a target Foo_bar in directory Foo which depends on quux.
MACRO(SETUP _dir _src _target_link)
  GET_FILENAME_COMPONENT(_out_name "${_src}" NAME_WE)
  SET(_target "tests-${_dir}_${_out_name}")
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL "${_dir}/${_src}")
  SET_TARGET_PROPERTIES(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
    ${_out_name}
    )
  TARGET_COMPILE_OPTIONS(${_target} PUBLIC -DSOURCE_DIR="${CMAKE_SOURCE_DIR}/tests/${_dir}/")
  TARGET_LINK_LIBRARIES(${_target} PRIVATE "${_target_link}")
  ADD_DEPENDENCIES("tests-${_dir}" ${_target})
ENDMACRO()

# Like the last one, but for 2d targets. For example, if the inputs are Foo and
# bar.cpp then we create a target Foo_bar_2d in directory Foo which depends on
# IBAMR2d.
MACRO(SETUP_2D _dir _src)
  GET_FILENAME_COMPONENT(_dest "${_src}" NAME_WE)
  SET(_out_name "${_dest}_2d")
  SET(_target "tests-${_dir}_${_out_name}")
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL "${_dir}/${_src}")
  SET_TARGET_PROPERTIES(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
    ${_out_name}
    )
  TARGET_COMPILE_OPTIONS(${_target} PUBLIC -DSOURCE_DIR="${CMAKE_SOURCE_DIR}/tests/${_dir}/")
  TARGET_LINK_LIBRARIES(${_target} PRIVATE IBAMR2d)
  ADD_DEPENDENCIES("tests-${_dir}" ${_target})
ENDMACRO()

# Like the last one, but for 3d targets. For example, if the inputs are Foo and
# bar.cpp then we create a target Foo_bar_3d in directory Foo which depends on
# IBAMR3d.
MACRO(SETUP_3D _dir _src)
  GET_FILENAME_COMPONENT(_dest "${_src}" NAME_WE)
  SET(_out_name "${_dest}_3d")
  SET(_target "tests-${_dir}_${_out_name}")
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL "${_dir}/${_src}")
  SET_TARGET_PROPERTIES(${_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/tests/${_dir}"
    OUTPUT_NAME
    ${_out_name}
    )
  TARGET_COMPILE_OPTIONS(${_target} PUBLIC -DSOURCE_DIR="${CMAKE_SOURCE_DIR}/tests/${_dir}/")
  TARGET_LINK_LIBRARIES(${_target} PRIVATE IBAMR3d)
  ADD_DEPENDENCIES(tests-${_dir} ${_target})
ENDMACRO()

# adv_diff:
SETUP_2D(adv_diff adv_diff_02.cpp)
SETUP_2D(adv_diff adv_diff_03.cpp)
SETUP_2D(adv_diff adv_diff_convec_opers.cpp)

SETUP_3D(adv_diff adv_diff_01.cpp)
SETUP_3D(adv_diff adv_diff_02.cpp)
SETUP_3D(adv_diff adv_diff_convec_opers.cpp)

# advect:
SETUP_2D(advect advect_01.cpp)
SETUP_3D(advect advect_01.cpp)

# complex_fluids:
SETUP(complex_fluids cf_four_roll_mill.cpp IBAMR2d)
SETUP_2D(complex_fluids cf_relaxation_op_01.cpp)
SETUP_2D(complex_fluids cf_forcing_op_01.cpp)

SETUP_3D(complex_fluids cf_relaxation_op_01.cpp)
SETUP_3D(complex_fluids cf_forcing_op_01.cpp)

# interpolate:
SETUP_2D(interpolate interpolate_01.cpp)
SETUP_3D(interpolate interpolate_01.cpp)

# level_set:
SETUP_2D(level_set fe_surface_distance.cpp)
SETUP_3D(level_set fe_surface_distance.cpp)

# multiphase_flow:
SETUP(multiphase_flow free_falling_cyl_cib.cpp IBAMR2d)
SETUP(multiphase_flow rotating_barge_cib.cpp IBAMR2d)
SETUP(multiphase_flow check_hydro_force.cpp IBAMR2d)
SETUP(multiphase_flow high_density_droplet.cpp IBAMR2d)

# navier_stokes:
SETUP_2D(navier_stokes navier_stokes_01.cpp)
SETUP_3D(navier_stokes navier_stokes_01.cpp)

# physical_boundary:
SETUP(physical_boundary extrapolation_01.cpp IBAMR2d)
SETUP_2D(physical_boundary 01.cpp)
SETUP_3D(physical_boundary 01.cpp)

# refine:
SETUP_2D(refine rt0_refine_01.cpp)
SETUP_3D(refine rt0_refine_01.cpp)

# spread:
IF(IBAMR_HAVE_LIBMESH)
  SETUP_2D(spread spread_01.cpp)
  SETUP_3D(spread spread_01.cpp)
  SETUP_2D(spread spread_02.cpp)
  SETUP_3D(spread spread_02.cpp)
ENDIF()

# vc_navier_stokes:
SETUP_2D(vc_navier_stokes vc_navier_stokes_01.cpp)
SETUP_3D(vc_navier_stokes vc_navier_stokes_01.cpp)

# wave_tank:
SETUP(wave_tank nwt_cylinder.cpp IBAMR2d)
SETUP(wave_tank nwt.cpp IBAMR2d)

# CIB:
SETUP(CIB cib_plate.cpp IBAMR2d)
SETUP(CIB cib_double_shell.cpp IBAMR3d)

# IB:
SETUP(IB explicit_ex0 IBAMR2d)
SETUP(IB explicit_ex1 IBAMR2d)

# IBFE:
IF(IBAMR_HAVE_LIBMESH)
  SETUP(IBFE interpolate_velocity_02.cpp IBAMR2d)

  SETUP_2D(IBFE explicit_ex0.cpp)
  SETUP_2D(IBFE explicit_ex1.cpp)
  SETUP_2D(IBFE explicit_ex4.cpp)
  SETUP_2D(IBFE explicit_ex5.cpp)
  SETUP_2D(IBFE explicit_ex8.cpp)
  SETUP_2D(IBFE interpolate_velocity_01.cpp)
  SETUP_2D(IBFE ib_partitioning_01.cpp)
  SETUP_2D(IBFE ib_partitioning_02.cpp)
  SETUP_2D(IBFE zero_exterior_values.cpp)

  SETUP_3D(IBFE explicit_ex2.cpp)
  SETUP_3D(IBFE explicit_ex4.cpp)
  SETUP_3D(IBFE explicit_ex5.cpp)
  SETUP_3D(IBFE ib_partitioning_01.cpp)
  SETUP_3D(IBFE ib_partitioning_02.cpp)
  SETUP_3D(IBFE interpolate_velocity_01.cpp)
  SETUP_3D(IBFE zero_exterior_values.cpp)
ENDIF()

# IBTK:
SETUP(IBTK equal_eps IBAMR2d)
SETUP(IBTK hierarchy_callbacks IBAMR2d)
SETUP(IBTK ibtk_init.cpp IBAMR2d)
SETUP(IBTK ibtk_mpi.cpp IBAMR2d)
SETUP(IBTK ldata_01.cpp IBAMR2d)
SETUP(IBTK mpi_type_wrappers.cpp IBAMR2d)

IF(IBAMR_HAVE_LIBMESH)
  SETUP(IBTK elem_hmax_01.cpp IBAMR2d)
  SETUP(IBTK elem_hmax_02.cpp IBAMR3d)
  SETUP(IBTK fe_values_01.cpp IBAMR2d)
  SETUP(IBTK fe_values_02.cpp IBAMR2d)
  SETUP(IBTK fischer_guess_01.cpp IBAMR2d)
  SETUP(IBTK jacobian_calc_01.cpp IBAMR2d)
  SETUP(IBTK mapping_01.cpp IBAMR2d)
  SETUP(IBTK subdomain_level_translation_01.cpp IBAMR2d)

  SETUP_2D(IBTK bounding_boxes_01.cpp)
  SETUP_2D(IBTK multilevel_fe_01.cpp)
ENDIF()
SETUP_2D(IBTK box_utilities_01.cpp)
SETUP_2D(IBTK ghost_accumulation_01.cpp)
SETUP_2D(IBTK ghost_indices_01.cpp)
SETUP_2D(IBTK laplace_01.cpp)
SETUP_2D(IBTK laplace_02.cpp)
SETUP_2D(IBTK laplace_03.cpp)
SETUP_2D(IBTK phys_boundary_ops.cpp)
SETUP_2D(IBTK poisson_01.cpp)
SETUP_2D(IBTK prolongation_mat.cpp)
SETUP_2D(IBTK samraidatacache_01.cpp)
SETUP_2D(IBTK vc_viscous_solver.cpp)

IF(IBAMR_HAVE_LIBMESH)
  SETUP_3D(IBTK bounding_boxes_01.cpp)
  SETUP_3D(IBTK multilevel_fe_01.cpp)
ENDIF()
SETUP_3D(IBTK box_utilities_01.cpp)
SETUP_3D(IBTK ghost_accumulation_01.cpp)
SETUP_3D(IBTK ghost_indices_01.cpp)
SETUP_3D(IBTK laplace_01.cpp)
SETUP_3D(IBTK laplace_02.cpp)
SETUP_3D(IBTK laplace_03.cpp)
SETUP_3D(IBTK phys_boundary_ops.cpp)
SETUP_3D(IBTK poisson_01.cpp)
SETUP_3D(IBTK prolongation_mat.cpp)
SETUP_3D(IBTK samraidatacache_01.cpp)
SETUP_3D(IBTK vc_viscous_solver.cpp)

ADD_CUSTOM_COMMAND(TARGET tests
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/attest ${CMAKE_BINARY_DIR}/attest)

# Set up the input and output files. Since the input and output files aren't
# really used by the build system we use a shell script to find them every time
# 'make test' is run rather than evaluating the glob when cmake generates the
# build system.
SET(TEST_DIRECTORIES CIB IB IBFE IBTK adv_diff advect complex_fluids interpolate
  level_set multiphase_flow navier_stokes physical_boundary refine spread
  vc_navier_stokes wave_tank)

FOREACH(_dir ${TEST_DIRECTORIES})
  ADD_CUSTOM_COMMAND(TARGET "tests-${_dir}"
    POST_BUILD
    COMMAND bash ${CMAKE_SOURCE_DIR}/tests/link-test-files.sh
    ${CMAKE_SOURCE_DIR}/tests/${_dir} ${CMAKE_BINARY_DIR}/tests/${_dir}
    VERBATIM)
ENDFOREACH()

# Find numdiff, if possible (we only need it for tests so its not essential that
# we find it now)
FIND_PROGRAM(NUMDIFF_EXECUTABLE NAMES numdiff HINTS ${NUMDIFF_ROOT} PATH_SUFFIXES bin)

IF ("${NUMDIFF_EXECUTABLE}" STREQUAL "NUMDIFF_EXECUTABLE-NOTFOUND")
  MESSAGE(WARNING "\
The configuration script was not able to locate numdiff. If you want to run \
the test suite you will need to either edit attest.conf, specify the path to \
numdiff to attest, or rerun CMake with the argument NUMDIFF_ROOT specifying \
numdiff's root installation directory.")
  # clear the value so that attest.conf doesn't contain an invalid path
  SET(NUMDIFF_EXECUTABLE "")
ENDIF()

# Set up the default attest configuration file:
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/tests/attest.conf.in
  ${CMAKE_BINARY_DIR}/attest.conf)
