# 2-D Tests
ibamr_add_tests(equal_eps
  LIBRARIES IBAMR2d
)
ibamr_add_tests(hierarchy_callbacks
  LIBRARIES IBAMR2d
)
ibamr_add_tests(ibtk_init
  LIBRARIES IBAMR2d
)
ibamr_add_tests(ibtk_mpi
  LIBRARIES IBAMR2d
  CASES
    "mpirun=8"
)
ibamr_add_tests(ldata_01
  LIBRARIES IBAMR2d
  CASES
    "<default>"
    "mpirun=4"
)
ibamr_add_tests(mpi_type_wrappers
  LIBRARIES IBAMR2d
)

ibamr_add_tests(child_integrators
  LIBRARIES IBAMR2d)

ibamr_add_tests(version_macros
  LIBRARIES IBAMR2d)

if (IBAMR_HAVE_LIBMESH)
  ibamr_add_tests(elem_hmax_01
    LIBRARIES IBAMR2d
  )
  ibamr_add_tests(fe_values_01
    LIBRARIES IBAMR2d
  )
  ibamr_add_tests(fe_values_02
    LIBRARIES IBAMR2d
  )
  ibamr_add_tests(fischer_guess_01
    LIBRARIES IBAMR2d
  )
  ibamr_add_tests(jacobian_calc_01
    LIBRARIES IBAMR2d
  )
  ibamr_add_tests(mapping_01
    LIBRARIES IBAMR2d
  )
  ibamr_add_tests(subdomain_level_translation_01
    LIBRARIES IBAMR2d
    CASES
      "<default>"
      "expect_error=true"
  )

  ibamr_add_2d_tests(bounding_boxes_01
    CASES
      "<default>"
      "amr"
  )
  ibamr_add_2d_tests(multilevel_fe_01
    CASES
      "<default>"
      "half"
  )
endif ()
ibamr_add_2d_tests(box_utilities_01)
ibamr_add_2d_tests(ghost_accumulation_01
  CASES
    "cell.a"
    "cell.a.mpirun=2"
    "cell.spread.a"
    "cell.spread.b"
    "cell.spread.b.mpirun=2"
    "cell.spread.b.mpirun=4"
    "side.b"
    "side"
    "side.mpirun=2"
    "side.spread.a"
    "side.spread.b"
    "side.spread.b.mpirun=2"
    "side.spread.b.mpirun=4"
)
ibamr_add_2d_tests(ghost_indices_01
  CASES
    "a.mpirun=2"
    "b.mpirun=2"
)
ibamr_add_2d_tests(laplace_01
  CASES
    "mpirun=6"
)
ibamr_add_2d_tests(laplace_02)
ibamr_add_2d_tests(laplace_03)
ibamr_add_2d_tests(phys_boundary_ops)
ibamr_add_2d_tests(poisson_01
  CASES
    "mpirun=4"
)
ibamr_add_2d_tests(prolongation_mat)
ibamr_add_2d_tests(samraidatacache_01)
ibamr_add_2d_tests(secondary_hierarchy_01
  CASES
    "<default>"
    "mpirun=4"
)
ibamr_add_2d_tests(vc_viscous_solver)
ibamr_add_2d_tests(helmholtz
  CASES
    "a"
    "b"
    "c"
    "d"
)

# 3-D Tests
if(IBAMR_HAVE_LIBMESH)
  ibamr_add_tests(elem_hmax_02
    LIBRARIES IBAMR3d
  )
  ibamr_add_3d_tests(bounding_boxes_01
    CASES
      "<default>"
      "amr"
  )
  ibamr_add_3d_tests(multilevel_fe_01
    CASES
      "<default>"
      "quarter"
      "quarter.mpirun=4"
  )
endif()
ibamr_add_3d_tests(box_utilities_01)
ibamr_add_3d_tests(ghost_accumulation_01
  CASES
    "cell.spread.a"
    "cell.spread.b"
    "cell.spread.b.mpirun=2"
    "cell.spread.b.mpirun=4"
    "side.spread.a"
    "side.spread.b"
    "side.spread.b.mpirun=2"
    "side.spread.b.mpirun=4"
)
ibamr_add_3d_tests(ghost_indices_01
  CASES
    "a.mpirun=4"
    "b.mpirun=4"
)
ibamr_add_3d_tests(laplace_01
  CASES
    "mpirun=8"
)
ibamr_add_3d_tests(laplace_02)
ibamr_add_3d_tests(laplace_03)
ibamr_add_3d_tests(phys_boundary_ops)
ibamr_add_3d_tests(poisson_01
  CASES
    "mpirun=4"
)
ibamr_add_3d_tests(prolongation_mat)
ibamr_add_3d_tests(samraidatacache_01)
ibamr_add_3d_tests(vc_viscous_solver)
ibamr_add_3d_tests(helmholtz
  CASES
    "a.mpirun=2"
    "b.mpirun=2"
    "c.mpirun=4"
    "d.mpirun=4"
)
