if(IBAMR_HAVE_LIBMESH)
  ibamr_add_2d_tests(fe_mechanics_ex0
    CASES
      "<default>"
      "mpirun=4"
  )
  ibamr_add_3d_tests(fe_mechanics_ex0
    CASES
      "<default>"
      "mpirun=4"
  )
endif()

