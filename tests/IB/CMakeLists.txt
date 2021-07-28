# 2-D Tests
ibamr_add_tests(explicit_ex0
  LIBRARIES IBAMR2d
  CASES
    "<default>"
    "restart=300.mpirun=4"
)
ibamr_add_tests(explicit_ex1
  LIBRARIES IBAMR2d
)
ibamr_add_tests(ib_body_force
  LIBRARIES IBAMR2d
)

# 3-D Tests
ibamr_add_tests(ib_body_force_kirchhoff
  LIBRARIES IBAMR3d
)
