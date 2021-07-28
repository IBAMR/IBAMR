# 2-D Tests
ibamr_add_tests(cib_plate
  LIBRARIES IBAMR2d
  CASES
    "mpirun=4"
    "restart=5.mpirun=4"
)

# 3-D Tests
ibamr_add_tests(cib_double_shell
  LIBRARIES IBAMR3d
  CASES
    "<default>"
    "mpirun=4"
)
