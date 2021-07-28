# 2-D Tests
ibamr_add_tests(nwt_cylinder LIBRARIES IBAMR2d
  CASES
    "<default>"
    "mpirun=8"
)
ibamr_add_tests(nwt LIBRARIES IBAMR2d
  CASES
    "mpirun=8"
)

# 3-D Tests
