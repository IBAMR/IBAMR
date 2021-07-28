# 2-D Tests
ibamr_add_2d_tests(fe_surface_distance
  CASES
    "<default>"
    "mpirun=8"
)

# 3-D Tests
ibamr_add_3d_tests(fe_surface_distance
  CASES
    "<default>"
    "mpirun=8"
)
