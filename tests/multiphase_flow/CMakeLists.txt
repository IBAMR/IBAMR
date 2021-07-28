# 2-D Tests
ibamr_add_tests(free_falling_cyl_cib
  LIBRARIES IBAMR2d
  CASES
    "mpirun=4"
)
ibamr_add_tests(rotating_barge_cib
  LIBRARIES IBAMR2d
  CASES
    "mpirun=4"
)
ibamr_add_tests(check_hydro_force
  LIBRARIES IBAMR2d
  CASES
    "mpirun=8"
)
ibamr_add_tests(high_density_droplet
  LIBRARIES IBAMR2d
)

# 3-D Tests
