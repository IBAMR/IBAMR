# 2-D Tests
ibamr_add_2d_tests(adv_diff_02)
ibamr_add_2d_tests(adv_diff_03)
ibamr_add_2d_tests(adv_diff_convec_opers
  CASES
    "constant"
    "linear"
    "none"
    "quadratic"
)
ibamr_add_2d_tests(bp_adv_diff_01)
ibamr_add_2d_tests(bp_adv_diff_02
  CASES
    "mpirun=8"
)

# 3-D Tests
ibamr_add_3d_tests(adv_diff_01)
ibamr_add_3d_tests(adv_diff_02
  CASES
    "mpirun=2"
    "mpirun=8"
)
ibamr_add_3d_tests(adv_diff_convec_opers
  CASES
    "constant"
    "linear"
    "none"
    "quadratic"
)
