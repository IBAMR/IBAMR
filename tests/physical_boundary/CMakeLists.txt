# 2-D Tests
ibamr_add_tests(extrapolation_01
  LIBRARIES IBAMR2d
  CASES
    "cell.constant"
    "cell.linear"
    "cell.quadratic"
    "face.constant"
    "face.linear"
    "face.quadratic"
    "node.constant"
    "node.linear"
    "node.quadratic"
    "side.constant"
    "side.linear"
    "side.quadratic"
)
ibamr_add_2d_tests(01
  CASES
    "cell.dirichlet.linear"
    "cell.dirichlet.quadratic"
    "cell.neumann.linear"
    "cell.neumann.quadratic"
    "side.dirichlet.linear"
    "side.dirichlet.quadratic"
    "side.neumann.linear"
    "side.neumann.quadratic"
)

# 3-D Tests
ibamr_add_3d_tests(01
  CASES
    "cell.dirichlet.linear"
    "cell.dirichlet.quadratic"
    "cell.neumann.linear"
    "cell.neumann.quadratic"
    "side.dirichlet.linear"
    "side.dirichlet.quadratic"
    "side.neumann.linear"
    "side.neumann.quadratic"
)
