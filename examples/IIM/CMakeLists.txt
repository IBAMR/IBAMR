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

FOREACH(_dir ex0 ex1 ex2 ex3 ex4 ex5 ex6 ex7 ex8 ex9 anvil_block ball-in-blocked-tube ball-in-narrowing-tube ball-in-tube concentric_rotating_cylinders eccentric_rotating_cylinders flow-past-star misaligned_sliding_bearing parallel-step-slider-bearing pressurized_circle pressurized_cylinder_3d shearing_parallel_plates)
  ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
