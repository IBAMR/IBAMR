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

FOREACH(_dir eel2d eel3d falling_sphere flow_past_cylinder
    impulsively_started_cylinder knifefish moving_plate
    oscillating_rigid_cylinder stokes_first_problem
    )
  ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
