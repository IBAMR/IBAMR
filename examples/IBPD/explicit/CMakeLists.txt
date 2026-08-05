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

FOREACH(_dir Anisotropic_Cooks Cooks_membrane Turek_hron Compressed_block Elastic_band Torsion)
  ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
