c ---------------------------------------------------------------------
c
c Copyright (c) 2019 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

      INTEGER ctu_only,minmod_limited,mc_limited,superbee_limited,
     &        muscl_limited,second_order,fourth_order,ppm,xsppm7
      parameter (ctu_only=1)
      parameter (minmod_limited=2)
      parameter (mc_limited=3)
      parameter (superbee_limited=4)
      parameter (muscl_limited=5)
      parameter (second_order=6)
      parameter (fourth_order=7)
      parameter (ppm=8)
      parameter (xsppm7=9)
