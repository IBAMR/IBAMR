c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2022 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the arithmetic average based on the four inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL function a_avg4(a0,a1,a2,a3)
      implicit none
      REAL a0,a1,a2,a3
      a_avg4 = 0.25d0*(a0+a1+a2+a3)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the arithmetic average based on the twelve inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL function a_avg12(a0,a1,a2,a3,
     &                      a4,a5,a6,a7,
     &                      a8,a9,a10,a11)
      implicit none
      REAL a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
      a_avg12 = (1.d0/12.d0)*(a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11)
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the harmonic average based on the two inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL function h_avg2(a0,a1)
      implicit none
      REAL a0,a1
      REAL dmr,nmr, eps
      nmr = 2.d0
      eps = sqrt(epsilon(0.d0))
c     If any of the arguments is zero, then harmonic average is zero
      if (abs(a0) .le. eps .or. abs(a1) .le. eps) then
        h_avg2 = 0.d0
      else
        dmr = 1.d0/a0+1.d0/a1
        h_avg2 = nmr/dmr
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the harmonic average based on the four inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL function h_avg4(a0,a1,a2,a3)
      implicit none
      REAL a0,a1,a2,a3
      REAL dmr,nmr,eps
      nmr = 4.d0
      eps = sqrt(epsilon(0.d0))

c     If any of the arguments is zero, then harmonic average is zero
      if (abs(a0) .le. eps .or. abs(a1) .le. eps .or.
     &    abs(a2) .le. eps .or. abs(a3) .le. eps) then
        h_avg4 = 0.0
      else
        dmr = 1.d0/a0+1.d0/a1+1.d0/a2+1.d0/a3
        h_avg4 = nmr/dmr
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the harmonic average based on the twelve inputs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL function h_avg12(a0,a1,a2,a3,
     &                      a4,a5,a6,a7,
     &                      a8,a9,a10,a11)
      implicit none
      REAL a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
      REAL dmr,nmr,eps
      nmr = 12.d0
      eps = sqrt(epsilon(0.d0))

c     If any of the arguments is zero, then harmonic average is zero
      if (abs(a0) .le. eps .or. abs(a1) .le. eps .or.
     &    abs(a2) .le. eps .or. abs(a3) .le. eps .or.
     &    abs(a4) .le. eps .or. abs(a5) .le. eps .or.
     &    abs(a6) .le. eps .or. abs(a7) .le. eps .or.
     &    abs(a8) .le. eps .or. abs(a9) .le. eps .or.
     &    abs(a10) .le. eps .or. abs(a11) .le. eps) then
          h_avg12 = 0.d0
      else
          dmr = 1.d0/a0+1.d0/a1+1.d0/a2+1.d0/a3+1.d0/a4+1.d0/a5
     &             +1.d0/a6+1.d0/a7+1.d0/a8+1.d0/a9+1.d0/a10+1.d0/a11
          h_avg12 = nmr/dmr
      endif
      return
      end
