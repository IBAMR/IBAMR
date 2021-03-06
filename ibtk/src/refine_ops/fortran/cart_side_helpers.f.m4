c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2020 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

dnl This file only exists to be included into
dnl cart_side_refine2d.f and cart_side_refine3d.f: i.e., it is
dnl never explicitly compiled on its own so that the delta function kernels
dnl can be inlined. This is necessary since Fortran compilers will (unless
dnl we use something fancy like LTO) only do inlining of functions and
dnl subroutines defined in the same translation unit.

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The minmod function of three arguments.
c
c     The minmod function is a function of two or more arguments that
c     takes the value of the argument with the smallest modulus if all
c     arguments have the same sign.  Otherwise it takes the value zero.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function minmod3(a,b,c)
c
      implicit none
c
c     Input.
c
      REAL a,b,c
c
c     minmod(a,b,c)
c
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0).and.(c.ge.0.0d0) ) then
         minmod3 = dmin1(a,b,c)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0).and.(c.le.0.0d0) ) then
         minmod3 = dmax1(a,b,c)
      else
         minmod3 = 0.d0
      endif
c
      return
      end
