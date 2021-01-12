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

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
define(coarsen_index,`dnl
if ($1.lt.0) then
            $2=($1+1)/$4-1
         else
            $2=$1/$4
         endif
         $3=$2*$4
')dnl'
define(d_dx0_C,0.5d0*($1($2+1,$3)-$1($2-1,$3)))dnl
define(d_dx0_L,($1($2,$3)-$1($2-1,$3)))dnl
define(d_dx0_R,($1($2+1,$3)-$1($2,$3)))dnl
define(d_dx0_MC,`dnl
minmod3(
     & d_dx0_C($1,$2,$3),
     & 2.d0*d_dx0_L($1,$2,$3),
     & 2.d0*d_dx0_R($1,$2,$3))')dnl'
define(d_dx1_C,0.5d0*($1($2,$3+1)-$1($2,$3-1)))dnl
define(d_dx1_L,($1($2,$3)-$1($2,$3-1)))dnl
define(d_dx1_R,($1($2,$3+1)-$1($2,$3)))dnl
define(d_dx1_MC,`dnl
minmod3(
     & d_dx1_C($1,$2,$3),
     & 2.d0*d_dx1_L($1,$2,$3),
     & 2.d0*d_dx1_R($1,$2,$3))')dnl'

c     this is an m4 include, not a Fortran include
include(CURRENT_SRCDIR/cart_side_helpers.f.m4)

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform side-centered refine operation based on RT0 interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_rt0_refine2d(
     &     u0_f,u1_f,u_f_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     u0_c,u1_c,u_c_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ratio)
c
      implicit none
c
c     Input.
c
      INTEGER u_f_gcw
      INTEGER flower0,fupper0
      INTEGER flower1,fupper1
      INTEGER u_c_gcw
      INTEGER clower0,cupper0
      INTEGER clower1,cupper1
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ratio(0:NDIM-1)

      REAL u0_c(SIDE2d0(clower,cupper,u_c_gcw))
      REAL u1_c(SIDE2d1(clower,cupper,u_c_gcw))
c
c     Output.
c
      REAL u0_f(SIDE2d0(flower,fupper,u_f_gcw))
      REAL u1_f(SIDE2d1(flower,fupper,u_f_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER i_c0,i_c1
      INTEGER i_f0,i_f1
      REAL w0,w1
c
c     Refine data.
c
      do i1=ilower1,iupper1
         coarsen_index(i1,i_c1,i_f1,ratio(1))
         do i0=ilower0,iupper0+1
            coarsen_index(i0,i_c0,i_f0,ratio(0))
            if ( i0 .eq. i_f0 ) then
               u0_f(i0,i1) = u0_c(i_c0,i_c1)
            else
               w1 = dble(i0-i_f0)/dble(ratio(0))
               w0 = 1.d0-w1
               u0_f(i0,i1) =
     &              w0*u0_c(i_c0  ,i_c1) +
     &              w1*u0_c(i_c0+1,i_c1)
            endif
         enddo
      enddo

      do i1=ilower1,iupper1+1
         coarsen_index(i1,i_c1,i_f1,ratio(1))
         do i0=ilower0,iupper0
            coarsen_index(i0,i_c0,i_f0,ratio(0))
            if ( i1 .eq. i_f1 ) then
               u1_f(i0,i1) = u1_c(i_c0,i_c1)
            else
               w1 = dble(i1-i_f1)/dble(ratio(1))
               w0 = 1.d0-w1
               u1_f(i0,i1) =
     &              w0*u1_c(i_c0,i_c1  ) +
     &              w1*u1_c(i_c0,i_c1+1)
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform specialized refine operation that employs linear
c     interpolation in the normal direction and MC-limited
c     piecewise-linear interpolation in the tangential direction.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_specialized_linear_refine2d(
     &     u0_f,u1_f,u_f_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     u0_c,u1_c,u_c_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ratio)
c
      implicit none
c
c     External functions.
c
      REAL minmod3
c
c     Input.
c
      INTEGER u_f_gcw
      INTEGER flower0,fupper0
      INTEGER flower1,fupper1
      INTEGER u_c_gcw
      INTEGER clower0,cupper0
      INTEGER clower1,cupper1
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ratio(0:NDIM-1)

      REAL u0_c(SIDE2d0(clower,cupper,u_c_gcw))
      REAL u1_c(SIDE2d1(clower,cupper,u_c_gcw))
c
c     Output.
c
      REAL u0_f(SIDE2d0(flower,fupper,u_f_gcw))
      REAL u1_f(SIDE2d1(flower,fupper,u_f_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER i_c0,i_c1
      INTEGER i_f0,i_f1
      REAL w0,w1
      REAL dx0,dx1
      REAL du0_dx1,du1_dx0
c
c     Refine data.
c
      dx0 = dble(ratio(0))      ! effective grid spacings on the coarse level
      dx1 = dble(ratio(1))

      do i1=ilower1,iupper1
         coarsen_index(i1,i_c1,i_f1,ratio(1))
         do i0=ilower0,iupper0+1
            coarsen_index(i0,i_c0,i_f0,ratio(0))
            i_f0 = i_c0*ratio(0)

            w0 = 1.d0-dble(i0-i_f0)/dx0
            w1 = dble(i1-i_f1)+0.5d0-0.5d0*dx1

            du0_dx1 = d_dx1_MC(u0_c,i_c0,i_c1)/dx1
            u0_f(i0,i1) = w0*(u0_c(i_c0,i_c1)+w1*du0_dx1)

            w0 = 1.d0-w0

            du0_dx1 = d_dx1_MC(u0_c,i_c0+1,i_c1)/dx1
            u0_f(i0,i1) = u0_f(i0,i1)+w0*(u0_c(i_c0+1,i_c1)+w1*du0_dx1)
         enddo
      enddo

      do i1=ilower1,iupper1+1
         coarsen_index(i1,i_c1,i_f1,ratio(1))
         do i0=ilower0,iupper0
            coarsen_index(i0,i_c0,i_f0,ratio(0))

            w0 = dble(i0-i_f0)+0.5d0-0.5d0*dx0
            w1 = 1.d0-dble(i1-i_f1)/dx1

            du1_dx0 = d_dx0_MC(u1_c,i_c0,i_c1)/dx0
            u1_f(i0,i1) = w1*(u1_c(i_c0,i_c1)+w0*du1_dx0)

            w1 = 1.d0-w1

            du1_dx0 = d_dx0_MC(u1_c,i_c0,i_c1+1)/dx0
            u1_f(i0,i1) = u1_f(i0,i1)+w1*(u1_c(i_c0,i_c1+1)+w0*du1_dx0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
