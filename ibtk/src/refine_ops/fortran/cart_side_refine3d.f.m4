c
c     Routines to refine side-centered values.
c
c     Created on 08 Feb 2011 by Boyce Griffith
c
c     Copyright (c) 2002-2014, Boyce Griffith
c     All rights reserved.
c
c     Redistribution and use in source and binary forms, with or without
c     modification, are permitted provided that the following conditions
c     are met:
c
c        * Redistributions of source code must retain the above
c          copyright notice, this list of conditions and the following
c          disclaimer.
c
c        * Redistributions in binary form must reproduce the above
c          copyright notice, this list of conditions and the following
c          disclaimer in the documentation and/or other materials
c          provided with the distribution.
c
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
c
c     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
c     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
c     BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
c     TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
c     ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
c     TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
define(coarsen_index,`dnl
         if ($1.lt.0) then
            $2=($1+1)/$3-1
         else
            $2=$1/$3
         endif
')dnl'
define(d_dx0_C,0.5d0*($1($2+1,$3,$4)-$1($2-1,$3,$4)))dnl
define(d_dx0_L,($1($2,$3,$4)-$1($2-1,$3,$4)))dnl
define(d_dx0_R,($1($2+1,$3,$4)-$1($2,$3,$4)))dnl
define(d_dx0_MC,`dnl
minmod3(
     & d_dx0_C($1,$2,$3,$4),
     & 2.d0*d_dx0_L($1,$2,$3,$4),
     & 2.d0*d_dx0_R($1,$2,$3,$4))')dnl'
define(d_dx1_C,0.5d0*($1($2,$3+1,$4)-$1($2,$3-1,$4)))dnl
define(d_dx1_L,($1($2,$3,$4)-$1($2,$3-1,$4)))dnl
define(d_dx1_R,($1($2,$3+1,$4)-$1($2,$3,$4)))dnl
define(d_dx1_MC,`dnl
minmod3(
     & d_dx1_C($1,$2,$3,$4),
     & 2.d0*d_dx1_L($1,$2,$3,$4),
     & 2.d0*d_dx1_R($1,$2,$3,$4))')dnl'
define(d_dx2_C,0.5d0*($1($2,$3,$4+1)-$1($2,$3,$4-1)))dnl
define(d_dx2_L,($1($2,$3,$4)-$1($2,$3,$4-1)))dnl
define(d_dx2_R,($1($2,$3,$4+1)-$1($2,$3,$4)))dnl
define(d_dx2_MC,`dnl
minmod3(
     & d_dx2_C($1,$2,$3,$4),
     & 2.d0*d_dx2_L($1,$2,$3,$4),
     & 2.d0*d_dx2_R($1,$2,$3,$4))')dnl'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform specialized refine operation that employs constant
c     prolongation followed by linear interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_specialized_constant_refine3d(
     &     u0_f,u1_f,u2_f,u_f_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     u0_c,u1_c,u2_c,u_c_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     fill_lower0,fill_upper0,
     &     fill_lower1,fill_upper1,
     &     fill_lower2,fill_upper2,
     &     ratio)
c
      implicit none
c
c     Input.
c
      INTEGER u_f_gcw
      INTEGER flower0,fupper0
      INTEGER flower1,fupper1
      INTEGER flower2,fupper2
      INTEGER u_c_gcw
      INTEGER clower0,cupper0
      INTEGER clower1,cupper1
      INTEGER clower2,cupper2
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER fill_lower0,fill_upper0
      INTEGER fill_lower1,fill_upper1
      INTEGER fill_lower2,fill_upper2
      INTEGER ratio(0:NDIM-1)

      REAL u0_c(SIDE3d0(clower,cupper,u_c_gcw))
      REAL u1_c(SIDE3d1(clower,cupper,u_c_gcw))
      REAL u2_c(SIDE3d2(clower,cupper,u_c_gcw))
c
c     Output.
c
      REAL u0_f(SIDE3d0(flower,fupper,u_f_gcw))
      REAL u1_f(SIDE3d1(flower,fupper,u_f_gcw))
      REAL u2_f(SIDE3d2(flower,fupper,u_f_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      INTEGER i_c0,i_c1,i_c2
      REAL w0,w1
c
c     Refine data.
c
      do i2=fill_lower2,fill_upper2
         coarsen_index(i2,i_c2,ratio(2))
         do i1=fill_lower1,fill_upper1
            coarsen_index(i1,i_c1,ratio(1))
            do i0=fill_lower0,fill_upper0+1
               coarsen_index(i0,i_c0,ratio(0))
               if ( i0 .ge. ilower0 .and. i0 .le. iupper0+1 .and.
     &              i1 .ge. ilower1 .and. i1 .le. iupper1   .and.
     &              i2 .ge. ilower2 .and. i2 .le. iupper2   ) then
                  w1 = dble(i0-ratio(0)*i_c0)/dble(ratio(0))
                  w0 = 1.d0-w1
                  u0_f(i0,i1,i2) =
     &                 w0*u0_c(i_c0  ,i_c1,i_c2) +
     &                 w1*u0_c(i_c0+1,i_c1,i_c2)
               endif
            enddo
         enddo
      enddo

      do i2=fill_lower2,fill_upper2
         coarsen_index(i2,i_c2,ratio(2))
         do i1=fill_lower1,fill_upper1+1
            coarsen_index(i1,i_c1,ratio(1))
            do i0=fill_lower0,fill_upper0
               coarsen_index(i0,i_c0,ratio(0))
               if ( i0 .ge. ilower0 .and. i0 .le. iupper0   .and.
     &              i1 .ge. ilower1 .and. i1 .le. iupper1+1 .and.
     &              i2 .ge. ilower2 .and. i2 .le. iupper2   ) then
                  w1 = dble(i1-ratio(1)*i_c1)/dble(ratio(1))
                  w0 = 1.d0-w1
                  u1_f(i0,i1,i2) =
     &                 w0*u1_c(i_c0,i_c1  ,i_c2) +
     &                 w1*u1_c(i_c0,i_c1+1,i_c2)
               endif
            enddo
         enddo
      enddo

      do i2=fill_lower2,fill_upper2+1
         coarsen_index(i2,i_c2,ratio(2))
         do i1=fill_lower1,fill_upper1
            coarsen_index(i1,i_c1,ratio(1))
            do i0=fill_lower0,fill_upper0
               coarsen_index(i0,i_c0,ratio(0))
               if ( i0 .ge. ilower0 .and. i0 .le. iupper0   .and.
     &              i1 .ge. ilower1 .and. i1 .le. iupper1   .and.
     &              i2 .ge. ilower2 .and. i2 .le. iupper2+1 ) then
                  w1 = dble(i2-ratio(2)*i_c2)/dble(ratio(2))
                  w0 = 1.d0-w1
                  u2_f(i0,i1,i2) =
     &                 w0*u2_c(i_c0,i_c1,i_c2  ) +
     &                 w1*u2_c(i_c0,i_c1,i_c2+1)
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform specialized refine operation that employs linear
c     interpolation in the normal direction and constant interpolation
c     in the tangential direction.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_specialized_linear_refine3d(
     &     u0_f,u1_f,u2_f,u_fine_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     u0_c,u1_c,u2_c,u_coarse_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
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
      INTEGER u_fine_gcw
      INTEGER flower0,fupper0
      INTEGER flower1,fupper1
      INTEGER flower2,fupper2
      INTEGER u_coarse_gcw
      INTEGER clower0,cupper0
      INTEGER clower1,cupper1
      INTEGER clower2,cupper2
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER ratio(0:NDIM-1)

      REAL u0_c(SIDE3d0(clower,cupper,u_coarse_gcw))
      REAL u1_c(SIDE3d1(clower,cupper,u_coarse_gcw))
      REAL u2_c(SIDE3d2(clower,cupper,u_coarse_gcw))
c
c     Input/Output.
c
      REAL u0_f(SIDE3d0(flower,fupper,u_fine_gcw))
      REAL u1_f(SIDE3d1(flower,fupper,u_fine_gcw))
      REAL u2_f(SIDE3d2(flower,fupper,u_fine_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      INTEGER i_c0,i_c1,i_c2
      INTEGER i_f0,i_f1,i_f2
      REAL w0,w1,w2
      REAL dx0,dx1,dx2
      REAL du0_dx1,du0_dx2,du1_dx0,du1_dx2,du2_dx0,du2_dx1
c
c     Refine data.
c
      dx0 = dble(ratio(0))      ! effective grid spacings on the coarse level
      dx1 = dble(ratio(1))
      dx2 = dble(ratio(2))

      do i2=ilower2,iupper2
         do i1=ilower1,iupper1
            do i0=ilower0,iupper0+1
               coarsen_index(i0,i_c0,ratio(0))
               coarsen_index(i1,i_c1,ratio(1))
               coarsen_index(i2,i_c2,ratio(2))

               i_f0 = i_c0*ratio(0)
               i_f1 = i_c1*ratio(1)
               i_f2 = i_c2*ratio(2)

               w0 = 1.d0-dble(i0-i_f0)/dx0
               w1 = dble(i1-i_f1)+0.5d0-0.5d0*dx1
               w2 = dble(i2-i_f2)+0.5d0-0.5d0*dx2

               du0_dx1 = d_dx1_MC(u0_c,i_c0,i_c1,i_c2)/dx1
               du0_dx2 = d_dx2_MC(u0_c,i_c0,i_c1,i_c2)/dx2
               u0_f(i0,i1,i2) =
     &              w0*(u0_c(i_c0,i_c1,i_c2)+w1*du0_dx1+w2*du0_dx2)

               w0 = 1.d0-w0

               du0_dx1 = d_dx1_MC(u0_c,i_c0+1,i_c1,i_c2)/dx1
               du0_dx2 = d_dx2_MC(u0_c,i_c0+1,i_c1,i_c2)/dx2
               u0_f(i0,i1,i2) = u0_f(i0,i1,i2)+
     &              w0*(u0_c(i_c0+1,i_c1,i_c2)+w1*du0_dx1+w2*du0_dx2)
            enddo
         enddo
      enddo

      do i2=ilower2,iupper2
         do i1=ilower1,iupper1+1
            do i0=ilower0,iupper0
               coarsen_index(i0,i_c0,ratio(0))
               coarsen_index(i1,i_c1,ratio(1))
               coarsen_index(i2,i_c2,ratio(2))

               i_f0 = i_c0*ratio(0)
               i_f1 = i_c1*ratio(1)
               i_f2 = i_c2*ratio(2)

               w0 = dble(i0-i_f0)+0.5d0-0.5d0*dx0
               w1 = 1.d0-dble(i1-i_f1)/dx1
               w2 = dble(i2-i_f2)+0.5d0-0.5d0*dx2

               du1_dx0 = d_dx0_MC(u1_c,i_c0,i_c1,i_c2)/dx0
               du1_dx2 = d_dx2_MC(u1_c,i_c0,i_c1,i_c2)/dx2
               u1_f(i0,i1,i2) =
     &              w1*(u1_c(i_c0,i_c1,i_c2)+w0*du1_dx0+w2*du1_dx2)

               w1 = 1.d0-w1

               du1_dx0 = d_dx0_MC(u1_c,i_c0,i_c1+1,i_c2)/dx0
               du1_dx2 = d_dx2_MC(u1_c,i_c0,i_c1+1,i_c2)/dx2
               u1_f(i0,i1,i2) = u1_f(i0,i1,i2)+
     &              w1*(u1_c(i_c0,i_c1+1,i_c2)+w0*du1_dx0+w2*du1_dx2)
            enddo
         enddo
      enddo

      do i2=ilower2,iupper2+1
         do i1=ilower1,iupper1
            do i0=ilower0,iupper0
               coarsen_index(i0,i_c0,ratio(0))
               coarsen_index(i1,i_c1,ratio(1))
               coarsen_index(i2,i_c2,ratio(2))

               i_f0 = i_c0*ratio(0)
               i_f1 = i_c1*ratio(1)
               i_f2 = i_c2*ratio(2)

               w0 = dble(i0-i_f0)+0.5d0-0.5d0*dx0
               w1 = dble(i1-i_f1)+0.5d0-0.5d0*dx1
               w2 = 1.d0-dble(i2-i_f2)/dx2

               du2_dx0 = d_dx0_MC(u2_c,i_c0,i_c1,i_c2)/dx0
               du2_dx1 = d_dx1_MC(u2_c,i_c0,i_c1,i_c2)/dx1
               u2_f(i0,i1,i2) =
     &              w2*(u2_c(i_c0,i_c1,i_c2)+w0*du2_dx0+w1*du2_dx1)

               w2 = 1.d0-w2

               du2_dx0 = d_dx0_MC(u2_c,i_c0,i_c1,i_c2+1)/dx0
               du2_dx1 = d_dx1_MC(u2_c,i_c0,i_c1,i_c2+1)/dx1
               u2_f(i0,i1,i2) = u2_f(i0,i1,i2)+
     &              w2*(u2_c(i_c0,i_c1,i_c2+1)+w0*du2_dx0+w1*du2_dx1)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
