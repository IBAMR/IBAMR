c
c     Routines to refine values via divergence- and curl-preserving
c     interpolation.
c
c     Created on 02 Jul 2009 by Boyce Griffith
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply the divergence- and gradient-preserving correction to values
c     refined from the next coarser level of the patch hierarchy.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine div_preserving_correction3d(
     &     u0,u1,u2,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     correction_box_ilower0,correction_box_iupper0,
     &     correction_box_ilower1,correction_box_iupper1,
     &     correction_box_ilower2,correction_box_iupper2,
     &     ratio,dx_fine)
c
      implicit none
c
c     Input.
c
      INTEGER u_gcw
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER correction_box_ilower0,correction_box_iupper0
      INTEGER correction_box_ilower1,correction_box_iupper1
      INTEGER correction_box_ilower2,correction_box_iupper2
      INTEGER ratio(0:NDIM-1)
      REAL    dx_fine(0:NDIM-1)
c
c     Input/Output.
c
      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER d,i0,i1,i2,i,j,k
      REAL u(-2:2,-1:1,-1:1),u_xx,u_xyz
      REAL v(-1:1,-2:2,-1:1),v_yy,v_xyz
      REAL w(-1:1,-1:1,-2:2),w_zz,w_xyz
      REAL dx,dy,dz
c
c     Apply the divergence- and curl-preserving corrections.
c
      do d = 0,NDIM-1
         if ( .not.(ratio(d).eq.2) ) then
            print *,'error: invalid refinement ratio'
            call abort
         endif
      enddo

      dx = dx_fine(0)
      dy = dx_fine(1)
      dz = dx_fine(2)

      do i2=correction_box_ilower2,correction_box_iupper2,ratio(2)
         do i1=correction_box_ilower1,correction_box_iupper1,ratio(1)
            do i0=correction_box_ilower0,correction_box_iupper0,ratio(0)

               u(-2,-1,-1) = u0(i0  ,i1  ,i2  )
               u( 2,-1,-1) = u0(i0+2,i1  ,i2  )
               u(-2, 1,-1) = u0(i0  ,i1+1,i2  )
               u( 2, 1,-1) = u0(i0+2,i1+1,i2  )
               u(-2,-1, 1) = u0(i0  ,i1  ,i2+1)
               u( 2,-1, 1) = u0(i0+2,i1  ,i2+1)
               u(-2, 1, 1) = u0(i0  ,i1+1,i2+1)
               u( 2, 1, 1) = u0(i0+2,i1+1,i2+1)

               v(-1,-2,-1) = u1(i0  ,i1  ,i2  )
               v( 1,-2,-1) = u1(i0+1,i1  ,i2  )
               v(-1, 2,-1) = u1(i0  ,i1+2,i2  )
               v( 1, 2,-1) = u1(i0+1,i1+2,i2  )
               v(-1,-2, 1) = u1(i0  ,i1  ,i2+1)
               v( 1,-2, 1) = u1(i0+1,i1  ,i2+1)
               v(-1, 2, 1) = u1(i0  ,i1+2,i2+1)
               v( 1, 2, 1) = u1(i0+1,i1+2,i2+1)

               w(-1,-1,-2) = u2(i0  ,i1  ,i2  )
               w( 1,-1,-2) = u2(i0+1,i1  ,i2  )
               w(-1, 1,-2) = u2(i0  ,i1+1,i2  )
               w( 1, 1,-2) = u2(i0+1,i1+1,i2  )
               w(-1,-1, 2) = u2(i0  ,i1  ,i2+2)
               w( 1,-1, 2) = u2(i0+1,i1  ,i2+2)
               w(-1, 1, 2) = u2(i0  ,i1+1,i2+2)
               w( 1, 1, 2) = u2(i0+1,i1+1,i2+2)

               u_xx = 0.d0;
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        u_xx = u_xx + 0.125d0*(
     &                       dble(i)*dble(j)*v(i,2*j,k) +
     &                       dble(i)*dble(k)*w(i,j,2*k) )
                     enddo
                  enddo
               enddo

               v_yy = 0.d0;
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        v_yy = v_yy + 0.125d0*(
     &                       dble(i)*dble(j)*u(2*i,j,k) +
     &                       dble(j)*dble(k)*w(i,j,2*k) )
                     enddo
                  enddo
               enddo

               w_zz = 0.d0;
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        w_zz = w_zz + 0.125d0*(
     &                       dble(i)*dble(k)*u(2*i,j,k) +
     &                       dble(j)*dble(k)*v(i,2*j,k) )
                     enddo
                  enddo
               enddo

               u_xyz = 0.d0
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        u_xyz = u_xyz + 0.125d0*(
     &                       dble(i)*dble(j)*dble(k)*u(2*i,j,k)/
     &                       (dy**2.d0 + dz**2.d0))
                     enddo
                  enddo
               enddo

               v_xyz = 0.d0
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        v_xyz = v_xyz + 0.125d0*(
     &                       dble(i)*dble(j)*dble(k)*v(i,2*j,k)/
     &                       (dx**2.d0 + dz**2.d0))
                     enddo
                  enddo
               enddo

               w_xyz = 0.d0
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        w_xyz = w_xyz + 0.125d0*(
     &                       dble(i)*dble(j)*dble(k)*w(i,j,2*k)/
     &                       (dx**2.d0 + dy**2.d0))
                     enddo
                  enddo
               enddo

               do k = -1,1,2
                  do j = -1,1,2
                     u(0,j,k) = 0.5d0*(u(2,j,k)+u(-2,j,k)) + u_xx +
     &                    dble(k)*dz**2.d0*v_xyz +
     &                    dble(j)*dy**2.d0*w_xyz
                  enddo
               enddo

               do k = -1,1,2
                  do i = -1,1,2
                     v(i,0,k) = 0.5d0*(v(i,2,k)+v(i,-2,k)) + v_yy +
     &                    dble(i)*dx**2.d0*w_xyz +
     &                    dble(k)*dz**2.d0*u_xyz
                  enddo
               enddo

               do j = -1,1,2
                  do i = -1,1,2
                     w(i,j,0) = 0.5d0*(w(i,j,2)+w(i,j,-2)) + w_zz +
     &                    dble(j)*dy**2.d0*u_xyz +
     &                    dble(i)*dx**2.d0*v_xyz
                  enddo
               enddo

               u0(i0+1,i1  ,i2  ) = u(0,-1,-1)
               u0(i0+1,i1+1,i2  ) = u(0, 1,-1)
               u0(i0+1,i1  ,i2+1) = u(0,-1, 1)
               u0(i0+1,i1+1,i2+1) = u(0, 1, 1)

               u1(i0  ,i1+1,i2  ) = v(-1,0,-1)
               u1(i0+1,i1+1,i2  ) = v( 1,0,-1)
               u1(i0  ,i1+1,i2+1) = v(-1,0, 1)
               u1(i0+1,i1+1,i2+1) = v( 1,0, 1)

               u2(i0  ,i1  ,i2+1) = w(-1,-1,0)
               u2(i0+1,i1  ,i2+1) = w( 1,-1,0)
               u2(i0  ,i1+1,i2+1) = w(-1, 1,0)
               u2(i0+1,i1+1,i2+1) = w( 1, 1,0)

            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
