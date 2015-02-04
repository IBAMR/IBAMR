c
c     Routines to compute discrete diffusive fluxes on patches.
c
c     Created on 12 Jun 2003 by Boyce Griffith
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
c     Computes g = alpha grad U.
c
c     Uses centered differences to compute the face centered diffusive
c     flux of a cell centered variable U.
c
c     GRID ALIGNED version.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofflux3d(
     &     g0,g1,g2,g_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER g_gcw,alpha_gcw,U_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw))

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(FACE3d0(ilower,iupper,g_gcw))
      REAL g1(FACE3d1(ilower,iupper,g_gcw))
      REAL g2(FACE3d2(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER d,i0,i1,i2
      REAL    dU_dx(0:NDIM-1),nfac(0:NDIM-1)
c
c     Compute the face centered diffusive flux of U.
c
      do d = 0,NDIM-1
         nfac(d) = 1.d0/dx(d)
      enddo

      d = 0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               dU_dx(d) = nfac(d)*(U(i0,i1,i2)-U(i0-1,i1,i2))
               g0(i0,i1,i2) = alpha0(i0,i1,i2)*dU_dx(d)
            enddo
         enddo
      enddo

      d = 1
      do i0 = ilower0,iupper0
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1+1
               dU_dx(d) = nfac(d)*(U(i0,i1,i2)-U(i0,i1-1,i2))
               g1(i1,i2,i0) = alpha1(i1,i2,i0)*dU_dx(d)
            enddo
         enddo
      enddo

      d = 2
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            do i2 = ilower2,iupper2+1
               dU_dx(d) = nfac(d)*(U(i0,i1,i2)-U(i0,i1,i2-1))
               g2(i2,i0,i1) = alpha2(i2,i0,i1)*dU_dx(d)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes g = alpha grad U.
c
c     Uses centered differences to compute the face centered diffusive
c     flux of a cell centered variable U.
c
c     NON-GRID ALIGNED version.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofanisoflux3d(
     &     g0,g1,g2,g_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER g_gcw,alpha_gcw,U_gcw

      REAL alpha0(FACE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(FACE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(FACE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(FACE3d0(ilower,iupper,g_gcw))
      REAL g1(FACE3d1(ilower,iupper,g_gcw))
      REAL g2(FACE3d2(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER d,i0,i1,i2
      REAL    dU_dx(0:NDIM-1),nfac(0:NDIM-1),tfac(0:NDIM-1)
c
c     Compute the face centered diffusive flux of U.
c
      do d = 0,NDIM-1
         nfac(d) = 1.d0/dx(d)
         tfac(d) = 1.d0/(4.d0*dx(d))
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               dU_dx(0) = nfac(0)*(U(i0,i1,i2)-U(i0-1,i1,i2))
               dU_dx(1) = tfac(1)*(
     &              U(i0  ,i1+1,i2)-U(i0  ,i1-1,i2)+
     &              U(i0-1,i1+1,i2)-U(i0-1,i1-1,i2))
               dU_dx(2) = tfac(2)*(
     &              U(i0  ,i1,i2+1)-U(i0  ,i1,i2-1)+
     &              U(i0-1,i1,i2+1)-U(i0-1,i1,i2-1))

               g0(i0,i1,i2) = 0.d0
               do d = 0,NDIM - 1
                  g0(i0,i1,i2) = alpha0(i0,i1,i2,d)*dU_dx(d)
     &                 + g0(i0,i1,i2)
               enddo
            enddo
         enddo
      enddo

      do i0 = ilower0,iupper0
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1+1
               dU_dx(0) = tfac(0)*(
     &              U(i0+1,i1  ,i2)-U(i0-1,i1  ,i2)+
     &              U(i0+1,i1-1,i2)-U(i0-1,i1-1,i2))
               dU_dx(1) = nfac(1)*(U(i0,i1,i2)-U(i0,i1-1,i2))
               dU_dx(2) = tfac(2)*(
     &              U(i0,i1  ,i2+1)-U(i0,i1  ,i2-1)+
     &              U(i0,i1-1,i2+1)-U(i0,i1-1,i2-1))

               g1(i1,i2,i0) = 0.d0
               do d = 0,NDIM - 1
                  g1(i1,i2,i0) = alpha1(i1,i2,i0,d)*dU_dx(d)
     &                 + g1(i1,i2,i0)
               enddo
            enddo
         enddo
      enddo

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            do i2 = ilower2,iupper2+1
               dU_dx(0) = tfac(0)*(
     &              U(i0+1,i1,i2  )-U(i0-1,i1,i2  )+
     &              U(i0+1,i1,i2-1)-U(i0-1,i1,i2-1))
               dU_dx(1) = tfac(1)*(
     &              U(i0,i1+1,i2  )-U(i0,i1-1,i2  )+
     &              U(i0,i1+1,i2-1)-U(i0,i1-1,i2-1))
               dU_dx(2) = nfac(2)*(U(i0,i1,i2)-U(i0,i1,i2-1))

               g2(i2,i0,i1) = 0.d0
               do d = 0,NDIM - 1
                  g2(i2,i0,i1) = alpha2(i2,i0,i1,d)*dU_dx(d)
     &                 + g2(i2,i0,i1)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes g = alpha grad U.
c
c     Uses centered differences to compute the side centered diffusive
c     flux of a cell centered variable U.
c
c     GRID ALIGNED version.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosflux3d(
     &     g0,g1,g2,g_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER g_gcw,alpha_gcw,U_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw))
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw))
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw))

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(SIDE3d0(ilower,iupper,g_gcw))
      REAL g1(SIDE3d1(ilower,iupper,g_gcw))
      REAL g2(SIDE3d2(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER d,i0,i1,i2
      REAL    dU_dx(0:NDIM-1),nfac(0:NDIM-1)
c
c     Compute the side centered diffusive flux of U.
c
      do d = 0,NDIM-1
         nfac(d) = 1.d0/dx(d)
      enddo

      d = 0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               dU_dx(d) = nfac(d)*(U(i0,i1,i2)-U(i0-1,i1,i2))
               g0(i0,i1,i2) = alpha0(i0,i1,i2)*dU_dx(d)
            enddo
         enddo
      enddo

      d = 1
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               dU_dx(d) = nfac(d)*(U(i0,i1,i2)-U(i0,i1-1,i2))
               g1(i0,i1,i2) = alpha1(i0,i1,i2)*dU_dx(d)
            enddo
         enddo
      enddo

      d = 2
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               dU_dx(d) = nfac(d)*(U(i0,i1,i2)-U(i0,i1,i2-1))
               g2(i0,i1,i2) = alpha2(i0,i1,i2)*dU_dx(d)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes g = alpha grad U.
c
c     Uses centered differences to compute the side centered diffusive
c     flux of a cell centered variable U.
c
c     NON-GRID ALIGNED version.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosanisoflux3d(
     &     g0,g1,g2,g_gcw,
     &     alpha0,alpha1,alpha2,alpha_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER g_gcw,alpha_gcw,U_gcw

      REAL alpha0(SIDE3d0(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha1(SIDE3d1(ilower,iupper,alpha_gcw),0:NDIM-1)
      REAL alpha2(SIDE3d2(ilower,iupper,alpha_gcw),0:NDIM-1)

      REAL U(CELL3d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(SIDE3d0(ilower,iupper,g_gcw))
      REAL g1(SIDE3d1(ilower,iupper,g_gcw))
      REAL g2(SIDE3d2(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER d,i0,i1,i2
      REAL    dU_dx(0:NDIM-1),nfac(0:NDIM-1),tfac(0:NDIM-1)
c
c     Compute the side centered diffusive flux of U.
c
      do d = 0,NDIM-1
         nfac(d) = 1.d0/dx(d)
         tfac(d) = 1.d0/(4.d0*dx(d))
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               dU_dx(0) = nfac(0)*(U(i0,i1,i2)-U(i0-1,i1,i2))
               dU_dx(1) = tfac(1)*(
     &              U(i0  ,i1+1,i2)-U(i0  ,i1-1,i2)+
     &              U(i0-1,i1+1,i2)-U(i0-1,i1-1,i2))
               dU_dx(2) = tfac(2)*(
     &              U(i0  ,i1,i2+1)-U(i0  ,i1,i2-1)+
     &              U(i0-1,i1,i2+1)-U(i0-1,i1,i2-1))

               g0(i0,i1,i2) = 0.d0
               do d = 0,NDIM - 1
                  g0(i0,i1,i2) = alpha0(i0,i1,i2,d)*dU_dx(d)
     &                 + g0(i0,i1,i2)
               enddo
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               dU_dx(0) = tfac(0)*(
     &              U(i0+1,i1  ,i2)-U(i0-1,i1  ,i2)+
     &              U(i0+1,i1-1,i2)-U(i0-1,i1-1,i2))
               dU_dx(1) = nfac(1)*(U(i0,i1,i2)-U(i0,i1-1,i2))
               dU_dx(2) = tfac(2)*(
     &              U(i0,i1  ,i2+1)-U(i0,i1  ,i2-1)+
     &              U(i0,i1-1,i2+1)-U(i0,i1-1,i2-1))

               g1(i0,i1,i2) = 0.d0
               do d = 0,NDIM - 1
                  g1(i0,i1,i2) = alpha1(i0,i1,i2,d)*dU_dx(d)
     &                 + g1(i0,i1,i2)
               enddo
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               dU_dx(0) = tfac(0)*(
     &              U(i0+1,i1,i2  )-U(i0-1,i1,i2  )+
     &              U(i0+1,i1,i2-1)-U(i0-1,i1,i2-1))
               dU_dx(1) = tfac(1)*(
     &              U(i0,i1+1,i2  )-U(i0,i1-1,i2  )+
     &              U(i0,i1+1,i2-1)-U(i0,i1-1,i2-1))
               dU_dx(2) = nfac(2)*(U(i0,i1,i2)-U(i0,i1,i2-1))

               g2(i0,i1,i2) = 0.d0
               do d = 0,NDIM - 1
                  g2(i0,i1,i2) = alpha2(i0,i1,i2,d)*dU_dx(d)
     &                 + g2(i0,i1,i2)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
