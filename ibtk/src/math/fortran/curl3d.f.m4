c
c     Routines to compute discrete curls on patches.
c
c     Created on 01 Jul 2004 by Boyce Griffith
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
c     Computes W = curl U.
c
c     Uses centered differences to compute the cell centered curl of a
c     cell centered vector field U=(U0,U1,U2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoccurl3d(
     &     W,W_gcw,
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
      INTEGER W_gcw,U_gcw

      REAL U(CELL3d(ilower,iupper,U_gcw),0:NDIM-1)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(CELL3d(ilower,iupper,W_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    dU0_dx1,dU0_dx2,fac01,fac02
      REAL    dU1_dx0,dU1_dx2,fac10,fac12
      REAL    dU2_dx0,dU2_dx1,fac20,fac21
c
c     Compute the cell centered curl of U=(U0,U1,U2).
c
      fac12 = 0.5d0/dx(2)
      fac21 = 0.5d0/dx(1)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               dU1_dx2 = fac12*(
     &              U(i0  ,i1  ,i2+1,1)-U(i0  ,i1  ,i2-1,1) )
               dU2_dx1 = fac21*(
     &              U(i0  ,i1+1,i2  ,2)-U(i0  ,i1-1,i2  ,2) )
               W(i0,i1,i2,0) = dU2_dx1-dU1_dx2
            enddo
         enddo
      enddo

      fac02 = 0.5d0/dx(2)
      fac20 = 0.5d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               dU0_dx2 = fac02*(
     &              U(i0  ,i1  ,i2+1,0)-U(i0  ,i1  ,i2-1,0) )
               dU2_dx0 = fac20*(
     &              U(i0+1,i1  ,i2  ,2)-U(i0-1,i1  ,i2  ,2) )
               W(i0,i1,i2,1) = dU0_dx2-dU2_dx0
            enddo
         enddo
      enddo

      fac01 = 0.5d0/dx(1)
      fac10 = 0.5d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               dU0_dx1 = fac01*(
     &              U(i0  ,i1+1,i2  ,0)-U(i0  ,i1-1,i2  ,0) )
               dU1_dx0 = fac10*(
     &              U(i0+1,i1  ,i2  ,1)-U(i0-1,i1  ,i2  ,1) )
               W(i0,i1,i2,2) = dU1_dx0-dU0_dx1
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u.
c
c     Uses centered differences to compute the cell centered curl of a
c     face centered vector field u=(u0,u1,u2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftoccurl3d(
     &     W,W_gcw,
     &     u0,u1,u2,u_gcw,
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
      INTEGER W_gcw,u_gcw

      REAL u0(FACE3d0(ilower,iupper,u_gcw))
      REAL u1(FACE3d1(ilower,iupper,u_gcw))
      REAL u2(FACE3d2(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(CELL3d(ilower,iupper,W_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    du0_dx1,du0_dx2,fac01,fac02
      REAL    du1_dx0,du1_dx2,fac10,fac12
      REAL    du2_dx0,du2_dx1,fac20,fac21
c
c     Compute the cell centered curl of u=(u0,u1,u2).
c
      fac12 = 0.25d0/dx(2)
      fac21 = 0.25d0/dx(1)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du1_dx2 = fac12*(
     &              +u1(i1  ,i2+1,i0  )+u1(i1+1,i2+1,i0  )
     &              -u1(i1  ,i2-1,i0  )-u1(i1+1,i2-1,i0  ) )
               du2_dx1 = fac21*(
     &              +u2(i2  ,i0  ,i1+1)+u2(i2+1,i0  ,i1+1)
     &              -u2(i2  ,i0  ,i1-1)-u2(i2+1,i0  ,i1-1) )
               W(i0,i1,i2,0) = du2_dx1-du1_dx2
            enddo
         enddo
      enddo

      fac02 = 0.25d0/dx(2)
      fac20 = 0.25d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du0_dx2 = fac02*(
     &              +u0(i0  ,i1  ,i2+1)+u0(i0+1,i1  ,i2+1)
     &              -u0(i0  ,i1  ,i2-1)-u0(i0+1,i1  ,i2-1) )
               du2_dx0 = fac20*(
     &              +u2(i2  ,i0+1,i1  )+u2(i2+1,i0+1,i1  )
     &              -u2(i2  ,i0-1,i1  )-u2(i2+1,i0-1,i1  ) )
               W(i0,i1,i2,1) = du0_dx2-du2_dx0
            enddo
         enddo
      enddo

      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du0_dx1 = fac01*(
     &              +u0(i0  ,i1+1,i2  )+u0(i0+1,i1+1,i2  )
     &              -u0(i0  ,i1-1,i2  )-u0(i0+1,i1-1,i2  ) )
               du1_dx0 = fac10*(
     &              +u1(i1  ,i2  ,i0+1)+u1(i1+1,i2  ,i0+1)
     &              -u1(i1  ,i2  ,i0-1)-u1(i1+1,i2  ,i0-1) )
               W(i0,i1,i2,2) = du1_dx0-du0_dx1
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes w = curl u.
c
c     Uses centered differences to compute the face centered curl of a
c     face centered vector field u=(u0,u1,u2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftofcurl3d(
     &     w0,w1,w2,w_gcw,
     &     u0,u1,u2,u_gcw,
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
      INTEGER w_gcw,u_gcw

      REAL u0(FACE3d0(ilower,iupper,u_gcw))
      REAL u1(FACE3d1(ilower,iupper,u_gcw))
      REAL u2(FACE3d2(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL w0(FACE3d0(ilower,iupper,w_gcw))
      REAL w1(FACE3d1(ilower,iupper,w_gcw))
      REAL w2(FACE3d2(ilower,iupper,w_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    du0_dx1,du0_dx2,fac01,fac02
      REAL    du1_dx0,du1_dx2,fac10,fac12
      REAL    du2_dx0,du2_dx1,fac20,fac21
c
c     Compute the face centered curl of u=(u0,u1,u2).
c
      fac12 = 0.125d0/dx(2)
      fac21 = 0.125d0/dx(1)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               du2_dx1 = fac21*(
     &              u2(i2  ,i0-1,i1+1) - u2(i2  ,i0-1,i1-1) +
     &              u2(i2  ,i0  ,i1+1) - u2(i2  ,i0  ,i1-1) +
     &              u2(i2+1,i0-1,i1+1) - u2(i2+1,i0-1,i1-1) +
     &              u2(i2+1,i0  ,i1+1) - u2(i2+1,i0  ,i1-1) )
               du1_dx2 = fac12*(
     &              u1(i1  ,i2+1,i0-1) - u1(i1  ,i2-1,i0-1) +
     &              u1(i1  ,i2+1,i0  ) - u1(i1  ,i2-1,i0  ) +
     &              u1(i1+1,i2+1,i0-1) - u1(i1+1,i2-1,i0-1) +
     &              u1(i1+1,i2+1,i0  ) - u1(i1+1,i2-1,i0  ) )
               w0(i0,i1,i2) = du2_dx1-du1_dx2
            enddo
         enddo
      enddo

      fac02 = 0.125d0/dx(2)
      fac20 = 0.125d0/dx(0)

      do i0 = ilower0,iupper0
         do i2 = ilower2,iupper2
            do i1 = ilower1,iupper1+1
               du0_dx2 = fac02*(
     &              u0(i0  ,i1-1,i2+1) - u0(i0  ,i1-1,i2-1) +
     &              u0(i0+1,i1-1,i2+1) - u0(i0+1,i1-1,i2-1) +
     &              u0(i0  ,i1  ,i2+1) - u0(i0  ,i1  ,i2-1) +
     &              u0(i0+1,i1  ,i2+1) - u0(i0+1,i1  ,i2-1) )
               du2_dx0 = fac20*(
     &              u2(i2  ,i0+1,i1-1) - u2(i2  ,i0-1,i1-1) +
     &              u2(i2  ,i0+1,i1  ) - u2(i2  ,i0-1,i1  ) +
     &              u2(i2+1,i0+1,i1-1) - u2(i2+1,i0-1,i1-1) +
     &              u2(i2+1,i0+1,i1  ) - u2(i2+1,i0-1,i1  ) )
               w1(i1,i2,i0) = du0_dx2-du2_dx0
            enddo
         enddo
      enddo

      fac01 = 0.125d0/dx(1)
      fac10 = 0.125d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            do i2 = ilower2,iupper2+1
               du1_dx0 = fac10*(
     &              u1(i1  ,i2-1,i0+1) - u1(i1  ,i2-1,i0-1) +
     &              u1(i1+1,i2-1,i0+1) - u1(i1+1,i2-1,i0-1) +
     &              u1(i1  ,i2  ,i0+1) - u1(i1  ,i2  ,i0-1) +
     &              u1(i1+1,i2  ,i0+1) - u1(i1+1,i2  ,i0-1) )
               du0_dx1 = fac01*(
     &              u0(i0  ,i1+1,i2-1) - u0(i0  ,i1-1,i2-1) +
     &              u0(i0+1,i1+1,i2-1) - u0(i0+1,i1-1,i2-1) +
     &              u0(i0  ,i1+1,i2  ) - u0(i0  ,i1-1,i2  ) +
     &              u0(i0+1,i1+1,i2  ) - u0(i0+1,i1-1,i2  ) )
               w2(i2,i0,i1) = du1_dx0-du0_dx1
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u.
c
c     Uses centered differences to compute the cell centered curl of a
c     side centered vector field u=(u0,u1,u2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoccurl3d(
     &     W,W_gcw,
     &     u0,u1,u2,u_gcw,
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
      INTEGER W_gcw,u_gcw

      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(CELL3d(ilower,iupper,W_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    du0_dx1,du0_dx2,fac01,fac02
      REAL    du1_dx0,du1_dx2,fac10,fac12
      REAL    du2_dx0,du2_dx1,fac20,fac21
c
c     Compute the cell centered curl of u=(u0,u1,u2).
c
      fac12 = 0.25d0/dx(2)
      fac21 = 0.25d0/dx(1)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du1_dx2 = fac12*(
     &              +u1(i0  ,i1  ,i2+1)+u1(i0  ,i1+1,i2+1)
     &              -u1(i0  ,i1  ,i2-1)-u1(i0  ,i1+1,i2-1) )
               du2_dx1 = fac21*(
     &              +u2(i0  ,i1+1,i2  )+u2(i0  ,i1+1,i2+1)
     &              -u2(i0  ,i1-1,i2  )-u2(i0  ,i1-1,i2+1) )
               W(i0,i1,i2,0) = du2_dx1-du1_dx2
            enddo
         enddo
      enddo

      fac02 = 0.25d0/dx(2)
      fac20 = 0.25d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du0_dx2 = fac02*(
     &              +u0(i0  ,i1  ,i2+1)+u0(i0+1,i1  ,i2+1)
     &              -u0(i0  ,i1  ,i2-1)-u0(i0+1,i1  ,i2-1) )
               du2_dx0 = fac20*(
     &              +u2(i0+1,i1  ,i2  )+u2(i0+1,i1  ,i2+1)
     &              -u2(i0-1,i1  ,i2  )-u2(i0-1,i1  ,i2+1) )
               W(i0,i1,i2,1) = du0_dx2-du2_dx0
            enddo
         enddo
      enddo

      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du0_dx1 = fac01*(
     &              +u0(i0  ,i1+1,i2  )+u0(i0+1,i1+1,i2  )
     &              -u0(i0  ,i1-1,i2  )-u0(i0+1,i1-1,i2  ) )
               du1_dx0 = fac10*(
     &              +u1(i0+1,i1  ,i2  )+u1(i0+1,i1+1,i2  )
     &              -u1(i0-1,i1  ,i2  )-u1(i0-1,i1+1,i2  ) )
               W(i0,i1,i2,2) = du1_dx0-du0_dx1
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes w = curl u.
c
c     Uses centered differences to compute the side centered curl of a
c     side centered vector field u=(u0,u1,u2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoscurl3d(
     &     w0,w1,w2,w_gcw,
     &     u0,u1,u2,u_gcw,
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
      INTEGER w_gcw,u_gcw

      REAL u0(SIDE3d0(ilower,iupper,u_gcw))
      REAL u1(SIDE3d1(ilower,iupper,u_gcw))
      REAL u2(SIDE3d2(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL w0(SIDE3d0(ilower,iupper,w_gcw))
      REAL w1(SIDE3d1(ilower,iupper,w_gcw))
      REAL w2(SIDE3d2(ilower,iupper,w_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    du0_dx1,du0_dx2,fac01,fac02
      REAL    du1_dx0,du1_dx2,fac10,fac12
      REAL    du2_dx0,du2_dx1,fac20,fac21
c
c     Compute the side centered curl of u=(u0,u1,u2).
c
      fac12 = 0.125d0/dx(2)
      fac21 = 0.125d0/dx(1)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               du2_dx1 = fac21*(
     &              u2(i0-1,i1+1,i2  ) - u2(i0-1,i1-1,i2  ) +
     &              u2(i0  ,i1+1,i2  ) - u2(i0  ,i1-1,i2  ) +
     &              u2(i0-1,i1+1,i2+1) - u2(i0-1,i1-1,i2+1) +
     &              u2(i0  ,i1+1,i2+1) - u2(i0  ,i1-1,i2+1) )
               du1_dx2 = fac12*(
     &              u1(i0-1,i1  ,i2+1) - u1(i0-1,i1  ,i2-1) +
     &              u1(i0  ,i1  ,i2+1) - u1(i0  ,i1  ,i2-1) +
     &              u1(i0-1,i1+1,i2+1) - u1(i0-1,i1+1,i2-1) +
     &              u1(i0  ,i1+1,i2+1) - u1(i0  ,i1+1,i2-1) )
               w0(i0,i1,i2) = du2_dx1-du1_dx2
            enddo
         enddo
      enddo

      fac02 = 0.125d0/dx(2)
      fac20 = 0.125d0/dx(0)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               du0_dx2 = fac02*(
     &              u0(i0  ,i1-1,i2+1) - u0(i0  ,i1-1,i2-1) +
     &              u0(i0+1,i1-1,i2+1) - u0(i0+1,i1-1,i2-1) +
     &              u0(i0  ,i1  ,i2+1) - u0(i0  ,i1  ,i2-1) +
     &              u0(i0+1,i1  ,i2+1) - u0(i0+1,i1  ,i2-1) )
               du2_dx0 = fac20*(
     &              u2(i0+1,i1-1,i2  ) - u2(i0-1,i1-1,i2  ) +
     &              u2(i0+1,i1  ,i2  ) - u2(i0-1,i1  ,i2  ) +
     &              u2(i0+1,i1-1,i2+1) - u2(i0-1,i1-1,i2+1) +
     &              u2(i0+1,i1  ,i2+1) - u2(i0-1,i1  ,i2+1) )
               w1(i0,i1,i2) = du0_dx2-du2_dx0
            enddo
         enddo
      enddo

      fac01 = 0.125d0/dx(1)
      fac10 = 0.125d0/dx(0)

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               du1_dx0 = fac10*(
     &              u1(i0+1,i1  ,i2-1) - u1(i0-1,i1  ,i2-1) +
     &              u1(i0+1,i1+1,i2-1) - u1(i0-1,i1+1,i2-1) +
     &              u1(i0+1,i1  ,i2  ) - u1(i0-1,i1  ,i2  ) +
     &              u1(i0+1,i1+1,i2  ) - u1(i0-1,i1+1,i2  ) )
               du0_dx1 = fac01*(
     &              u0(i0  ,i1+1,i2-1) - u0(i0  ,i1-1,i2-1) +
     &              u0(i0+1,i1+1,i2-1) - u0(i0+1,i1-1,i2-1) +
     &              u0(i0  ,i1+1,i2  ) - u0(i0  ,i1-1,i2  ) +
     &              u0(i0+1,i1+1,i2  ) - u0(i0+1,i1-1,i2  ) )
               w2(i0,i1,i2) = du1_dx0-du0_dx1
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes w = curl u.
c
c     Uses centered differences to compute the edge centered curl of a
c     side centered vector field u=(u0,u1,u2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoecurl3d(
     &     W0,W1,W2,W_gcw,
     &     U0,U1,U2,U_gcw,
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
      INTEGER W_gcw,U_gcw

      REAL U0(SIDE3d0(ilower,iupper,U_gcw))
      REAL U1(SIDE3d1(ilower,iupper,U_gcw))
      REAL U2(SIDE3d2(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(EDGE3d0(ilower,iupper,W_gcw))
      REAL W1(EDGE3d1(ilower,iupper,W_gcw))
      REAL W2(EDGE3d2(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the edge centered curl of u=(u0,u1,u2).
c

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               W0(i0,i1,i2) = (U2(i0,i1,i2) - U2(i0,i1-1,i2)  )/dx(1)  
     &                      - (U1(i0,i1,i2) - U1(i0,i1  ,i2-1))/dx(2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               W1(i0,i1,i2) = (U0(i0,i1,i2) - U0(i0,i1,i2-1))/dx(2)  
     &                      - (U2(i0,i1,i2) - U2(i0-1,i1,i2))/dx(0)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               W2(i0,i1,i2) = (U1(i0,i1,i2) - U1(i0-1,i1,i2))/dx(0) 
     &                      - (U0(i0,i1,i2) - U0(i0,i1-1,i2))/dx(1)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
