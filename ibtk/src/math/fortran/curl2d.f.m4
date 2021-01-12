c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2019 by the IBAMR developers
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl U = dU1/dx0 - dU0/dx1.
c
c     Uses centered differences to compute the cell centered curl of a
c     cell centered vector field U=(U0,U1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoccurl2d(
     &     W,W_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER W_gcw,U_gcw

      REAL U(CELL2d(ilower,iupper,U_gcw),0:NDIM-1)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(CELL2d(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    dU0_dx1,dU1_dx0,fac01,fac10
c
c     Compute the cell centered curl of U=(U0,U1).
c
      fac01 = 0.5d0/dx(1)
      fac10 = 0.5d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            dU0_dx1 = fac01*(U(i0  ,i1+1,0)-U(i0  ,i1-1,0))
            dU1_dx0 = fac10*(U(i0+1,i1  ,1)-U(i0-1,i1  ,1))
            W(i0,i1) = dU1_dx0-dU0_dx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u = du1/dx0 - du0/dx1.
c
c     Uses centered differences to compute the cell centered curl of a
c     face centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftoccurl2d(
     &     W,W_gcw,
     &     u0,u1,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER W_gcw,u_gcw

      REAL u0(FACE2d0(ilower,iupper,u_gcw))
      REAL u1(FACE2d1(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(CELL2d(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    du0_dx1,du1_dx0,fac01,fac10
c
c     Compute the cell centered curl of u=(u0,u1).
c
      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            du0_dx1 = fac01*(
     &           +u0(i0  ,i1+1)+u0(i0+1,i1+1)
     &           -u0(i0  ,i1-1)-u0(i0+1,i1-1) )
            du1_dx0 = fac10*(
     &           +u1(i1  ,i0+1)+u1(i1+1,i0+1)
     &           -u1(i1  ,i0-1)-u1(i1+1,i0-1) )
            W(i0,i1) = du1_dx0-du0_dx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u = du1/dx0 - du0/dx1.
c
c     Uses centered differences to compute the cell centered curl of a
c     side centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoccurl2d(
     &     W,W_gcw,
     &     u0,u1,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER W_gcw,u_gcw

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(CELL2d(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    du0_dx1,du1_dx0,fac01,fac10
c
c     Compute the cell centered curl of u=(u0,u1).
c
      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            du0_dx1 = fac01*(
     &           +u0(i0  ,i1+1)+u0(i0+1,i1+1)
     &           -u0(i0  ,i1-1)-u0(i0+1,i1-1) )
            du1_dx0 = fac10*(
     &           +u1(i0+1,i1  )+u1(i0+1,i1+1)
     &           -u1(i0-1,i1  )-u1(i0-1,i1+1) )
            W(i0,i1) = du1_dx0-du0_dx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u = du1/dx0 - du0/dx1.
c
c     Uses centered differences to compute the node centered curl of a
c     side centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoncurl2d(
     &     W,W_gcw,
     &     U0,U1,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER W_gcw,U_gcw

      REAL U0(SIDE2d0(ilower,iupper,U_gcw))
      REAL U1(SIDE2d1(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W(NODE2d(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the node centered curl of U.
c

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0+1
            W(i0,i1) = (U1(i0,i1) - U1(i0-1,i1))/dx(0) -
     &                 (U0(i0,i1) - U0(i0,i1-1))/dx(1) 
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
