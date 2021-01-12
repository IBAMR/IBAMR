c ---------------------------------------------------------------------
c
c Copyright (c) 2014 - 2019 by the IBAMR developers
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
c     Computes W = rot U = [dU/dx1, -dU/dx0].
c
c     Uses centered differences to compute the side centered rot of a
c     node centered scalar field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ntosrot2d(
     &     W0,W1,W_gcw,
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

      REAL U(NODE2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(SIDE2d0(ilower,iupper,W_gcw))
      REAL W1(SIDE2d1(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the side centered rot of node centered U.
c

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            W0(i0,i1) = (U(i0,i1+1)-U(i0,i1))/dx(1)
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            W1(i0,i1) = -(U(i0+1,i1)-U(i0,i1))/dx(0)
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = rot U = [dU/dy, -dU/dx]
c
c     Uses centered differences to compute the side centered rot of a
c     cell centered scalar field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosrot2d(
     &     W0,W1,W_gcw,
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

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(SIDE2d0(ilower,iupper,W_gcw))
      REAL W1(SIDE2d1(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
c
c     Compute the side centered rot of cell centered U.
c

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            W0(i0,i1) = (U(i0-1,i1+1)-U(i0-1,i1-1) + 
     &                   U(i0,i1+1)-U(i0,i1-1))/(0.4d1*dx(1))
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            W1(i0,i1) = -(U(i0+1,i1-1)-U(i0-1,i1-1) + 
     &                    U(i0+1,i1)-U(i0-1,i1))/(0.4d1*dx(0))
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
