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

define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = rot U = [du2/dx1 -du1/dx2, du0/dx2 - du2/dx0, du1/dx0 - du0/dx1].
c
c     Uses centered differences to compute the side centered rot of a
c     edge centered vector field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine etosrot3d(
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

      REAL U0(EDGE3d0(ilower,iupper,U_gcw))
      REAL U1(EDGE3d1(ilower,iupper,U_gcw))
      REAL U2(EDGE3d2(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(SIDE3d0(ilower,iupper,W_gcw))
      REAL W1(SIDE3d1(ilower,iupper,W_gcw))
      REAL W2(SIDE3d2(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the side centered rot of edge centered U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               W0(i0,i1,i2) = (U2(i0,i1+1,i2)-U2(i0,i1,i2))/dx(1) -
     &                        (U1(i0,i1,i2+1)-U1(i0,i1,i2))/dx(2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               W1(i0,i1,i2) = (U0(i0,i1,i2+1)-U0(i0,i1,i2))/dx(2) -
     &                        (U2(i0+1,i1,i2)-U2(i0,i1,i2))/dx(0)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               W2(i0,i1,i2) = (U1(i0+1,i1,i2)-U1(i0,i1,i2))/dx(0) -
     &                        (U0(i0,i1+1,i2)-U0(i0,i1,i2))/dx(1)
            enddo
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

