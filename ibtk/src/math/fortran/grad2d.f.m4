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
c     Computes G = alpha grad U.
c
c     Uses centered differences to compute the cell centered total
c     gradient of a cell centered variable U.
c
c     This is a total gradient in the sense that each component of the
c     gradient is computed for each cell center.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgrad2d(
     &     G,G_gcw,
     &     alpha,
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
      INTEGER G_gcw,U_gcw

      REAL alpha

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL G(CELL2d(ilower,iupper,G_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1
c
c     Compute the cell centered total gradient of U.
c
      fac0 = alpha/(2.d0*dx(0))
      fac1 = alpha/(2.d0*dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            G(i0,i1,0) = fac0*(U(i0+1,i1)-U(i0-1,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            G(i0,i1,1) = fac1*(U(i0,i1+1)-U(i0,i1-1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes G = alpha grad U + beta V.
c
c     Uses centered differences to compute the cell centered total
c     gradient of a cell centered variable U.
c
c     This is a total gradient in the sense that each component of the
c     gradient is computed for each cell center.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocgradadd2d(
     &     G,G_gcw,
     &     alpha,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
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
      INTEGER G_gcw,U_gcw,V_gcw

      REAL alpha

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL beta

      REAL V(CELL2d(ilower,iupper,V_gcw),0:NDIM-1)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL G(CELL2d(ilower,iupper,G_gcw),0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1
c
c     Compute the cell centered total gradient of U.
c
      fac0 = alpha/(2.d0*dx(0))
      fac1 = alpha/(2.d0*dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            G(i0,i1,0) = fac0*(U(i0+1,i1)-U(i0-1,i1)) + beta*V(i0,i1,0)
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            G(i0,i1,1) = fac1*(U(i0,i1+1)-U(i0,i1-1)) + beta*V(i0,i1,1)
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
c     Uses centered differences to compute the face centered partial
c     gradient of a cell centered variable U.
c
c     This is a partial gradient in the sense that only the normal
c     component of the gradient is computed on each face.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofgrad2d(
     &     g0,g1,g_gcw,
     &     alpha,
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
      INTEGER g_gcw,U_gcw

      REAL alpha

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(FACE2d0(ilower,iupper,g_gcw))
      REAL g1(FACE2d1(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1
c
c     Compute the face centered partial gradient of U.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            g0(i0,i1) = fac0*(U(i0,i1)-U(i0-1,i1))
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1+1
            g1(i1,i0) = fac1*(U(i0,i1)-U(i0,i1-1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes g = alpha grad U + beta v.
c
c     Uses centered differences to compute the face centered partial
c     gradient of a cell centered variable U.
c
c     This is a partial gradient in the sense that only the normal
c     component of the gradient is computed on each face.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofgradadd2d(
     &     g0,g1,g_gcw,
     &     alpha,
     &     U,U_gcw,
     &     beta,
     &     v0,v1,v_gcw,
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
      INTEGER g_gcw,U_gcw,v_gcw

      REAL alpha

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL beta

      REAL v0(FACE2d0(ilower,iupper,v_gcw))
      REAL v1(FACE2d1(ilower,iupper,v_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(FACE2d0(ilower,iupper,g_gcw))
      REAL g1(FACE2d1(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1
c
c     Compute the face centered partial gradient of U.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            g0(i0,i1) = fac0*(U(i0,i1)-U(i0-1,i1)) + beta*v0(i0,i1)
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1+1
            g1(i1,i0) = fac1*(U(i0,i1)-U(i0,i1-1)) + beta*v1(i1,i0)
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
c     Uses centered differences to compute the side centered partial
c     gradient of a cell centered variable U.
c
c     This is a partial gradient in the sense that only the normal
c     component of the gradient is computed on each side.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosgrad2d(
     &     g0,g1,g_gcw,
     &     alpha,
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
      INTEGER g_gcw,U_gcw

      REAL alpha

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(SIDE2d0(ilower,iupper,g_gcw))
      REAL g1(SIDE2d1(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1
c
c     Compute the side centered partial gradient of U.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            g0(i0,i1) = fac0*(U(i0,i1)-U(i0-1,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            g1(i0,i1) = fac1*(U(i0,i1)-U(i0,i1-1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes g = alpha grad U + beta v.
c
c     Uses centered differences to compute the side centered partial
c     gradient of a cell centered variable U.
c
c     This is a partial gradient in the sense that only the normal
c     component of the gradient is computed on each side.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosgradadd2d(
     &     g0,g1,g_gcw,
     &     alpha,
     &     U,U_gcw,
     &     beta,
     &     v0,v1,v_gcw,
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
      INTEGER g_gcw,U_gcw,v_gcw

      REAL alpha

      REAL U(CELL2d(ilower,iupper,U_gcw))

      REAL beta

      REAL v0(SIDE2d0(ilower,iupper,v_gcw))
      REAL v1(SIDE2d1(ilower,iupper,v_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL g0(SIDE2d0(ilower,iupper,g_gcw))
      REAL g1(SIDE2d1(ilower,iupper,g_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1
c
c     Compute the side centered partial gradient of U.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            g0(i0,i1) = fac0*(U(i0,i1)-U(i0-1,i1)) + beta*v0(i0,i1)
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            g1(i0,i1) = fac1*(U(i0,i1)-U(i0,i1-1)) + beta*v1(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
