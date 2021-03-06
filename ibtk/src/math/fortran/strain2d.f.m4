c ---------------------------------------------------------------------
c
c Copyright (c) 2016 - 2019 by the IBAMR developers
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
c     Computes E_diag    = (du0/dx0, du1/dx1).
c     Computes E_offDiag = 0.5*(du0/dx1 + du1/dx0)
c
c     Uses centered differences to compute the cell centered strain of a
c     side centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocstrain2d(
     &     E_diag,E_diag_gcw,
     &     E_off,E_off_gcw,
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
      INTEGER E_diag_gcw,E_off_gcw,u_gcw

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL E_diag(CELL2d(ilower,iupper,E_diag_gcw),0:NDIM-1)
      REAL E_off(CELL2d(ilower,iupper,E_off_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    du0_dx0,du1_dx1,du0_dx1,du1_dx0
      REAL    fac00,fac11,fac01,fac10
c
c     Compute the cell centered stress of u=(u0,u1).
c
      fac00 = 1.0d0/dx(0)
      fac11 = 1.0d0/dx(1)
      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            du0_dx0 = fac00*(u0(i0+1,i1) - u0(i0,i1))
            du1_dx1 = fac11*(u1(i0,i1+1) - u1(i0,i1))
            
            du0_dx1 = fac01*(
     &           +u0(i0  ,i1+1)+u0(i0+1,i1+1)
     &           -u0(i0  ,i1-1)-u0(i0+1,i1-1) )
            du1_dx0 = fac10*(
     &           +u1(i0+1,i1  )+u1(i0+1,i1+1)
     &           -u1(i0-1,i1  )-u1(i0-1,i1+1) )
     
            E_diag(i0,i1,0) = du0_dx0
            E_diag(i0,i1,1) = du1_dx1
            
            E_off(i0,i1) = 0.5d0*(du1_dx0+du0_dx1)
         enddo
      enddo
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
