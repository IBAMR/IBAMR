c
c     Routines to compute discrete strain rate on patches.
c
c     Created on 08 Dec 2016 by Nishant Nangia
c
c     Copyright (c) 2002-2017, Boyce Griffith
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
