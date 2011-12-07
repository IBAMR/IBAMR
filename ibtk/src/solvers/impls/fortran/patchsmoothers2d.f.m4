c
c     Patch smoothers for use in multigrid preconditioners.
c
c     Created on 25 Aug 2010 by Boyce Griffith
c
c     Copyright (c) 2002-2010, Boyce Griffith
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
c        * Neither the name of New York University nor the names of its
c          contributors may be used to endorse or promote products
c          derived from this software without specific prior written
c          permission.
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
c     Gauss-Seidel sweeps for F = alpha div grad U + beta U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine gssmooth2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     sweeps)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw
      INTEGER sweeps

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,l
      REAL    fac0,fac1,fac
c
c     Perform one or more Gauss-Seidel sweeps.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do l = 1,sweeps
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Red-black Gauss-Seidel sweeps for F = alpha div grad U + beta U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rbgssmooth2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     sweeps)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw
      INTEGER sweeps

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,l,m
      REAL    fac0,fac1,fac
c
c     Perform one or more red-black Gauss-Seidel sweeps.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do l = 1,sweeps
         do m = 0,1
            do i1 = ilower1,iupper1
               do i0 = ilower0,iupper0
                  if ( mod(i0+i1,2) .eq. m ) then
                     U(i0,i1) = fac*(
     &                    fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &                    fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &                    F(i0,i1))
                  endif
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
