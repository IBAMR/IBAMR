c
c     Routines to compute level set operations on patches.
c
c     Created on 28 Sep 2017 by Nishant Nangia and Amneet Bhalla
c
c     Copyright (c) 2002-2014, Nishant Nangia and Amneet Bhalla
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
c     Compute the Godunov Hamiltonian.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function HG(a,b,c,d,e,f,sgn)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
      REAL a,b,c,d,e,f,sgn
      REAL am,ap,bm,bp,cm,cp,dm,dp,em,ep,fm,fp
      if (sgn .ge. zero) then
        am = dmin1(a,zero)
        bp = dmax1(b,zero)
        cm = dmin1(c,zero)
        dp = dmax1(d,zero)
        em = dmin1(e,zero)
        fp = dmax1(f,zero)
        HG = sqrt(dmax1(am**two,bp**two) 
     &          + dmax1(cm**two,dp**two)
     &          + dmax1(em**two,fp**two))
      else
        ap = dmax1(a,zero)
        bm = dmin1(b,zero)
        cp = dmax1(c,zero)
        dm = dmin1(d,zero)
        ep = dmax1(e,zero)
        fm = dmin1(f,zero)
        HG = sqrt(dmax1(ap**two,bm**two)
     &          + dmax1(cp**two,dm**two)
     &          + dmax1(ep**two,fm**two))
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out first-order accurate fast sweeping algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fastsweep1storder3d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dlower0,dupper0,
     &     dlower1,dupper1,
     &     dlower2,dupper2,
     &     dx,
     &     patch_touches_bdry,
     &     consider_bdry_wall)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER dlower0,dupper0
      INTEGER dlower1,dupper1
      INTEGER dlower2,dupper2
      INTEGER U_gcw
      INTEGER patch_touches_bdry
      INTEGER consider_bdry_wall

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    U_wall_coef
      REAL    h_wall_coef

      if (consider_bdry_wall .eq. 1) then
        U_wall_coef   = zero
        h_wall_coef   = half
      else
        U_wall_coef   = one
        h_wall_coef   = one
      endif
     
      
c     Do the eight sweeping directions.
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)

            enddo
         enddo
      enddo
      
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      do i2 = iupper2,ilower2,-1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      do i2 = iupper2,ilower2,-1
         do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      do i2 = iupper2,ilower2,-1
         do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      do i2 = iupper2,ilower2,-1
         do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
              call evalsweep1storder3d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 ilower2,iupper2,
     &                                 i0,i1,i2,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dlower2,dupper2,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
            enddo
         enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out single first order sweep
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine evalsweep1storder3d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     dlower0,dupper0,
     &     dlower1,dupper1,
     &     dlower2,dupper2,
     &     dx,
     &     patch_touches_bdry,
     &     U_wall_coef,h_wall_coef)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER dlower0,dupper0
      INTEGER dlower1,dupper1
      INTEGER dlower2,dupper2
      INTEGER U_gcw
      INTEGER patch_touches_bdry
      REAL    U_wall_coef,h_wall_coef

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      REAL    sgn
      INTEGER i0,i1,i2
      REAL    a,b,c
      REAL    a1,a2,a3
      REAL    hx,hy,hz
      REAL    h1,h2,h3
      REAL    Q,R,S
      REAL    dtil,dbar
      
c     Carry out a single sweep
      if (U(i0,i1,i2) .eq. zero) then
        sgn = zero
      else
        sgn = sign(one,U(i0,i1,i2))
      endif

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      a  = sgn*dmin1(sgn*U(i0-1,i1,i2),sgn*U(i0+1,i1,i2))
      b  = sgn*dmin1(sgn*U(i0,i1-1,i2),sgn*U(i0,i1+1,i2))
      c  = sgn*dmin1(sgn*U(i0,i1,i2-1),sgn*U(i0,i1,i2+1))

c     Take care of physical boundaries.
      if (patch_touches_bdry .eq. 1) then
        if (i0 .eq. dlower0) then 
          a  = U(i0+1,i1,i2)*U_wall_coef   
          hx = hx*h_wall_coef                
        elseif (i0 .eq. dupper0) then
          a  = U(i0-1,i1,i2)*U_wall_coef
          hx = hx*h_wall_coef
        endif
        if (i1 .eq. dlower1) then 
          b  = U(i0,i1+1,i2)*U_wall_coef   
          hy = hy*h_wall_coef                
        elseif (i1 .eq. dupper1) then
          b  = U(i0,i1-1,i2)*U_wall_coef
          hy = hy*h_wall_coef
        endif
        if (i2 .eq. dlower2) then
          c = U(i0,i1,i2+1)*U_wall_coef
          hz = hz*h_wall_coef
        elseif (i2 .eq. dupper2) then
          c = U(i0,i1,i2-1)*U_wall_coef
          hz = hz*h_wall_coef
        endif
      endif

c     Additional sorting step for 3D
      if (sgn*a .le. sgn*b) then
        if (sgn*a .le. sgn*c) then
          if(sgn*b .le. sgn*c) then
            a1 = a; a2 = b; a3 = c
            h1 = hx; h2 = hy; h3 = hz
          else
            a1 = a; a2 = c; a3 = b
            h1 = hx; h2 = hz; h3 = hy
          endif
        else
          a1 = c; a2 = a; a3 = b
          h1 = hz; h2 = hx; h3 = hy
        endif
      else 
        if (sgn*b .le. sgn*c) then
          if (sgn*a .le. sgn*c) then
            a1 = b; a2 = a; a3 = c
            h1 = hy; h2 = hx; h3 = hz
          else
            a1 = b; a2 = c; a3 = a
            h1 = hy; h2 = hz; h3 = hx
          endif
        else
          a1 = c; a2 = b; a3 = a
          h1 = hz; h2 = hy; h3 = hx
        endif
      endif

c     Algorithm to find dbar
      dtil = a1 + sgn*h1
      if (sgn*dtil .le. sgn*a2) then
        dbar = dtil
      else
        Q = h1*h1 + h2*h2
        R = -2.d0*(h2*h2*a1 + h1*h1*a2)
        S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
        dtil = (-R + sgn*sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
      endif
    
      if (sgn*dtil .lt. sgn*a3) then
        dbar = dtil
      else
        Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
        R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
        S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
        dtil = (-R + sgn*sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
        dbar = dtil
      endif

      U(i0,i1,i2) = sgn*dmin1(sgn*U(i0,i1,i2),sgn*dbar)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out first order relaxation scheme using Gauss Seidel updates
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxation1storder3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     dir)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw
      INTEGER dir

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2

      if (dir .eq. 0) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 3) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 4) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 5) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 6) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 7) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax1storder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      endif

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute relaxation solution at a given grid cell
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine evalrelax1storder3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL HG

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    hx,hy,hz
      REAL    sgn
      REAL    dt
      REAL    G
      REAL    Dxp,Dxm
      REAL    Dyp,Dym
      REAL    Dzp,Dzm

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      dt = third*dmin1(hx,hy,hz)
      if (V(i0,i1,i2) .eq. zero) then
         sgn = zero
      else
         sgn = sign(one,V(i0,i1,i2))
      endif

      Dxm = one/hx*(U(i0,i1,i2)-U(i0-1,i1,i2))
      Dxp = one/hx*(U(i0+1,i1,i2)-U(i0,i1,i2))
      Dym = one/hy*(U(i0,i1,i2)-U(i0,i1-1,i2))
      Dyp = one/hy*(U(i0,i1+1,i2)-U(i0,i1,i2))
      Dzm = one/hz*(U(i0,i1,i2)-U(i0,i1,i2-1))
      Dzp = one/hz*(U(i0,i1,i2+1)-U(i0,i1,i2))
      G = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)

      U(i0,i1,i2) = U(i0,i1,i2) - dt*sgn*(G-one)

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out third order relaxation scheme using Gauss Seidel updates
c     NOTE: this scheme between third and fourth s
c     order near the interface and second order everywhere else
c
c     Uses second order WENO for spatial discretization with a subcell
c     fix near the interface
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxation3rdorder3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     dir)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw
      INTEGER dir

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2

      if (dir .eq. 0) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 3) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 4) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 5) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 6) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      elseif (dir .eq. 7) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorder3d(U,U_gcw,V,V_gcw,
     &                                   ilower0,iupper0,
     &                                   ilower1,iupper1,
     &                                   ilower2,iupper2,
     &                                   i0,i1,i2,dx)
            enddo
          enddo
        enddo
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out single third order sweep
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine evalrelax3rdorder3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL minmod, HG

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    hx,hxp,hxm
      REAL    hy,hyp,hym
      REAL    hz,hzp,hzm
      REAL    Dxm,Dxp,Dym,Dyp,Dzm,Dzp
      REAL    Dxx,Dxxp,Dxxm
      REAL    Dyy,Dyyp,Dyym
      REAL    Dzz,Dzzp,Dzzm
      REAL    Dxx0,Dyy0,Dzz0
      REAL    H,dt,sgn,cfl,eps,D,diff

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      cfl = 0.3d0
      eps = 1.d-10

      if (V(i0,i1,i2) .eq. zero) then
         sgn = zero
      else
         sgn = sign(one,V(i0,i1,i2))
      endif

c     Compute all the required finite differences
      Dxx  = (U(i0-1,i1,i2) - two*U(i0,i1,i2) + U(i0+1,i1,i2))/(hx**two)
      Dxxp = (U(i0,i1,i2) - two*U(i0+1,i1,i2) + U(i0+2,i1,i2))/(hx**two)
      Dxxm = (U(i0-2,i1,i2) - two*U(i0-1,i1,i2) + U(i0,i1,i2))/(hx**two)
      Dyy  = (U(i0,i1-1,i2) - two*U(i0,i1,i2) + U(i0,i1+1,i2))/(hy**two)
      Dyyp = (U(i0,i1,i2) - two*U(i0,i1+1,i2) + U(i0,i1+2,i2))/(hy**two)
      Dyym = (U(i0,i1-2,i2) - two*U(i0,i1-1,i2) + U(i0,i1,i2))/(hy**two)
      Dzz  = (U(i0,i1,i2-1) - two*U(i0,i1,i2) + U(i0,i1,i2+1))/(hz**two)
      Dzzp = (U(i0,i1,i2) - two*U(i0,i1,i2+1) + U(i0,i1,i2+2))/(hz**two)
      Dzzm = (U(i0,i1,i2-2) - two*U(i0,i1,i2-1) + U(i0,i1,i2))/(hz**two)

c     Set dummy values for hxp,hxm,hyp,hym
      hxp = 1.d12;hxm = 1.d12
      hyp = 1.d12;hym = 1.d12
      hzp = 1.d12;hzm = 1.d12

c     Compute ENO differences with subcell fix
      if (V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero) then
        Dxx0 = minmod(V(i0-1,i1,i2)-two*V(i0,i1,i2)+V(i0+1,i1,i2),
     &                V(i0,i1,i2)-two*V(i0+1,i1,i2)+V(i0+2,i1,i2))
        diff = V(i0,i1,i2)-V(i0+1,i1,i2)
        if (abs(Dxx0) .gt. eps) then
          D = (Dxx0/two-V(i0,i1,i2)-V(i0+1,i1,i2))**two
     &        -four*V(i0,i1,i2)*V(i0+1,i1,i2)
          hxp = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
        else
          hxp = hx*V(i0,i1,i2)/diff
        endif
        Dxp = (zero-U(i0,i1,i2))/hxp - hxp/two*minmod(Dxx,Dxxp)
      else
        Dxp = (U(i0+1,i1,i2)-U(i0,i1,i2))/hx - hx/two*minmod(Dxx,Dxxp)
      endif

      if (V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero) then
        Dxx0 = minmod(V(i0-1,i1,i2)-two*V(i0,i1,i2)+V(i0+1,i1,i2),
     &                V(i0,i1,i2)-two*V(i0-1,i1,i2)+V(i0-2,i1,i2))
        diff = V(i0,i1,i2)-V(i0-1,i1,i2)
        if (abs(Dxx0) .gt. eps) then
          D = (Dxx0/two-V(i0,i1,i2)-V(i0-1,i1,i2))**two
     &        -four*V(i0,i1,i2)*V(i0-1,i1,i2)
          hxm = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
        else
          hxm = hx*V(i0,i1,i2)/diff
        endif
        Dxm = (U(i0,i1,i2)-zero)/hxm + hxm/two*minmod(Dxx,Dxxm)
      else
        Dxm = (U(i0,i1,i2)-U(i0-1,i1,i2))/hx + hx/two*minmod(Dxx,Dxxm)
      endif

      if (V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero) then
        Dyy0 = minmod(V(i0,i1-1,i2)-two*V(i0,i1,i2)+V(i0,i1+1,i2),
     &                V(i0,i1,i2)-two*V(i0,i1+1,i2)+V(i0,i1+2,i2))
        diff = V(i0,i1,i2)-V(i0,i1+1,i2)
        if (abs(Dyy0) .gt. eps) then
          D = (Dyy0/two-V(i0,i1,i2)-V(i0,i1+1,i2))**two
     &        -four*V(i0,i1,i2)*V(i0,i1+1,i2)
          hyp = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
        else
          hyp = hy*V(i0,i1,i2)/diff
        endif
        Dyp = (zero-U(i0,i1,i2))/hyp - hyp/two*minmod(Dyy,Dyyp)
      else
        Dyp = (U(i0,i1+1,i2)-U(i0,i1,i2))/hy - hy/two*minmod(Dyy,Dyyp)
      endif

      if (V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero) then
        Dyy0 = minmod(V(i0,i1-1,i2)-two*V(i0,i1,i2)+V(i0,i1+1,i2),
     &                V(i0,i1,i2)-two*V(i0,i1-1,i2)+V(i0,i1-2,i2))
        diff = V(i0,i1,i2)-V(i0,i1-1,i2)
        if (abs(Dyy0) .gt. eps) then
          D = (Dyy0/two-V(i0,i1,i2)-V(i0,i1-1,i2))**two
     &        -four*V(i0,i1,i2)*V(i0,i1-1,i2)
          hym = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
        else
          hym = hy*V(i0,i1,i2)/diff
        endif
        Dym = (U(i0,i1,i2)-zero)/hym + hym/two*minmod(Dyy,Dyym)
      else
        Dym = (U(i0,i1,i2)-U(i0,i1-1,i2))/hy + hy/two*minmod(Dyy,Dyym)
      endif

      if (V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero) then
        Dzz0 = minmod(V(i0,i1,i2-1)-two*V(i0,i1,i2)+V(i0,i1,i2+1),
     &                V(i0,i1,i2)-two*V(i0,i1,i2+1)+V(i0,i1,i2+2))
        diff = V(i0,i1,i2)-V(i0,i1,i2+1)
        if (abs(Dzz0) .gt. eps) then
          D = (Dzz0/two-V(i0,i1,i2)-V(i0,i1,i2+1))**two
     &        -four*V(i0,i1,i2)*V(i0,i1,i2+1)
          hzp = hz*(half + (diff-sign(one,diff)*sqrt(D))/Dzz0)
        else
          hzp = hz*V(i0,i1,i2)/diff
        endif
        Dzp = (zero-U(i0,i1,i2))/hzp - hzp/two*minmod(Dzz,Dzzp)
      else
        Dzp = (U(i0,i1,i2+1)-U(i0,i1,i2))/hz - hz/two*minmod(Dzz,Dzzp)
      endif

      if (V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero) then
        Dzz0 = minmod(V(i0,i1,i2-1)-two*V(i0,i1,i2)+V(i0,i1,i2+1),
     &                V(i0,i1,i2)-two*V(i0,i1,i2-1)+V(i0,i1,i2-2))
        diff = V(i0,i1,i2)-V(i0,i1,i2-1)
        if (abs(Dzz0) .gt. eps) then
          D = (Dzz0/two-V(i0,i1,i2)-V(i0,i1,i2-1))**two
     &        -four*V(i0,i1,i2)*V(i0,i1,i2-1)
          hzm = hz*(half + (diff-sign(one,diff)*sqrt(D))/Dzz0)
        else
          hzm = hz*V(i0,i1,i2)/diff
        endif
        Dzm = (U(i0,i1,i2)-zero)/hzm + hzm/two*minmod(Dzz,Dzzm)
      else
        Dzm = (U(i0,i1,i2)-U(i0,i1,i2-1))/hz + hz/two*minmod(Dzz,Dzzm)
      endif

      H = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)
      dt = cfl*dmin1(hx,hy,hz,hxp,hxm,hyp,hym,hzp,hzm)

      U(i0,i1,i2) = U(i0,i1,i2) - dt*sgn*(H-one)

      return
      end

 