c
c     Routines to compute level set operations on patches.
c
c     Created on 28 Sep 2017 by Nishant Nangia and Amneet Bhalla
c
c     Copyright (c) 2002-2017, Nishant Nangia and Amneet Bhalla
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
c     Compute the smoothed Heaviside
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function H_eps(x,eps)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
      REAL x,eps
      if (x .lt. -eps) then
        H_eps = zero
      else if (abs(x) .le. eps) then
        H_eps = half*(one + x/eps + sin(pi*x/eps)/pi)
      else
        H_eps = one
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the smoothed delta function
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function D_eps(x,eps)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
      REAL x,eps
      if (abs(x) .le. eps) then
        D_eps = one/(two*eps)*(one + cos(pi*x/eps))
      else
        D_eps = zero
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the smoothed sgn function
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function S_eps(x,eps)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
      REAL x,eps
c
c     Prevent compiler warning about unused variables.
c
      eps = eps

C       S_eps = two*H_eps(x,eps) - one
C       S_eps = x/sqrt(x**2+eps**2)

      S_eps = sign(one,x)
      if (x.eq.zero) then
        S_eps = zero
      endif

      return
      end
c
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the WENO5 approximation based on Jiang and Peng 2000
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL function WENO5(Q)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
      REAL Q(-2:1)
      REAL a,b,c,d
      REAL w0,w2
      REAL eps,a0,a1,a2
      REAL IS0,IS1,IS2
      a = Q(-2); b = Q(-1); c = Q(0); d = Q(1)

      eps = 1.d-6
      IS0 = 13.d0*(a-b)**2+three*(a-three*b)**2
      IS1 = 13.d0*(b-c)**2+three*(b+c)**2
      IS2 = 13.d0*(c-d)**2+three*(three*c-d)**2
      a0 = one/(eps+IS0)**2
      a1 = 6.d0/(eps+IS1)**2
      a2 = three/(eps+IS2)**2
      w0 = a0/(a0+a1+a2); w2 = a2/(a0+a1+a2)
      WENO5 = third*w0*(a-two*b+c)+sixth*(w2-half)*(b-two*c+d)
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
     &     touches_wall_loc_idx)
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

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
      INTEGER touches_wall_loc_idx(0:2*NDIM - 1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      
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
     &                                 touches_wall_loc_idx)

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
     &                                 touches_wall_loc_idx)
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
     &                                 touches_wall_loc_idx)
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
     &                                 touches_wall_loc_idx)
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
     &                                 touches_wall_loc_idx)
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
     &                                 touches_wall_loc_idx)
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
     &                                 touches_wall_loc_idx)
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
     &                                 touches_wall_loc_idx)
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
     &     touches_wall_loc_idx)
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

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
      INTEGER touches_wall_loc_idx(0:2*NDIM - 1)
      
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
c     The grid spacing to the boundary will be h/2
c     The distance value imposed at the boundary should be zero
      if (patch_touches_bdry .eq. 1) then
        if (i0 .eq. dlower0 .and. 
     &      touches_wall_loc_idx(0) .eq. 1) then 
          a  = zero   
          hx = hx*half                
        elseif (i0 .eq. dupper0 .and. 
     &      touches_wall_loc_idx(1) .eq. 1) then
          a  = zero
          hx = hx*half
        endif
        if (i1 .eq. dlower1 .and. 
     &      touches_wall_loc_idx(2) .eq. 1) then 
          b  = zero  
          hy = hy*half               
        elseif (i1 .eq. dupper1 .and. 
     &      touches_wall_loc_idx(3) .eq. 1) then
          b  = zero
          hy = hy*half
        endif
        if (i2 .eq. dlower2 .and. 
     &      touches_wall_loc_idx(4) .eq. 1) then
          c = zero
          hz = hz*half
        elseif (i2 .eq. dupper2 .and. 
     &      touches_wall_loc_idx(5) .eq. 1) then
          c = zero
          hz = hz*half
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
      subroutine relaxationls1storder3d(
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
c     NOTE: this scheme is between third and fourth
c     order near the interface and second order everywhere else
c
c     Uses second order ENO for spatial discretization with a subcell
c     fix near the interface
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxationls3rdordereno3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     dir,use_subcell,
     &     use_sign_fix)
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
      INTEGER use_subcell,use_sign_fix

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
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 3) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 4) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 5) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 6) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 7) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdordereno3d(U,U_gcw,V,V_gcw,
     &                                      ilower0,iupper0,
     &                                      ilower1,iupper1,
     &                                      ilower2,iupper2,
     &                                      i0,i1,i2,dx,
     &                                      use_subcell,use_sign_fix)
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
      subroutine evalrelax3rdordereno3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     dx,
     &     use_subcell,
     &     use_sign_fix)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL minmod, HG, S_eps

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw
      INTEGER use_subcell,use_sign_fix

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
      REAL    hmin

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      cfl = 0.3d0
      eps = 1.d-10

      hmin = dmin1(hx,hy,hz)
      sgn = S_eps(V(i0,i1,i2),hmin)

c     Sign fix
      if (use_sign_fix .ne. 0) then
        if (V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0+1,i1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0-1,i1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1+1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1-1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1,i2+1))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1,i2-1))) then
          sgn = zero
        endif
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
      if (use_subcell .ne. 0 .and. 
     &    V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero) then
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

      if (use_subcell .ne. 0 .and.
     &    V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero) then
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

      if (use_subcell .ne. 0 .and.
     &    V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero) then
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

      if (use_subcell .ne. 0 .and.
     &    V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero) then
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

      if (use_subcell .ne. 0 .and.
     &    V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero) then
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

      if (use_subcell .ne. 0 .and.
     &    V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero) then
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

      if (dt .gt. zero) then
        U(i0,i1,i2) = U(i0,i1,i2) - dt*sgn*(H-one)
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Godunov Hamiltonian of the indicator field |grad phi_0|
c     
c     Uses second order ENO for spatial discretization with a subcell
c     fix near the interface
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine godunovhamiltonianeno3d(
     &     H,H_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     use_subcell)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL minmod, HG, S_eps

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER H_gcw,V_gcw
      INTEGER use_subcell

c
c     Input/Output.
c
      REAL H(CELL3d(ilower,iupper,H_gcw))
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
      REAL    sgn,eps,D,diff
      REAL    hmin

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      eps = 1.d-10
      hmin = dmin1(hx,hy,hz)
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
          sgn = S_eps(V(i0,i1,i2),hmin)

c         Compute all the required finite differences
          Dxx  = (V(i0-1,i1,i2) - two*V(i0,i1,i2) + V(i0+1,i1,i2))
     &            /(hx**two)
          Dxxp = (V(i0,i1,i2) - two*V(i0+1,i1,i2) + V(i0+2,i1,i2))
     &            /(hx**two)
          Dxxm = (V(i0-2,i1,i2) - two*V(i0-1,i1,i2) + V(i0,i1,i2))
     &            /(hx**two)
          Dyy  = (V(i0,i1-1,i2) - two*V(i0,i1,i2) + V(i0,i1+1,i2))
     &            /(hy**two)
          Dyyp = (V(i0,i1,i2) - two*V(i0,i1+1,i2) + V(i0,i1+2,i2))
     &            /(hy**two)
          Dyym = (V(i0,i1-2,i2) - two*V(i0,i1-1,i2) + V(i0,i1,i2))
     &            /(hy**two)
          Dzz  = (V(i0,i1,i2-1) - two*V(i0,i1,i2) + V(i0,i1,i2+1))
     &           /(hz**two)
          Dzzp = (V(i0,i1,i2) - two*V(i0,i1,i2+1) + V(i0,i1,i2+2))
     &           /(hz**two)
          Dzzm = (V(i0,i1,i2-2) - two*V(i0,i1,i2-1) + V(i0,i1,i2))
     &           /(hz**two)

c         Set dummy values for hxp,hxm,hyp,hym
          hxp = 1.d12;hxm = 1.d12
          hyp = 1.d12;hym = 1.d12
          hzp = 1.d12;hzm = 1.d12

c         Compute ENO differences with subcell fix
          if (use_subcell .ne. 0 .and. 
     &        V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero) then
            Dxx0 = minmod(V(i0-1,i1,i2)-two*V(i0,i1,i2)+V(i0+1,i1,i2),
     &                    V(i0,i1,i2)-two*V(i0+1,i1,i2)+V(i0+2,i1,i2))
            diff = V(i0,i1,i2)-V(i0+1,i1,i2)
            if (abs(Dxx0) .gt. eps) then
              D = (Dxx0/two-V(i0,i1,i2)-V(i0+1,i1,i2))**two
     &            -four*V(i0,i1,i2)*V(i0+1,i1,i2)
              hxp = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
            else
              hxp = hx*V(i0,i1,i2)/diff
            endif
            hxp = dmax1(hxp,sqrt(smallr))
            Dxp = (zero-V(i0,i1,i2))/hxp - hxp/two*minmod(Dxx,Dxxp)
          else
            Dxp = (V(i0+1,i1,i2)-V(i0,i1,i2))/hx-hx/two*minmod(Dxx,Dxxp)
          endif

          if (use_subcell .ne. 0 .and.
     &        V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero) then
           Dxx0 = minmod(V(i0-1,i1,i2)-two*V(i0,i1,i2)+V(i0+1,i1,i2),
     &                    V(i0,i1,i2)-two*V(i0-1,i1,i2)+V(i0-2,i1,i2))
            diff = V(i0,i1,i2)-V(i0-1,i1,i2)
            if (abs(Dxx0) .gt. eps) then
              D = (Dxx0/two-V(i0,i1,i2)-V(i0-1,i1,i2))**two
     &            -four*V(i0,i1,i2)*V(i0-1,i1,i2)
              hxm = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
            else
              hxm = hx*V(i0,i1,i2)/diff
            endif
            hxm = dmax1(hxm,sqrt(smallr))
            Dxm = (V(i0,i1,i2)-zero)/hxm + hxm/two*minmod(Dxx,Dxxm)
          else
            Dxm = (V(i0,i1,i2)-V(i0-1,i1,i2))/hx+hx/two*minmod(Dxx,Dxxm)
          endif

          if (use_subcell .ne. 0 .and.
     &        V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero) then
            Dyy0 = minmod(V(i0,i1-1,i2)-two*V(i0,i1,i2)+V(i0,i1+1,i2),
     &                    V(i0,i1,i2)-two*V(i0,i1+1,i2)+V(i0,i1+2,i2))
            diff = V(i0,i1,i2)-V(i0,i1+1,i2)
            if (abs(Dyy0) .gt. eps) then
              D = (Dyy0/two-V(i0,i1,i2)-V(i0,i1+1,i2))**two
     &            -four*V(i0,i1,i2)*V(i0,i1+1,i2)
              hyp = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
            else
              hyp = hy*V(i0,i1,i2)/diff
            endif
            hyp = dmax1(hyp,sqrt(smallr))
            Dyp = (zero-V(i0,i1,i2))/hyp - hyp/two*minmod(Dyy,Dyyp)
          else
            Dyp = (V(i0,i1+1,i2)-V(i0,i1,i2))/hy-hy/two*minmod(Dyy,Dyyp)
          endif

          if (use_subcell .ne. 0 .and.
     &        V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero) then
            Dyy0 = minmod(V(i0,i1-1,i2)-two*V(i0,i1,i2)+V(i0,i1+1,i2),
     &                    V(i0,i1,i2)-two*V(i0,i1-1,i2)+V(i0,i1-2,i2))
            diff = V(i0,i1,i2)-V(i0,i1-1,i2)
            if (abs(Dyy0) .gt. eps) then
              D = (Dyy0/two-V(i0,i1,i2)-V(i0,i1-1,i2))**two
     &            -four*V(i0,i1,i2)*V(i0,i1-1,i2)
              hym = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
            else
              hym = hy*V(i0,i1,i2)/diff
            endif
            hym = dmax1(hym,sqrt(smallr))
            Dym = (V(i0,i1,i2)-zero)/hym + hym/two*minmod(Dyy,Dyym)
          else
            Dym = (V(i0,i1,i2)-V(i0,i1-1,i2))/hy+hy/two*minmod(Dyy,Dyym)
          endif

          if (use_subcell .ne. 0 .and.
     &        V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero) then
            Dzz0 = minmod(V(i0,i1,i2-1)-two*V(i0,i1,i2)+V(i0,i1,i2+1),
     &                    V(i0,i1,i2)-two*V(i0,i1,i2+1)+V(i0,i1,i2+2))
            diff = V(i0,i1,i2)-V(i0,i1,i2+1)
            if (abs(Dzz0) .gt. eps) then
              D = (Dzz0/two-V(i0,i1,i2)-V(i0,i1,i2+1))**two
     &            -four*V(i0,i1,i2)*V(i0,i1,i2+1)
              hzp = hz*(half + (diff-sign(one,diff)*sqrt(D))/Dzz0)
            else
              hzp = hz*V(i0,i1,i2)/diff
            endif
            hzp = dmax1(hzp,sqrt(smallr))
            Dzp = (zero-V(i0,i1,i2))/hzp - hzp/two*minmod(Dzz,Dzzp)
          else
            Dzp = (V(i0,i1,i2+1)-V(i0,i1,i2))/hz-hz/two*minmod(Dzz,Dzzp)
          endif

          if (use_subcell .ne. 0 .and.
     &        V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero) then
            Dzz0 = minmod(V(i0,i1,i2-1)-two*V(i0,i1,i2)+V(i0,i1,i2+1),
     &                    V(i0,i1,i2)-two*V(i0,i1,i2-1)+V(i0,i1,i2-2))
            diff = V(i0,i1,i2)-V(i0,i1,i2-1)
            if (abs(Dzz0) .gt. eps) then
              D = (Dzz0/two-V(i0,i1,i2)-V(i0,i1,i2-1))**two
     &            -four*V(i0,i1,i2)*V(i0,i1,i2-1)
              hzm = hz*(half + (diff-sign(one,diff)*sqrt(D))/Dzz0)
            else
              hzm = hz*V(i0,i1,i2)/diff
            endif
            hzm = dmax1(hzm,sqrt(smallr))
            Dzm = (V(i0,i1,i2)-zero)/hzm + hzm/two*minmod(Dzz,Dzzm)
          else
            Dzm = (V(i0,i1,i2)-V(i0,i1,i2-1))/hz+hz/two*minmod(Dzz,Dzzm)
          endif

          H(i0,i1,i2) = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)
          enddo
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out third order relaxation scheme using Gauss Seidel updates
c     
c     Uses third order WENO for spatial discretization with a subcell
c     fix near the interface
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxationls3rdorderweno3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     dir,
     &     use_subcell,
     &     use_sign_fix)
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
      INTEGER use_subcell,use_sign_fix

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
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 3) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 4) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 5) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 6) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 7) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax3rdorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_subcell,use_sign_fix)
            enddo
          enddo
        enddo
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out single third order sweep using a WENO stencil with
c     subcell fix
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine evalrelax3rdorderweno3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     dx,
     &     use_subcell,
     &     use_sign_fix)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL minmod, HG, S_eps

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,V_gcw
      INTEGER use_subcell,use_sign_fix

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
      REAL    hmin
      REAL    Dxm,Dxp,Dym,Dyp,Dzm,Dzp
      REAL    Dxc,Dyc,Dzc
      REAL    Dxl,Dxr
      REAL    Dyb,Dyt
      REAL    Dzd,Dzu
      REAL    rxm,wxm
      REAL    rxp,wxp
      REAL    rym,wym
      REAL    ryp,wyp
      REAL    rzm,wzm
      REAL    rzp,wzp
      REAL    h1,h2
      REAL    Dxx0,Dyy0,Dzz0
      REAL    H,dt,sgn,cfl,eps,D,diff

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      cfl = 0.3d0
      eps = 1.d-10

      hmin = dmin1(hx,hy,hz)
      sgn = S_eps(V(i0,i1,i2),hmin)

c     Sign fix
      if (use_sign_fix .ne. 0) then
        if (V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0+1,i1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0-1,i1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1+1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1-1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1,i2+1))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1,i2-1))) then
          sgn = zero
        endif
      endif

c     Compute all the required finite differences
      Dxc  = (U(i0+1,i1,i2) - U(i0-1,i1,i2))/(two*hx) 
      Dyc  = (U(i0,i1+1,i2) - U(i0,i1-1,i2))/(two*hy)
      Dzc  = (U(i0,i1,i2+1) - U(i0,i1,i2-1))/(two*hz)
      Dxl  = (three*U(i0,i1,i2) - four*U(i0-1,i1,i2) + U(i0-2,i1,i2))
     &       /(two*hx)
      Dxr = (-three*U(i0,i1,i2) + four*U(i0+1,i1,i2) - U(i0+2,i1,i2))
     &      /(two*hx)
      Dyb  = (three*U(i0,i1,i2) - four*U(i0,i1-1,i2) + U(i0,i1-2,i2))
     &       /(two*hy)
      Dyt = (-three*U(i0,i1,i2) + four*U(i0,i1+1,i2) - U(i0,i1+2,i2))
     &      /(two*hy)
      Dzd  = (three*U(i0,i1,i2) - four*U(i0,i1,i2-1) + U(i0,i1,i2-2))
     &       /(two*hz)
      Dzu = (-three*U(i0,i1,i2) + four*U(i0,i1,i2+1) - U(i0,i1,i2+2))
     &      /(two*hz)

      rxm = (eps + (U(i0,i1,i2) -two*U(i0-1,i1,i2) + U(i0-2,i1,i2))**two
     &     )/(eps + (U(i0+1,i1,i2) -two*U(i0,i1,i2) + U(i0-1,i1,i2))**
     &     two)
      wxm = one/(one + two*rxm**two) 

      rxp = (eps + (U(i0+2,i1,i2) -two*U(i0+1,i1,i2) + U(i0,i1,i2))**
     &     two)/(eps + (U(i0+1,i1,i2) -two*U(i0,i1,i2) + U(i0-1,i1,i2))
     &     **two)
      wxp = one/(one + two*rxp**two)

      rym = (eps + (U(i0,i1,i2) -two*U(i0,i1-1,i2) + U(i0,i1-2,i2))**
     &     two)/(eps + (U(i0,i1+1,i2) -two*U(i0,i1,i2) + U(i0,i1-1,i2))
     &     **two)
      wym = one/(one + two*rym**two)

      ryp = (eps + (U(i0,i1+2,i2) -two*U(i0,i1+1,i2) + U(i0,i1,i2))**
     &     two)/(eps + (U(i0,i1+1,i2) -two*U(i0,i1,i2) + U(i0,i1-1,i2))
     &     **two)
      wyp = one/(one + two*ryp**two)

      rzm = (eps + (U(i0,i1,i2) -two*U(i0,i1,i2-1) + U(i0,i1,i2-2))**
     &     two)/(eps + (U(i0,i1,i2+1) -two*U(i0,i1,i2) + U(i0,i1,i2-1))
     &     **two)
      wzm = one/(one + two*rzm**two)

      rzp = (eps + (U(i0,i1,i2+2) -two*U(i0,i1,i2+1) + U(i0,i1,i2))**
     &     two)/(eps + (U(i0,i1,i2+1) -two*U(i0,i1,i2) + U(i0,i1,i2-1))
     &     **two)
      wzp = one/(one + two*rzp**two)
 

c     Set dummy values for hxp,hxm,hyp,hym
      hxp = 1.d12;hxm = 1.d12
      hyp = 1.d12;hym = 1.d12
      hzp = 1.d12;hzm = 1.d12

c     Compute ENO differences with subcell fix
      if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero) 
     &  then
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
        
        h1 = hx; h2 = hxp
        Dxc = (-U(i0,i1,i2)*h1**2 + 
     &        (-U(i0-1,i1,i2) + U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

        h1 = hxp; h2 = hx - hxp
        Dxr = (-U(i0+1,i1,i2)*h1**2 + two*(-U(i0,i1,i2))*h1*h2 + 
     &          (-U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2)) 
      endif

      if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero)
     &  then
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
       
        h1 = hxm; h2 = hx
        Dxc = ((-U(i0,i1,i2) + U(i0+1, i1,i2))*h1**2 + 
     &         (U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

        h1 = hx - hxm; h2 = hxm
        Dxl = (U(i0,i1,i2)*h1**2 + two*U(i0,i1,i2)*h1*h2 + 
     &         U(i0-1,i1,i2)*h2**2)/(h1*h2*(h1 + h2))
      endif

      if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero)
     &  then
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
 
        h1 = hy; h2 = hyp
        Dyc = (-U(i0,i1,i2)*h1**2 +
     &        (-U(i0,i1-1,i2) + U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

        h1 = hyp; h2 = hy - hyp
        Dyt = (-U(i0,i1+1,i2)*h1**2 + two*(-U(i0,i1,i2))*h1*h2 +
     &          (-U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))
      endif

      if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero)
     &  then
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

        h1 = hym; h2 = hy
        Dyc = ((-U(i0,i1,i2) + U(i0, i1+1,i2))*h1**2 +
     &         (U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

        h1 = hy - hym; h2 = hym
        Dyb = (U(i0,i1,i2)*h1**2 + two*U(i0,i1,i2)*h1*h2 +
     &         U(i0,i1-1,i2)*h2**2)/(h1*h2*(h1 + h2))  
      endif

      if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero)
     &  then
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
 
        h1 = hz; h2 = hzp
        Dzc = (-U(i0,i1,i2)*h1**2 +
     &        (-U(i0,i1,i2-1) + U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

        h1 = hzp; h2 = hz - hzp
        Dzu = (-U(i0,i1,i2+1)*h1**2 + two*(-U(i0,i1,i2))*h1*h2 +
     &          (-U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))
      endif

      if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero)
     &  then
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

        h1 = hzm; h2 = hz
        Dzc = ((-U(i0,i1,i2) + U(i0, i1,i2+1))*h1**2 +
     &         (U(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

        h1 = hz - hzm; h2 = hzm
        Dzd = (U(i0,i1,i2)*h1**2 + two*U(i0,i1,i2)*h1*h2 +
     &         U(i0,i1,i2-1)*h2**2)/(h1*h2*(h1 + h2))  
      endif

c     Compute first order derivatives

      Dxm = (one - wxm)*Dxc + wxm*Dxl
      Dxp = (one - wxp)*Dxc + wxp*Dxr

      Dym = (one - wym)*Dyc + wym*Dyb
      Dyp = (one - wyp)*Dyc + wyp*Dyt

      Dzm = (one - wzm)*Dzc + wzm*Dzd
      Dzp = (one - wzp)*Dzc + wzp*Dzu

      H = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)
      dt = cfl*dmin1(hx,hy,hz,
     &               hxp,hxm,hyp,hym,hzp,hzm,
     &               abs(hx-hxm),abs(hy-hym),abs(hz-hzm),
     &               abs(hx-hxp),abs(hy-hyp),abs(hz-hzp))

      if (dt .gt. zero) then
        U(i0,i1,i2) = U(i0,i1,i2) - dt*sgn*(H-one)
      endif

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Godunov Hamiltonian of the indicator field |grad phi_0|
c     
c     Uses second order WENO for spatial discretization with a subcell
c     fix near the interface
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine godunovhamiltonianweno3d(
     &     H,H_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     use_subcell)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL minmod, HG, S_eps

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER H_gcw,V_gcw
      INTEGER use_subcell

c
c     Input/Output.
c
      REAL H(CELL3d(ilower,iupper,H_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    hx,hxp,hxm
      REAL    hy,hyp,hym
      REAL    hz,hzp,hzm
      REAL    hmin
      REAL    Dxm,Dxp,Dym,Dyp,Dzm,Dzp
      REAL    Dxc,Dyc,Dzc
      REAL    Dxl,Dxr
      REAL    Dyb,Dyt
      REAL    Dzd,Dzu
      REAL    rxm,wxm
      REAL    rxp,wxp
      REAL    rym,wym
      REAL    ryp,wyp
      REAL    rzm,wzm
      REAL    rzp,wzp
      REAL    h1,h2
      REAL    Dxx0,Dyy0,Dzz0
      REAL    sgn,eps,D,diff

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      eps = 1.d-10

      hmin = dmin1(hx,hy,hz)
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            sgn = S_eps(V(i0,i1,i2),hmin)

c           Compute all the required finite differences
            Dxc  = (V(i0+1,i1,i2) - V(i0-1,i1,i2))/(two*hx) 
            Dyc  = (V(i0,i1+1,i2) - V(i0,i1-1,i2))/(two*hy)
            Dzc  = (V(i0,i1,i2+1) - V(i0,i1,i2-1))/(two*hz)
            Dxl  = (three*V(i0,i1,i2) - four*V(i0-1,i1,i2) 
     &              + V(i0-2,i1,i2))/(two*hx)
            Dxr = (-three*V(i0,i1,i2) + four*V(i0+1,i1,i2) 
     &              - V(i0+2,i1,i2))/(two*hx)
            Dyb  = (three*V(i0,i1,i2) - four*V(i0,i1-1,i2) 
     &              + V(i0,i1-2,i2))/(two*hy)
            Dyt = (-three*V(i0,i1,i2) + four*V(i0,i1+1,i2) 
     &              - V(i0,i1+2,i2))/(two*hy)
            Dzd  = (three*V(i0,i1,i2) - four*V(i0,i1,i2-1)
     &              + V(i0,i1,i2-2))/(two*hz)
            Dzu = (-three*V(i0,i1,i2) + four*V(i0,i1,i2+1)
     &              - V(i0,i1,i2+2))/(two*hz)

            rxm = (eps + (V(i0,i1,i2) -two*V(i0-1,i1,i2) + 
     &             V(i0-2,i1,i2))**two)/(eps + (V(i0+1,i1,i2) -
     &             two*V(i0,i1,i2) + V(i0-1,i1,i2))**two)
            wxm = one/(one + two*rxm**two) 
            rxp = (eps + (V(i0+2,i1,i2) -two*V(i0+1,i1,i2) + 
     &             V(i0,i1,i2))**two)/(eps + (V(i0+1,i1,i2) - 
     &             two*V(i0,i1,i2) + V(i0-1,i1,i2))**two)
            wxp = one/(one + two*rxp**two)

            rym = (eps + (V(i0,i1,i2) -two*V(i0,i1-1,i2) + 
     &             V(i0,i1-2,i2))**two)/(eps + (V(i0,i1+1,i2) - 
     &             two*V(i0,i1,i2) + V(i0,i1-1,i2))**two)
            wym = one/(one + two*rym**two)
            ryp = (eps + (V(i0,i1+2,i2) -two*V(i0,i1+1,i2) +
     &             V(i0,i1,i2))**two)/(eps + (V(i0,i1+1,i2) - 
     &             two*V(i0,i1,i2) + V(i0,i1-1,i2))**two)
            wyp = one/(one + two*ryp**two)

            rzm = (eps + (V(i0,i1,i2) -two*V(i0,i1,i2-1) +
     &             V(i0,i1,i2-2))**two)/(eps + (V(i0,i1,i2+1) -
     &             two*V(i0,i1,i2) + V(i0,i1,i2-1))**two)
            wzm = one/(one + two*rzm**two)
            rzp = (eps + (V(i0,i1,i2+2) -two*V(i0,i1,i2+1) + 
     &             V(i0,i1,i2))**two)/(eps + (V(i0,i1,i2+1) -
     &             two*V(i0,i1,i2) + V(i0,i1,i2-1))**two)
            wzp = one/(one + two*rzp**two)
 

c           Set dummy values for hxp,hxm,hyp,hym
            hxp = 1.d12;hxm = 1.d12
            hyp = 1.d12;hym = 1.d12
            hzp = 1.d12;hzm = 1.d12

c           Compute ENO differences with subcell fix
            if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0+1,i1,i2)
     &          .lt. zero) then
              Dxx0 = minmod(V(i0-1,i1,i2)-two*V(i0,i1,i2)+V(i0+1,i1,i2),
     &                      V(i0,i1,i2)-two*V(i0+1,i1,i2)+V(i0+2,i1,i2))
              diff = V(i0,i1,i2)-V(i0+1,i1,i2)
              if (abs(Dxx0) .gt. eps) then
                D = (Dxx0/two-V(i0,i1,i2)-V(i0+1,i1,i2))**two
     &              -four*V(i0,i1,i2)*V(i0+1,i1,i2)
                hxp = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
              else
                hxp = hx*V(i0,i1,i2)/diff
              endif
        
              h1 = hx; h2 = hxp
              Dxc = (-V(i0,i1,i2)*h1**2 + (-V(i0-1,i1,i2) + V(i0,i1,i2))
     &               *h2**2)/(h1*h2*(h1 + h2))

              h1 = hxp; h2 = hx - hxp
              Dxr = (-V(i0+1,i1,i2)*h1**2 + two*(-V(i0,i1,i2))*h1*h2 + 
     &                (-V(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2)) 
            endif

            if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0-1,i1,i2)
     &         .lt. zero) then
              Dxx0 = minmod(V(i0-1,i1,i2)-two*V(i0,i1,i2)+V(i0+1,i1,i2),
     &                      V(i0,i1,i2)-two*V(i0-1,i1,i2)+V(i0-2,i1,i2))
              diff = V(i0,i1,i2)-V(i0-1,i1,i2)
              if (abs(Dxx0) .gt. eps) then
                D = (Dxx0/two-V(i0,i1,i2)-V(i0-1,i1,i2))**two
     &              -four*V(i0,i1,i2)*V(i0-1,i1,i2)
                hxm = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
              else
                hxm = hx*V(i0,i1,i2)/diff
              endif
       
              h1 = hxm; h2 = hx
              Dxc = ((-V(i0,i1,i2) + V(i0+1, i1,i2))*h1**2 + 
     &              (V(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

              h1 = hx - hxm; h2 = hxm
              Dxl = (V(i0,i1,i2)*h1**2 + two*V(i0,i1,i2)*h1*h2 + 
     &              V(i0-1,i1,i2)*h2**2)/(h1*h2*(h1 + h2))
            endif

            if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1+1,i2) 
     &         .lt. zero) then
              Dyy0 = minmod(V(i0,i1-1,i2)-two*V(i0,i1,i2)+V(i0,i1+1,i2),
     &                      V(i0,i1,i2)-two*V(i0,i1+1,i2)+V(i0,i1+2,i2))
              diff = V(i0,i1,i2)-V(i0,i1+1,i2)
              if (abs(Dyy0) .gt. eps) then
                D = (Dyy0/two-V(i0,i1,i2)-V(i0,i1+1,i2))**two
     &              -four*V(i0,i1,i2)*V(i0,i1+1,i2)
                hyp = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
              else
                hyp = hy*V(i0,i1,i2)/diff
              endif
 
              h1 = hy; h2 = hyp
              Dyc = (-V(i0,i1,i2)*h1**2 + (-V(i0,i1-1,i2) + V(i0,i1,i2))
     &              *h2**2)/(h1*h2*(h1 + h2))

              h1 = hyp; h2 = hy - hyp
              Dyt = (-V(i0,i1+1,i2)*h1**2 + two*(-V(i0,i1,i2))*h1*h2 +
     &                (-V(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))
            endif

            if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1-1,i2)
     &         .lt. zero) then
              Dyy0 = minmod(V(i0,i1-1,i2)-two*V(i0,i1,i2)+V(i0,i1+1,i2),
     &                      V(i0,i1,i2)-two*V(i0,i1-1,i2)+V(i0,i1-2,i2))
              diff = V(i0,i1,i2)-V(i0,i1-1,i2)
              if (abs(Dyy0) .gt. eps) then
                D = (Dyy0/two-V(i0,i1,i2)-V(i0,i1-1,i2))**two
     &              -four*V(i0,i1,i2)*V(i0,i1-1,i2)
                hym = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
              else
                hym = hy*V(i0,i1,i2)/diff
              endif

              h1 = hym; h2 = hy
              Dyc = ((-V(i0,i1,i2) + V(i0, i1+1,i2))*h1**2 +
     &              (V(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

              h1 = hy - hym; h2 = hym
              Dyb = (V(i0,i1,i2)*h1**2 + two*V(i0,i1,i2)*h1*h2 +
     &              V(i0,i1-1,i2)*h2**2)/(h1*h2*(h1 + h2))  
            endif

            if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1,i2+1)
     &         .lt. zero) then
              Dzz0 = minmod(V(i0,i1,i2-1)-two*V(i0,i1,i2)+V(i0,i1,i2+1),
     &                      V(i0,i1,i2)-two*V(i0,i1,i2+1)+V(i0,i1,i2+2))
              diff = V(i0,i1,i2)-V(i0,i1,i2+1)
              if (abs(Dzz0) .gt. eps) then
                D = (Dzz0/two-V(i0,i1,i2)-V(i0,i1,i2+1))**two
     &              -four*V(i0,i1,i2)*V(i0,i1,i2+1)
                hzp = hz*(half + (diff-sign(one,diff)*sqrt(D))/Dzz0)
              else
                hzp = hz*V(i0,i1,i2)/diff
              endif
 
              h1 = hz; h2 = hzp
              Dzc = (-V(i0,i1,i2)*h1**2 +
     &              (-V(i0,i1,i2-1) + V(i0,i1,i2))*h2**2)/
     &               (h1*h2*(h1 + h2))

              h1 = hzp; h2 = hz - hzp
              Dzu = (-V(i0,i1,i2+1)*h1**2 + two*(-V(i0,i1,i2))*h1*h2 +
     &                (-V(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))
            endif

            if (use_subcell .ne. 0 .and. V(i0,i1,i2)*V(i0,i1,i2-1)
     &         .lt. zero) then
              Dzz0 = minmod(V(i0,i1,i2-1)-two*V(i0,i1,i2)+V(i0,i1,i2+1),
     &                      V(i0,i1,i2)-two*V(i0,i1,i2-1)+V(i0,i1,i2-2))
              diff = V(i0,i1,i2)-V(i0,i1,i2-1)
              if (abs(Dzz0) .gt. eps) then
                D = (Dzz0/two-V(i0,i1,i2)-V(i0,i1,i2-1))**two
     &              -four*V(i0,i1,i2)*V(i0,i1,i2-1)
                hzm = hz*(half + (diff-sign(one,diff)*sqrt(D))/Dzz0)
              else
                hzm = hz*V(i0,i1,i2)/diff
              endif

              h1 = hzm; h2 = hz
              Dzc = ((-V(i0,i1,i2) + V(i0, i1,i2+1))*h1**2 +
     &              (V(i0,i1,i2))*h2**2)/(h1*h2*(h1 + h2))

              h1 = hz - hzm; h2 = hzm
              Dzd = (V(i0,i1,i2)*h1**2 + two*V(i0,i1,i2)*h1*h2 +
     &              V(i0,i1,i2-1)*h2**2)/(h1*h2*(h1 + h2))  
            endif

c           Compute first order derivatives

            Dxm = (one - wxm)*Dxc + wxm*Dxl
            Dxp = (one - wxp)*Dxc + wxp*Dxr

            Dym = (one - wym)*Dyc + wym*Dyb
            Dyp = (one - wyp)*Dyc + wyp*Dyt

            Dzm = (one - wzm)*Dzc + wzm*Dzd
            Dzp = (one - wzp)*Dzc + wzp*Dzu

            H(i0,i1,i2) = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)
          enddo
        enddo
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out fifth order relaxation scheme using Gauss Seidel updates
c     
c     Uses fifth order WENO for spatial discretization
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxationls5thorderweno3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     dir,
     &     use_sign_fix)
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
      INTEGER use_sign_fix

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
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i2 = ilower2,iupper2
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 3) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 4) then
        do i2 = ilower2,iupper2
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 5) then
        do i2 = iupper2,ilower2,-1
          do i1 = ilower1,iupper1
            do i0 = iupper0,ilower0,-1
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 6) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = ilower0,iupper0
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      elseif (dir .eq. 7) then
        do i2 = iupper2,ilower2,-1
          do i1 = iupper1,ilower1,-1
            do i0 = iupper0,ilower0,-1
                call evalrelax5thorderweno3d(U,U_gcw,V,V_gcw,
     &                                       ilower0,iupper0,
     &                                       ilower1,iupper1,
     &                                       ilower2,iupper2,
     &                                       i0,i1,i2,dx,
     &                                       use_sign_fix)
            enddo
          enddo
        enddo
      endif
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out single fifth order sweep using a WENO stencil
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine evalrelax5thorderweno3d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     i0,i1,i2,
     &     dx,
     &     use_sign_fix)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL HG, WENO5, S_eps

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
      INTEGER i0,i1,i2,k,np,nm
      REAL    hx,hy,hz,hmin
      REAL    Dxm,Dxp,Dym,Dyp,Dzm,Dzp
      REAL    H,dt,sgn,cfl
      REAL    Qx(-2:1),Qy(-2:1),Qz(-2:1)
      REAL    Qxxp(-2:1),Qxxm(-2:1)
      REAL    Qyyp(-2:1),Qyym(-2:1)
      REAL    Qzzp(-2:1),Qzzm(-2:1)
      REAL    Ex,Ey,Ez
      INTEGER use_sign_fix

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      cfl = 0.3d0
      hmin = dmin1(hx,hy,hz)
      sgn = S_eps(V(i0,i1,i2),hmin)

c     Sign fix
      if (use_sign_fix .ne. 0) then
        if (V(i0,i1,i2)*V(i0+1,i1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0+1,i1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0-1,i1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0-1,i1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1+1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1+1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1-1,i2) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1-1,i2))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1,i2+1) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1,i2+1))) then
          sgn = zero
        endif
        if (V(i0,i1,i2)*V(i0,i1,i2-1) .lt. zero .and.
     &      abs(V(i0,i1,i2)) .le. abs(V(i0,i1,i2-1))) then
          sgn = zero
        endif
      endif

c     Compute all the required finite differences and their WENO5 interpolation
      do k = -2,1
        np = 1-k
        nm = k+1
        Qx(k) = (U(i0+k+1,i1,i2) - U(i0+k,i1,i2))/hx
        Qy(k) = (U(i0,i1+k+1,i2) - U(i0,i1+k,i2))/hy
        Qz(k) = (U(i0,i1,i2+k+1) - U(i0,i1,i2+k))/hy
        Qxxp(k) = (U(i0+np,i1,i2)-two*U(i0+np-1,i1,i2)+U(i0+np-2,i1,i2))
     &             /hx
        Qyyp(k) = (U(i0,i1+np,i2)-two*U(i0,i1+np-1,i2)+U(i0,i1+np-2,i2))
     &             /hy
        Qzzp(k) = (U(i0,i1,i2+np)-two*U(i0,i1,i2+np-1)+U(i0,i1,i2+np-2))
     &             /hz
        Qxxm(k) = (U(i0+nm,i1,i2)-two*U(i0+nm-1,i1,i2)+U(i0+nm-2,i1,i2))
     &             /hx
        Qyym(k) = (U(i0,i1+nm,i2)-two*U(i0,i1+nm-1,i2)+U(i0,i1+nm-2,i2))
     &             /hy
        Qzzm(k) = (U(i0,i1,i2+nm)-two*U(i0,i1,i2+nm-1)+U(i0,i1,i2+nm-2))
     &             /hz
      enddo

      Ex = (-Qx(-2)+7.d0*(Qx(-1)+Qx(0))-Qx(1))/12.d0
      Ey = (-Qy(-2)+7.d0*(Qy(-1)+Qy(0))-Qy(1))/12.d0
      Ez = (-Qz(-2)+7.d0*(Qz(-1)+Qz(0))-Qz(1))/12.d0
      Dxp = Ex+WENO5(Qxxp)
      Dxm = Ex-WENO5(Qxxm)
      Dyp = Ey+WENO5(Qyyp)
      Dym = Ey-WENO5(Qyym)
      Dzp = Ez+WENO5(Qzzp)
      Dzm = Ez-WENO5(Qzzm)
    
      H = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)
      dt = cfl*hmin

      U(i0,i1,i2) = U(i0,i1,i2) - dt*sgn*(H-one)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Godunov Hamiltonian of the indicator field |grad phi_0|
c     
c     Uses fifth order WENO for spatial discretization
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine godunovhamiltonian5thorderweno3d(
     &     H,H_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL HG, WENO5, S_eps

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER H_gcw,V_gcw

c
c     Input/Output.
c
      REAL H(CELL3d(ilower,iupper,H_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2,k,np,nm
      REAL    hx,hy,hz,hmin
      REAL    Dxm,Dxp,Dym,Dyp,Dzm,Dzp
      REAL    Qx(-2:1),Qy(-2:1),Qz(-2:1)
      REAL    Qxxp(-2:1),Qxxm(-2:1)
      REAL    Qyyp(-2:1),Qyym(-2:1)
      REAL    Qzzp(-2:1),Qzzm(-2:1)
      REAL    Ex,Ey,Ez
      REAL    sgn

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      hmin = dmin1(hx,hy,hz)

      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
            sgn = S_eps(V(i0,i1,i2),hmin)

c           Compute all the required finite differences and their WENO5 interpolation
            do k = -2,1
              np = 1-k
              nm = k+1
              Qx(k) = (V(i0+k+1,i1,i2) - V(i0+k,i1,i2))/hx
              Qy(k) = (V(i0,i1+k+1,i2) - V(i0,i1+k,i2))/hy
              Qz(k) = (V(i0,i1,i2+k+1) - V(i0,i1,i2+k))/hz
              Qxxp(k) = (V(i0+np,i1,i2)-two*V(i0+np-1,i1,i2)
     &                  +V(i0+np-2,i1,i2))/hx
              Qyyp(k) = (V(i0,i1+np,i2)-two*V(i0,i1+np-1,i2)
     &                  +V(i0,i1+np-2,i2))/hy
              Qzzp(k) = (V(i0,i1,i2+np)-two*V(i0,i1,i2+np-1)
     &                  +V(i0,i1,i2+np-2))/hz
              Qxxm(k) = (V(i0+nm,i1,i2)-two*V(i0+nm-1,i1,i2)
     &                  +V(i0+nm-2,i1,i2))/hx
              Qyym(k) = (V(i0,i1+nm,i2)-two*V(i0,i1+nm-1,i2)
     &                  +V(i0,i1+nm-2,i2))/hy
              Qzzm(k) = (V(i0,i1,i2+nm)-two*V(i0,i1,i2+nm-1))
     &                  +V(i0,i1,i2+nm-2)/hz
            enddo

            Ex = (-Qx(-2)+7.d0*(Qx(-1)+Qx(0))-Qx(1))/12.d0
            Ey = (-Qy(-2)+7.d0*(Qy(-1)+Qy(0))-Qy(1))/12.d0
            Ez = (-Qz(-2)+7.d0*(Qz(-1)+Qz(0))-Qz(1))/12.d0
            Dxp = Ex+WENO5(Qxxp)
            Dxm = Ex-WENO5(Qxxm)
            Dyp = Ey+WENO5(Qyyp)
            Dym = Ey-WENO5(Qyym)
            Dzp = Ez+WENO5(Qzzp)
            Dzm = Ez-WENO5(Qzzm)
    
            H(i0,i1,i2) = HG(Dxp,Dxm,Dyp,Dym,Dzp,Dzm,sgn)
          enddo
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Mass constraint on level set to ensure that it does not lose volume
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine projectlsmassconstraint3d(
     &     U,U_gcw,
     &     C,C_gcw,
     &     V,V_gcw,
     &     H,H_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Functions.
c
      REAL D_eps

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,C_gcw,V_gcw,H_gcw

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL C(CELL3d(ilower,iupper,C_gcw))
      REAL V(CELL3d(ilower,iupper,V_gcw))
      REAL H(CELL3d(ilower,iupper,H_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    hx,hy,hz,hmin
      REAL    lambda
      REAL    dijk,dlmn,phi0
      INTEGER l,m,n
      REAL    nmr,dnr
      REAL    w(-1:1,-1:1,-1:1)
      LOGICAL near_interface

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      hmin = dmin1(hx,hy,hz)

c     Compute integration weights based on Simpson's rule
      do l = -1,1
        do m = -1,1
          do n = -1,1
            if (abs(l) + abs(m) + abs(n) .eq. 1) then
              w(l,m,n) = 100.d0
            else if (abs(l) + abs(m) + abs(n) .eq. 2) then
              w(l,m,n) = 10.d0
            else if (abs(l) + abs(m) + abs(n) .eq. 3) then
              w(l,m,n) = 1.d0
            else
              w(l,m,n) = 0.d0
            endif
          enddo
        enddo
      enddo

      do i2 = ilower2,iupper2
       do i1 = ilower1,iupper1
        do i0 = ilower0,iupper0
c           If the point to be updated is not near the interface, then do not attempt to update it,
c           as this can cause the level set variable to blow up
            near_interface = (V(i0,i1,i2)*V(i0+1,i1,i2) .le. zero .or.
     &                        V(i0,i1,i2)*V(i0,i1+1,i2) .le. zero .or.
     &                        V(i0,i1,i2)*V(i0,i1,i2+1) .le. zero .or.
     &                        V(i0,i1,i2)*V(i0-1,i1,i2) .le. zero .or.
     &                        V(i0,i1,i2)*V(i0,i1-1,i2) .le. zero .or.
     &                        V(i0,i1,i2)*V(i0,i1,i2-1) .le. zero)
            if (.not. near_interface) then
              cycle
            endif
            phi0 = V(i0,i1,i2)
            dijk = D_eps(phi0, hmin)
            nmr = 1000.d0*dijk*(C(i0,i1,i2) - phi0)
            dnr = 1000.d0*(dijk**2)*H(i0,i1,i2)
            do l = -1,1
              do m = -1,1
                do n = -1,1
                  phi0 = V(i0+l,i1+m,i2+n)
                  dlmn = D_eps(phi0, hmin)
                  nmr = nmr + w(l,m,n)*dlmn*(C(i0+l,i1+m,i2+n) - phi0)
                  dnr = dnr + w(l,m,n)*(dlmn**2)*H(i0+l,i1+m,i2+n)
                enddo
              enddo
            enddo
            nmr = hx*hy*hz/1728.d0*nmr
            dnr = hx*hy*hz/1728.d0*dnr
            lambda = -nmr/dnr

            if (dnr .gt. zero) then
              U(i0,i1,i2) = C(i0,i1,i2) + lambda*dijk*H(i0,i1,i2)
            endif

        enddo
       enddo
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Volume constraint on level set to ensure that it does not lose volume
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine applylsvolumeshift3d(
     &     U,U_gcw,
     &     C,C_gcw,
     &     dV,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER U_gcw,C_gcw

c
c     Input/Output.
c
      REAL U(CELL3d(ilower,iupper,U_gcw))
      REAL C(CELL3d(ilower,iupper,C_gcw))
      REAL dV
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1,i2
      REAL    hx,hy,hz,hmin
      REAL    dt

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)
      hmin = dmin1(hx,hy,hz)
      dt = one

      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
                U(i0,i1,i2) = C(i0,i1,i2) + dt * dV
          enddo
        enddo
      enddo
      return
      end

 
