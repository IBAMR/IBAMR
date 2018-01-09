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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Godunov Hamiltonian.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      REAL function HG(a,b,c,d,sgn)
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
      REAL a,b,c,d,sgn
      REAL am,ap,bm,bp,cm,cp,dm,dp
      if (sgn .ge. zero) then
        am = dmin1(a,zero)
        bp = dmax1(b,zero)
        cm = dmin1(c,zero)
        dp = dmax1(d,zero)
        HG = sqrt(dmax1(am**two,bp**two) + dmax1(cm**two,dp**two))
      else
        ap = dmax1(a,zero)
        bm = dmin1(b,zero)
        cp = dmax1(c,zero)
        dm = dmin1(d,zero)
        HG = sqrt(dmax1(ap**two,bm**two) + dmax1(cp**two,dm**two))
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out first-order accurate fast sweeping algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fastsweep1storder2d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dlower0,dupper0,
     &     dlower1,dupper1,
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
      INTEGER dlower0,dupper0
      INTEGER dlower1,dupper1
      INTEGER U_gcw
      INTEGER patch_touches_bdry

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
      INTEGER touches_wall_loc_idx(0:2*NDIM - 1)
c
c     Local variables.
c
      INTEGER i0,i1
    
      
c     Do the four sweeping directions.
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            call evalsweep1storder2d(U,U_gcw,
     &                               ilower0,iupper0,
     &                               ilower1,iupper1,
     &                               i0,i1,
     &                               dlower0,dupper0,
     &                               dlower1,dupper1,
     &                               dx,
     &                               patch_touches_bdry,
     &                               touches_wall_loc_idx)
         enddo
      enddo
      
      do i1 = ilower1,iupper1
         do i0 = iupper0,ilower0,-1
            call evalsweep1storder2d(U,U_gcw,
     &                               ilower0,iupper0,
     &                               ilower1,iupper1,
     &                               i0,i1,
     &                               dlower0,dupper0,
     &                               dlower1,dupper1,
     &                               dx,
     &                               patch_touches_bdry,
     &                               touches_wall_loc_idx)     
         enddo
      enddo
      
      do i1 = iupper1,ilower1,-1
         do i0 = iupper0,ilower0,-1
            call evalsweep1storder2d(U,U_gcw,
     &                               ilower0,iupper0,
     &                               ilower1,iupper1,
     &                               i0,i1,
     &                               dlower0,dupper0,
     &                               dlower1,dupper1,
     &                               dx,
     &                               patch_touches_bdry,
     &                               touches_wall_loc_idx)
         enddo
      enddo
      
      do i1 = iupper1,ilower1,-1
         do i0 = ilower0,iupper0
            call evalsweep1storder2d(U,U_gcw,
     &                               ilower0,iupper0,
     &                               ilower1,iupper1,
     &                               i0,i1,
     &                               dlower0,dupper0,
     &                               dlower1,dupper1,
     &                               dx,
     &                               patch_touches_bdry,
     &                               touches_wall_loc_idx)
         enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute fast sweep solution at a given grid cell
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine evalsweep1storder2d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     i0,i1,
     &     dlower0,dupper0,
     &     dlower1,dupper1,
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
      INTEGER dlower0,dupper0
      INTEGER dlower1,dupper1
      INTEGER U_gcw
      INTEGER patch_touches_bdry
     
c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
      INTEGER touches_wall_loc_idx(0:2*NDIM - 1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    a,b,sgn
      REAL    hx,hy
      REAL    Q,R,S
      REAL    dbar
      REAL    U_wall_coef 
      REAL    h_wall_coef
      
      if (patch_touches_bdry .eq. 1) then  
        if ((i0 .eq. dlower0) .and. 
     &     touches_wall_loc_idx(0) .eq. 1) then
           U_wall_coef   = zero
           h_wall_coef   = half
        elseif ((i0 .eq. dupper0) .and. 
     &     touches_wall_loc_idx(1) .eq. 1) then
           U_wall_coef   = zero
           h_wall_coef   = half
        elseif ((i1 .eq. dlower1) .and. 
     &     touches_wall_loc_idx(2) .eq. 1) then
           U_wall_coef   = zero
           h_wall_coef   = half
        elseif ((i1 .eq. dupper1) .and. 
     &     touches_wall_loc_idx(3) .eq. 1) then
           U_wall_coef   = zero
           h_wall_coef   = half
        else
          U_wall_coef   = one
          h_wall_coef   = one
        endif
      endif

c     Carry out a single sweep
      if (U(i0,i1) .eq. zero) then
        sgn = zero
      else
        sgn = sign(one,U(i0,i1))
      endif

      hx = dx(0)
      hy = dx(1)
      a  = sgn*dmin1(sgn*U(i0-1,i1),sgn*U(i0+1,i1))
      b  = sgn*dmin1(sgn*U(i0,i1-1),sgn*U(i0,i1+1))
 
c     Take care of physical boundaries.
      if (patch_touches_bdry .eq. 1) then
         if (i0 .eq. dlower0) then 
            a  = U(i0+1,i1)*U_wall_coef   
            hx = hx*h_wall_coef                
         elseif (i0 .eq. dupper0) then
            a  = U(i0-1,i1)*U_wall_coef
            hx = hx*h_wall_coef
         elseif (i1 .eq. dlower1) then 
            b  = U(i0,i1+1)*U_wall_coef   
            hy = hy*h_wall_coef                
         elseif (i1 .eq. dupper1) then
            b  = U(i0,i1-1)*U_wall_coef
            hy = hy*h_wall_coef
         endif
      endif

      if (sgn*(b-a) .gt. hx) then
        dbar = a + sgn*hx
      elseif (sgn*(a-b) .gt. hy) then
        dbar = b + sgn*hy
      else
        Q = hx*hx + hy*hy
        R = -2.d0*(hy*hy*a + hx*hx*b)
        S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
        dbar = (-R + sgn*sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
      endif
                        
      U(i0,i1) = sgn*dmin1(sgn*U(i0,i1),sgn*dbar)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out first order relaxation scheme using Gauss Seidel updates
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxationls1storder2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
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
      INTEGER U_gcw,V_gcw
      INTEGER dir

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1

      if (dir .eq. 0) then
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
              call evalrelax1storder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i1 = ilower1,iupper1
          do i0 = iupper0,ilower0,-1
              call evalrelax1storder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i1 = iupper1,ilower1,-1
          do i0 = ilower0,iupper0
              call evalrelax1storder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
          enddo
        enddo
      elseif (dir .eq. 3 )then
        do i1 = iupper1,ilower1,-1
          do i0 = iupper0,ilower0,-1
              call evalrelax1storder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
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
      subroutine evalrelax1storder2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     i0,i1,
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
      INTEGER U_gcw,V_gcw

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    hx,hy
      REAL    sgn
      REAL    dt
      REAL    G
      REAL    Dxp,Dxm
      REAL    Dyp,Dym

      hx = dx(0)
      hy = dx(1)
      dt = half*dmin1(hx,hy)
      if (V(i0,i1) .eq. zero) then
         sgn = zero
      else
         sgn = sign(one,V(i0,i1))
      endif

      Dxm = one/hx*(U(i0,i1)-U(i0-1,i1))
      Dxp = one/hx*(U(i0+1,i1)-U(i0,i1))
      Dym = one/hy*(U(i0,i1)-U(i0,i1-1))
      Dyp = one/hy*(U(i0,i1+1)-U(i0,i1))
      G = HG(Dxp,Dxm,Dyp,Dym,sgn)

      U(i0,i1) = U(i0,i1) - dt*sgn*(G-one)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out third order relaxation scheme using Gauss Seidel updates
c     NOTE: this scheme is between between third and fourth
c     order near the interface and second order everywhere else
c
c     Uses second order WENO for spatial discretization with a subcell
c     fix near the interface
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxationls3rdorder2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
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
      INTEGER U_gcw,V_gcw
      INTEGER dir

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1

      if (dir .eq. 0) then
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0
              call evalrelax3rdorder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
          enddo
        enddo
      elseif (dir .eq. 1) then
        do i1 = ilower1,iupper1
          do i0 = iupper0,ilower0,-1
              call evalrelax3rdorder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
          enddo
        enddo
      elseif (dir .eq. 2) then
        do i1 = iupper1,ilower1,-1
          do i0 = ilower0,iupper0
              call evalrelax3rdorder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
          enddo
        enddo
      elseif (dir .eq. 3 )then
        do i1 = iupper1,ilower1,-1
          do i0 = iupper0,ilower0,-1
              call evalrelax3rdorder2d(U,U_gcw,V,V_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,dx)
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
      subroutine evalrelax3rdorder2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     i0,i1,
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
      INTEGER U_gcw,V_gcw

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL V(CELL2d(ilower,iupper,V_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    hx,hxp,hxm
      REAL    hy,hyp,hym
      REAL    Dxm,Dxp,Dym,Dyp
      REAL    Dxx,Dxxp,Dxxm,Dyy,Dyyp,Dyym
      REAL    Dxx0,Dyy0
      REAL    H,dt,sgn,cfl,eps,D,diff

      hx = dx(0)
      hy = dx(1)
      cfl = 0.45d0
      eps = 1.d-10

      if (V(i0,i1) .eq. zero) then
         sgn = zero
      else
         sgn = sign(one,V(i0,i1))
      endif

c     Compute all the required finite differences
      Dxx  = (U(i0-1,i1) - two*U(i0,i1) + U(i0+1,i1))/(hx**two)
      Dxxp = (U(i0,i1) - two*U(i0+1,i1) + U(i0+2,i1))/(hx**two)
      Dxxm = (U(i0-2,i1) - two*U(i0-1,i1) + U(i0,i1))/(hx**two)
      Dyy  = (U(i0,i1-1) - two*U(i0,i1) + U(i0,i1+1))/(hy**two)
      Dyyp = (U(i0,i1) - two*U(i0,i1+1) + U(i0,i1+2))/(hy**two)
      Dyym = (U(i0,i1-2) - two*U(i0,i1-1) + U(i0,i1))/(hy**two)

c     Set dummy values for hxp,hxm,hyp,hym
      hxp = 1.d12;hxm = 1.d12;hyp = 1.d12;hym = 1.d12

c     Compute ENO differences with subcell fix
      if (V(i0,i1)*V(i0+1,i1) .lt. zero) then
        Dxx0 = minmod(V(i0-1,i1)-two*V(i0,i1)+V(i0+1,i1),
     &                V(i0,i1)-two*V(i0+1,i1)+V(i0+2,i1))
        diff = V(i0,i1)-V(i0+1,i1)
        if (abs(Dxx0) .gt. eps) then
          D = (Dxx0/two-V(i0,i1)-V(i0+1,i1))**two
     &        -four*V(i0,i1)*V(i0+1,i1)
          hxp = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
        else
          hxp = hx*V(i0,i1)/diff
        endif

        Dxp = (zero-U(i0,i1))/hxp - hxp/two*minmod(Dxx,Dxxp)
      else
        Dxp = (U(i0+1,i1)-U(i0,i1))/hx - hx/two*minmod(Dxx,Dxxp)
      endif

      if (V(i0,i1)*V(i0-1,i1) .lt. zero) then
        Dxx0 = minmod(V(i0-1,i1)-two*V(i0,i1)+V(i0+1,i1),
     &                V(i0,i1)-two*V(i0-1,i1)+V(i0-2,i1))
        diff = V(i0,i1)-V(i0-1,i1)
        if (abs(Dxx0) .gt. eps) then
          D = (Dxx0/two-V(i0,i1)-V(i0-1,i1))**two
     &        -four*V(i0,i1)*V(i0-1,i1)
          hxm = hx*(half + (diff-sign(one,diff)*sqrt(D))/Dxx0)
        else
          hxm = hx*V(i0,i1)/diff
        endif

        Dxm = (U(i0,i1)-zero)/hxm + hxm/two*minmod(Dxx,Dxxm)
      else
        Dxm = (U(i0,i1)-U(i0-1,i1))/hx + hx/two*minmod(Dxx,Dxxm)
      endif

      if (V(i0,i1)*V(i0,i1+1) .lt. zero) then
        Dyy0 = minmod(V(i0,i1-1)-two*V(i0,i1)+V(i0,i1+1),
     &                V(i0,i1)-two*V(i0,i1+1)+V(i0,i1+2))
        diff = V(i0,i1)-V(i0,i1+1)
        if (abs(Dyy0) .gt. eps) then
          D = (Dyy0/two-V(i0,i1)-V(i0,i1+1))**two
     &        -four*V(i0,i1)*V(i0,i1+1)
          hyp = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
        else
          hyp = hy*V(i0,i1)/diff
        endif

        Dyp = (zero-U(i0,i1))/hyp - hyp/two*minmod(Dyy,Dyyp)
      else
        Dyp = (U(i0,i1+1)-U(i0,i1))/hy - hy/two*minmod(Dyy,Dyyp)
      endif

      if (V(i0,i1)*V(i0,i1-1) .lt. zero) then
        Dyy0 = minmod(V(i0,i1-1)-two*V(i0,i1)+V(i0,i1+1),
     &                V(i0,i1)-two*V(i0,i1-1)+V(i0,i1-2))
        diff = V(i0,i1)-V(i0,i1-1)
        if (abs(Dyy0) .gt. eps) then
          D = (Dyy0/two-V(i0,i1)-V(i0,i1-1))**two
     &        -four*V(i0,i1)*V(i0,i1-1)
          hym = hy*(half + (diff-sign(one,diff)*sqrt(D))/Dyy0)
        else
          hym = hy*V(i0,i1)/diff
        endif

        Dym = (U(i0,i1)-zero)/hym + hym/two*minmod(Dyy,Dyym)
      else
        Dym = (U(i0,i1)-U(i0,i1-1))/hy + hy/two*minmod(Dyy,Dyym)
      endif

      H = HG(Dxp,Dxm,Dyp,Dym,sgn)
      dt = cfl*dmin1(hx,hy,hxp,hxm,hyp,hym)

      if (dt .gt. zero) then
        U(i0,i1) = U(i0,i1) - dt*sgn*(H-one)
      endif

      return
      end
