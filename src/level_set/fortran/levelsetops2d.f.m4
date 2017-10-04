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
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
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
     &     consider_bdry_wall)
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
      INTEGER consider_bdry_wall

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    U_wall_coef
      REAL    h_wall_coef

      if (consider_bdry_wall .eq. 1) then
        U_wall_coef   = zero
        h_wall_coef   = half
      else
        U_wall_coef   = one
        h_wall_coef   = one
      endif
     
      
c     Do the four sweeping directions.
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            call singlesweep1storder2d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
         enddo
      enddo
      
      do i1 = ilower1,iupper1
         do i0 = iupper0,ilower0,-1
            call singlesweep1storder2d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
         enddo
      enddo
      
      do i1 = iupper1,ilower1,-1
         do i0 = iupper0,ilower0,-1
            call singlesweep1storder2d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
         enddo
      enddo
      
      do i1 = iupper1,ilower1,-1
         do i0 = ilower0,iupper0
            call singlesweep1storder2d(U,U_gcw,
     &                                 ilower0,iupper0,
     &                                 ilower1,iupper1,
     &                                 i0,i1,
     &                                 dlower0,dupper0,
     &                                 dlower1,dupper1,
     &                                 dx,
     &                                 patch_touches_bdry,
     &                                 U_wall_coef,h_wall_coef)
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
      subroutine singlesweep1storder2d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     i0,i1,
     &     dlower0,dupper0,
     &     dlower1,dupper1,
     &     dx,
     &     patch_touches_bdry,
     &     U_wall_coef,h_wall_coef)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER dlower0,dupper0
      INTEGER dlower1,dupper1
      INTEGER U_gcw
      INTEGER patch_touches_bdry
      REAL    U_wall_coef,h_wall_coef

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    a,b,sgn
      REAL    hx,hy
      REAL    Q,R,S
      REAL    dbar
      
c     Carry out a single sweep
      if (U(i0,i1) .eq. 0.d0) then
        sgn = 0.d0
      else
        sgn = sign(1.d0,U(i0,i1))
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

      if (sgn .gt. 0.d0) then 
        if (b-a .gt. hx) then
           dbar = a + hx
        elseif (a-b .gt. hy) then
           dbar = b + hy
        else
           Q = hx*hx + hy*hy
           R = -2.d0*(hy*hy*a + hx*hx*b)
           S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
           dbar = (-R + sgn*sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
        endif
      elseif (sgn .lt. 0.d0) then 
        if (b-a .lt. -hx) then
           dbar = a - hx
        elseif (a-b .lt. -hy) then
           dbar = b - hy
        else
           Q = hx*hx + hy*hy
           R = -2.d0*(hy*hy*a + hx*hx*b)
           S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
           dbar = (-R + sgn*sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
        endif
      else
        dbar = 0.d0
      endif
                        
      U(i0,i1) = sgn*dmin1(sgn*U(i0,i1),sgn*dbar)

      return
      end


c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Carry out second-order accurate fast sweeping algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fastsweep2ndorder2d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dlower0,dupper0,
     &     dlower1,dupper1,
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
      INTEGER dlower0,dupper0
      INTEGER dlower1,dupper1
      INTEGER U_gcw
      INTEGER patch_touches_bdry
      INTEGER consider_bdry_wall

c
c     Input/Output.
c
      REAL U(CELL2d(ilower,iupper,U_gcw))
      REAL dx(0:NDIM-1)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    a,b
      REAL    aplus,bplus
      REAL    aminus,bminus
      REAL    hx,hy
      REAL    Q,R,S
      REAL    dbar
      REAL    U_wall_coef
      REAL    h_wall_coef


      if (consider_bdry_wall .eq. 1) then
        U_wall_coef   = zero
        h_wall_coef   = threefourth
      else
        U_wall_coef   = one
        h_wall_coef   = one
      endif
     
      
c     Do the four sweeping directions.
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0

            hx = twothird*dx(0)
            hy = twothird*dx(1)
            aminus = fourthird*U(i0-1,i1) - third*U(i0-2,i1)
            aplus  = fourthird*U(i0+1,i1) - third*U(i0+2,i1)
            bminus = fourthird*U(i0,i1-1) - third*U(i0,i1-2)
            bplus  = fourthird*U(i0,i1+1) - third*U(i0,i1+2)
            
            if (U(i0-1,i1) .le. U(i0+1,i1)) then
              a  = aminus
            else
              a  = aplus
            endif

            if (U(i0,i1-1) .le. U(i0,i1+1)) then
              b  = bminus
            else
              b  = bplus
            endif
       
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then 
                a  = aplus*U_wall_coef   
                hx = hx*h_wall_coef
              elseif (i0 .eq. dlower0+1) then
                a  = aplus
                hx = hx                
              elseif (i0 .eq. dupper0) then
                a  = aminus*U_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0-1) then
                a  = aminus
                hx = hx
              elseif (i1 .eq. dlower1) then 
                b  = bplus*U_wall_coef   
                hy = hy*h_wall_coef
              elseif (i1 .eq. dlower1+1) then 
                b  = bplus
                hy = hy                 
              elseif (i1 .eq. dupper1) then
                b  = bminus*U_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1-1) then
                b  = bminus
                hy = hy
              endif
            endif

            if (b-a .gt. hx) then
              dbar = a + hx
            elseif (a-b .gt. hy) then
              dbar = b + hy
            else
              Q = hx*hx + hy*hy
              R = -2.d0*(hy*hy*a + hx*hx*b)
              S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
              dbar = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
            endif
                        
            U(i0,i1) = dmin1(U(i0,i1),dbar)

            enddo
      enddo
      
      do i1 = ilower1,iupper1
         do i0 = iupper0,ilower0,-1

            hx = twothird*dx(0)
            hy = twothird*dx(1)
            aminus = fourthird*U(i0-1,i1) - third*U(i0-2,i1)
            aplus  = fourthird*U(i0+1,i1) - third*U(i0+2,i1)
            bminus = fourthird*U(i0,i1-1) - third*U(i0,i1-2)
            bplus  = fourthird*U(i0,i1+1) - third*U(i0,i1+2)
            
            if (U(i0-1,i1) .le. U(i0+1,i1)) then
              a  = aminus
            else
              a  = aplus
            endif

            if (U(i0,i1-1) .le. U(i0,i1+1)) then
              b  = bminus
            else
              b  = bplus
            endif
       
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then 
                a  = aplus*U_wall_coef   
                hx = hx*h_wall_coef
              elseif (i0 .eq. dlower0+1) then
                a  = aplus
                hx = hx                
              elseif (i0 .eq. dupper0) then
                a  = aminus*U_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0-1) then
                a  = aminus
                hx = hx
              elseif (i1 .eq. dlower1) then 
                b  = bplus*U_wall_coef   
                hy = hy*h_wall_coef
              elseif (i1 .eq. dlower1+1) then 
                b  = bplus
                hy = hy                 
              elseif (i1 .eq. dupper1) then
                b  = bminus*U_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1-1) then
                b  = bminus
                hy = hy
              endif
            endif

            
            if (b-a .gt. hx) then
              dbar = a + hx
            elseif (a-b .gt. hy) then
              dbar = b + hy
            else
              Q = hx*hx + hy*hy
              R = -2.d0*(hy*hy*a + hx*hx*b)
              S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
              dbar = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
            endif
            
            U(i0,i1) = dmin1(U(i0,i1),dbar)

         enddo
      enddo
      
      do i1 = iupper1,ilower1,-1
         do i0 = iupper0,ilower0,-1

            hx = twothird*dx(0)
            hy = twothird*dx(1)
            aminus = fourthird*U(i0-1,i1) - third*U(i0-2,i1)
            aplus  = fourthird*U(i0+1,i1) - third*U(i0+2,i1)
            bminus = fourthird*U(i0,i1-1) - third*U(i0,i1-2)
            bplus  = fourthird*U(i0,i1+1) - third*U(i0,i1+2)
            
            if (U(i0-1,i1) .le. U(i0+1,i1)) then
              a  = aminus
            else
              a  = aplus
            endif

            if (U(i0,i1-1) .le. U(i0,i1+1)) then
              b  = bminus
            else
              b  = bplus
            endif
       
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then 
                a  = aplus*U_wall_coef   
                hx = hx*h_wall_coef
              elseif (i0 .eq. dlower0+1) then
                a  = aplus
                hx = hx                
              elseif (i0 .eq. dupper0) then
                a  = aminus*U_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0-1) then
                a  = aminus
                hx = hx
              elseif (i1 .eq. dlower1) then 
                b  = bplus*U_wall_coef   
                hy = hy*h_wall_coef
              elseif (i1 .eq. dlower1+1) then 
                b  = bplus
                hy = hy                 
              elseif (i1 .eq. dupper1) then
                b  = bminus*U_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1-1) then
                b  = bminus
                hy = hy
              endif
            endif

            
            if (b-a .gt. hx) then
              dbar = a + hx
            elseif (a-b .gt. hy) then
              dbar = b + hy
            else
              Q = hx*hx + hy*hy
              R = -2.d0*(hy*hy*a + hx*hx*b)
              S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
              dbar = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
            endif
            
            U(i0,i1) = dmin1(U(i0,i1),dbar)

         enddo
      enddo
      
      do i1 = iupper1,ilower1,-1
         do i0 = ilower0,iupper0

            hx = twothird*dx(0)
            hy = twothird*dx(1)
            aminus = fourthird*U(i0-1,i1) - third*U(i0-2,i1)
            aplus  = fourthird*U(i0+1,i1) - third*U(i0+2,i1)
            bminus = fourthird*U(i0,i1-1) - third*U(i0,i1-2)
            bplus  = fourthird*U(i0,i1+1) - third*U(i0,i1+2)
            
            if (U(i0-1,i1) .le. U(i0+1,i1)) then
              a  = aminus
            else
              a  = aplus
            endif

            if (U(i0,i1-1) .le. U(i0,i1+1)) then
              b  = bminus
            else
              b  = bplus
            endif
       
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then 
                a  = aplus*U_wall_coef   
                hx = hx*h_wall_coef
              elseif (i0 .eq. dlower0+1) then
                a  = aplus
                hx = hx                
              elseif (i0 .eq. dupper0) then
                a  = aminus*U_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0-1) then
                a  = aminus
                hx = hx
              elseif (i1 .eq. dlower1) then 
                b  = bplus*U_wall_coef   
                hy = hy*h_wall_coef
              elseif (i1 .eq. dlower1+1) then 
                b  = bplus
                hy = hy                 
              elseif (i1 .eq. dupper1) then
                b  = bminus*U_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1-1) then
                b  = bminus
                hy = hy
              endif
            endif

           
            if (b-a .gt. hx) then
              dbar = a + hx
            elseif (a-b .gt. hy) then
              dbar = b + hy
            else
              Q = hx*hx + hy*hy
              R = -2.d0*(hy*hy*a + hx*hx*b)
              S = hy*hy*a*a + hx*hx*b*b - hx*hx*hy*hy
              dbar = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
            endif
            
            U(i0,i1) = dmin1(U(i0,i1),dbar)

         enddo
      enddo

      return
      end
c
c