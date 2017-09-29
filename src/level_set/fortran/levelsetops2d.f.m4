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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine fastsweep1storder2d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dlower0,dupper0
     &     dlower1,dupper1,
     &     dx,
     &     patch_touches_bdry,
     &     consider_bdry_wall)
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
      REAL    hx,hy
      REAL    Q,R,S
      REAL    dbar
      REAL    u_wall_coef
      REAL    h_wall_coef

      if (consider_bdry_wall .eq. 1) then
        u_wall_coef   = 0.d0
        h_wall_coef   = 0.5d0
      else
        u_wall_coef   = 1.d0
        h_wall_coef   = 1.0d0
      endif
     
      
c     Do the four sweeping directions.
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            hx = dx(0)
            hy = dx(1)
            a  = dmin1(U(i0-1,i1),U(i0+1,i1))
            b  = dmin1(U(i0,i1-1),U(i0,i1+1))
 
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then 
                a  = U(i0+1,i1)*u_wall_coef   
                hx = hx*h_wall_coef                
              elseif (i0 .eq. dupper0) then
                a  = U(i0-1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i1 .eq. dlower1) then 
                b  = U(i0,i1+1)*u_wall_coef   
                hy = hy*h_wall_coef                
              elseif (i1 .eq. dupper1) then
                b  = U(i0,i1-1)*u_wall_coef
                hy = hy*h_wall_coef
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
            hx = dx(0)
            hy = dx(1)
            a  = dmin1(U(i0-1,i1),U(i0+1,i1))
            b  = dmin1(U(i0,i1-1),U(i0,i1+1))
            
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then
                a  = U(i0+1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0) then
                a  = U(i0-1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i1 .eq. dlower1) then
                b  = U(i0,i1+1)*u_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1) then
                b  = U(i0,i1-1)*u_wall_coef
                hy = hy*h_wall_coef
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
            hx = dx(0)
            hy = dx(1)
            a  = dmin1(U(i0-1,i1),U(i0+1,i1))
            b  = dmin1(U(i0,i1-1),U(i0,i1+1))
            
c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then
                a  = U(i0+1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0) then
                a  = U(i0-1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i1 .eq. dlower1) then
                b  = U(i0,i1+1)*u_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1) then
                b  = U(i0,i1-1)*u_wall_coef
                hy = hy*h_wall_coef
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
            hx = dx(0)
            hy = dx(1)
            a  = dmin1(U(i0-1,i1),U(i0+1,i1))
            b  = dmin1(U(i0,i1-1),U(i0,i1+1))

c           Take care of physical boundaries.
            if (patch_touches_bdry .eq. 1) then
              if (i0 .eq. dlower0) then
                a  = U(i0+1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i0 .eq. dupper0) then
                a  = U(i0-1,i1)*u_wall_coef
                hx = hx*h_wall_coef
              elseif (i1 .eq. dlower1) then
                b  = U(i0,i1+1)*u_wall_coef
                hy = hy*h_wall_coef
              elseif (i1 .eq. dupper1) then
                b  = U(i0,i1-1)*u_wall_coef
                hy = hy*h_wall_coef
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

