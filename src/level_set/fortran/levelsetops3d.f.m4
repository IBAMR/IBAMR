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
      REAL    a,b,c
      REAL    a1,a2,a3
      REAL    hx,hy,hz
      REAL    h1,h2,h3
      REAL    Q,R,S
      REAL    dtil,dbar
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

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo
      
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0,-1

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0,-1

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0,-1

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0,-1

               hx = dx(0)
               hy = dx(1)
               hz = dx(2)
               a  = dmin1(U(i0-1,i1,i2),U(i0+1,i1,i2))
               b  = dmin1(U(i0,i1-1,i2),U(i0,i1+1,i2))
               c  = dmin1(U(i0,i1,i2-1),U(i0,i1,i2+1))
 
c              Take care of physical boundaries.
               if (patch_touches_bdry .eq. 1) then
                 if (i0 .eq. dlower0) then 
                   a  = U(i0+1,i1,i2)*U_wall_coef   
                   hx = hx*h_wall_coef                
                 elseif (i0 .eq. dupper0) then
                   a  = U(i0-1,i1,i2)*U_wall_coef
                   hx = hx*h_wall_coef
                 elseif (i1 .eq. dlower1) then 
                   b  = U(i0,i1+1,i2)*U_wall_coef   
                   hy = hy*h_wall_coef                
                 elseif (i1 .eq. dupper1) then
                   b  = U(i0,i1-1,i2)*U_wall_coef
                   hy = hy*h_wall_coef
                 elseif (i2 .eq. dlower2) then
                   c = U(i0,i1,i2+1)*U_wall_coef
                   hz = hz*h_wall_coef
                 elseif (i2 .eq. dupper2) then
                   c = U(i0,i1,i2-1)*U_wall_coef
                   hz = hz*h_wall_coef
                 endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c             Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)

            enddo
         enddo
      enddo

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
      subroutine fastsweep2ndorder3d(
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
      REAL    a,b,c
      REAL    a1,a2,a3
      REAL    aplus,bplus,cplus
      REAL    aminus,bminus,cminus
      REAL    hx,hy,hz
      REAL    h1,h2,h3
      REAL    Q,R,S
      REAL    dtil,dbar
      REAL    U_wall_coef
      REAL    h_wall_coef


      if (consider_bdry_wall .eq. 1) then
        U_wall_coef   = zero
        h_wall_coef   = threefourth
      else
        U_wall_coef   = one
        h_wall_coef   = one
      endif
     
      
c     Do the eight sweeping directions.
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0,-1

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0,-1

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0,-1

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2,-1
         do i1 = ilower1,iupper1,-1
            do i0 = ilower0,iupper0,-1

               hx = twothird*dx(0)
               hy = twothird*dx(1)
               hz = twothird*dx(2)
               aminus = fourthird*U(i0-1,i1,i2) - third*U(i0-2,i1,i2)
               aplus  = fourthird*U(i0+1,i1,i2) - third*U(i0+2,i1,i2)
               bminus = fourthird*U(i0,i1-1,i2) - third*U(i0,i1-2,i2)
               bplus  = fourthird*U(i0,i1+1,i2) - third*U(i0,i1+2,i2)
               cminus = fourthird*U(i0,i1,i2-1) - third*U(i0,i1,i2-2)
               cplus  = fourthird*U(i0,i1,i2+1) - third*U(i0,i1,i2+2)
            
               if (U(i0-1,i1,i2) .le. U(i0+1,i1,i2)) then
                  a  = aminus
               else
                  a  = aplus
               endif

               if (U(i0,i1-1,i2) .le. U(i0,i1+1,i2)) then
                  b  = bminus
               else
                  b  = bplus
               endif

               if (U(i0,i1,i2-1) .le. U(i0,i1,i2+1)) then
                  c = cminus
               else
                  c = cplus
               endif
       
c              Take care of physical boundaries.
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
                  elseif (i2 .eq. dlower2) then 
                     c  = cplus*U_wall_coef   
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dlower2+1) then 
                     c  = cplus
                     hz = hz                 
                  elseif (i2 .eq. dupper2) then
                     c  = cminus*U_wall_coef
                     hz = hz*h_wall_coef
                  elseif (i2 .eq. dupper2-1) then
                     c  = cminus
                     hz = hz
                  endif
               endif

c              Additional sorting step for 3D
               if (a .le. b) then
                  if (a .le. c) then
                    if(b .le. c) then
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
                  if (b .le. c) then
                     if (a .le. c) then
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

c              Algorithm to find dbar
              dtil = a1 + h1
              if (dtil .le. a2) then
                 dbar = dtil
              else
                 Q = h1*h1 + h2*h2
                 R = -2.d0*(h2*h2*a1 + h1*h1*a2)
                 S = h2*h2*a1*a1 + h1*h1*a2*a2 - h1*h1*h2*h2
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
              endif
              
              if (dtil .lt. a3) then
                 dbar = dtil
              else
                 Q = 1.d0/(h1*h1) + 1.d0/(h2*h2) + 1.d0/(h3*h3)
                 R = -2.d0*(a1/(h1*h1)+a2/(h2*h2)+a3/(h3*h3))
                 S = a1*a1/(h1*h1)+a2*a2/(h2*h2)+a3*a3/(h3*h3)-1.d0
                 dtil = (-R + sqrt(R*R-4.d0*Q*S))/(2.d0*Q)
                 dbar = dtil
              endif

              U(i0,i1,i2) = dmin1(U(i0,i1,i2),dbar)
            enddo
         enddo
      enddo

      return
      end
c
c

 