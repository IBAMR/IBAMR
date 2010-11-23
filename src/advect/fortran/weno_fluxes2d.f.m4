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
dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute face-centered WENO fluxes corresponding to the advection
c     equation for a cell-centered advection velocity and a
c     cell-centered scalar-valued quantity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine weno_convective_fluxes2d(
     &     fface0,fface1,fface_gcw,
     &     ffwrd0,ffwrd1,ffwrd_gcw,
     &     frevr0,frevr1,frevr_gcw,
     &     fplus0,fplus1,fplus_gcw,
     &     fminus0,fminus1,fminus_gcw,
     &     F,F_gcw,
     &     Q,Q_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/wenocoefs.i)dnl
c
c     Functions.
c
      REAL minmod2,minmod3
c
c     Input.
c
      INTEGER fface_gcw,ffwrd_gcw,frevr_gcw,fplus_gcw,fminus_gcw
      INTEGER F_gcw,Q_gcw,U_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL Q(CELL2d(ilower,iupper,Q_gcw))
      REAL U(CELL2d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Input/Output.
c
      REAL fface0(FACE2d0(ilower,iupper,fface_gcw))
      REAL fface1(FACE2d1(ilower,iupper,fface_gcw))

      REAL ffwrd0(FACE2d0(ilower,iupper,ffwrd_gcw))
      REAL ffwrd1(FACE2d1(ilower,iupper,ffwrd_gcw))

      REAL frevr0(FACE2d0(ilower,iupper,frevr_gcw))
      REAL frevr1(FACE2d1(ilower,iupper,frevr_gcw))

      REAL fplus0(FACE2d0(ilower,iupper,fplus_gcw))
      REAL fplus1(FACE2d1(ilower,iupper,fplus_gcw))

      REAL fminus0(FACE2d0(ilower,iupper,fminus_gcw))
      REAL fminus1(FACE2d1(ilower,iupper,fminus_gcw))

      REAL F(CELL2d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER dim
      INTEGER ic0,ic1
      REAL alpha,alphaj
      PARAMETER(alpha=33.d0)
      REAL cL,cR
      REAL fL,fR
c
c     Compute the x-face fluxes on the "upper" side of each cell
c     (ic0,ic1).
c
      dim = 0

      do ic1 = ilower1,iupper1
         do ic0 = ilower0-R-1,iupper0+R+1
            F(ic0,ic1) = U(ic0,ic1,dim)*Q(ic0,ic1)
         enddo
      enddo

      call weno_flux_plus02d(
     &     fplus0,fplus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)

      call weno_flux_minus02d(
     &     fminus0,fminus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)

      do ic1 = ilower1,iupper1
         do ic0 = ilower0-2,iupper0+1

            fL = 0.25d0*fplus0 (ic0+1,ic1)
            fR = 0.25d0*fminus0(ic0+1,ic1)

            ffwrd0(ic0+1,ic1) =
     &           (fL+fR)+sign(1.d0,+U(ic0  ,ic1,dim))*(fL-fR) +
     &           (fL+fR)+sign(1.d0,+U(ic0+1,ic1,dim))*(fL-fR)

            frevr0(ic0+1,ic1) =
     &           (fL+fR)+sign(1.d0,-U(ic0  ,ic1,dim))*(fL-fR) +
     &           (fL+fR)+sign(1.d0,-U(ic0+1,ic1,dim))*(fL-fR)

         enddo
      enddo

      do ic1 = ilower1,iupper1
         do ic0 = ilower0-1,iupper0

            alphaj = alpha*((
     &           dabs(F(ic0+1,ic1)-2.d0*F(ic0  ,ic1)+F(ic0-1,ic1))/(
     &           dabs(F(ic0+1,ic1)-F(ic0  ,ic1))+
     &           dabs(F(ic0  ,ic1)-F(ic0-1,ic1))))**2.d0)

            cL = 0.25d0*minmod3(
     &           0.5d0*alphaj*minmod2(
     &           frevr0(ic0+1,ic1)-ffwrd0(ic0+1,ic1),
     &           frevr0(ic0  ,ic1)-ffwrd0(ic0  ,ic1)),
     &           F(ic0+1,ic1)-ffwrd0(ic0+1,ic1),
     &           frevr0(ic0  ,ic1)-F(ic0-1,ic1))

            alphaj = alpha*((
     &           dabs(F(ic0  ,ic1)-2.d0*F(ic0+1,ic1)+F(ic0+2,ic1))/(
     &           dabs(F(ic0  ,ic1)-F(ic0+1,ic1))+
     &           dabs(F(ic0+1,ic1)-F(ic0+2,ic1))))**2.d0)

            cR = 0.25d0*minmod3(
     &           0.5d0*alphaj*minmod2(
     &           frevr0(ic0+1,ic1)-ffwrd0(ic0+1,ic1),
     &           frevr0(ic0+2,ic1)-ffwrd0(ic0+2,ic1)),
     &           F(ic0  ,ic1)-ffwrd0(ic0+1,ic1),
     &           frevr0(ic0+2,ic1)-F(ic0+2,ic1))

            fface0(ic0+1,ic1) = ffwrd0(ic0+1,ic1) +
     &           (cL+cR)+sign(1.d0,+U(ic0  ,ic1,dim))*(cL-cR) +
     &           (cL+cR)+sign(1.d0,+U(ic0+1,ic1,dim))*(cL-cR)

         enddo
      enddo
c
c     Compute the y-face fluxes on the "upper" side of each cell
c     (ic0,ic1).
c
      dim = 1

      do ic1 = ilower1-R-1,iupper1+R+1
         do ic0 = ilower0,iupper0
            F(ic0,ic1) = U(ic0,ic1,dim)*Q(ic0,ic1)
         enddo
      enddo

      call weno_flux_plus12d(
     &     fplus1,fplus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)

      call weno_flux_minus12d(
     &     fminus1,fminus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)

      do ic0 = ilower0,iupper0
         do ic1 = ilower1-2,iupper1+1

            fL = 0.25d0*fplus1 (ic1+1,ic0)
            fR = 0.25d0*fminus1(ic1+1,ic0)

            ffwrd1(ic1+1,ic0) =
     &           (fL+fR)+sign(1.d0,+U(ic0,ic1  ,dim))*(fL-fR) +
     &           (fL+fR)+sign(1.d0,+U(ic0,ic1+1,dim))*(fL-fR)

            frevr1(ic1+1,ic0) =
     &           (fL+fR)+sign(1.d0,-U(ic0,ic1  ,dim))*(fL-fR) +
     &           (fL+fR)+sign(1.d0,-U(ic0,ic1+1,dim))*(fL-fR)

         enddo
      enddo

      do ic0 = ilower0,iupper0
         do ic1 = ilower1-1,iupper1

            alphaj = alpha*((
     &           dabs(F(ic0,ic1+1)-2.d0*F(ic0,ic1  )+F(ic0,ic1-1))/(
     &           dabs(F(ic0,ic1+1)-F(ic0,ic1  ))+
     &           dabs(F(ic0,ic1  )-F(ic0,ic1-1))))**2.d0)

            cL = 0.25d0*minmod3(
     &           0.5d0*alphaj*minmod2(
     &           frevr1(ic1+1,ic0)-ffwrd1(ic1+1,ic0),
     &           frevr1(ic1  ,ic0)-ffwrd1(ic1  ,ic0)),
     &           F(ic0,ic1+1)-ffwrd1(ic1+1,ic0),
     &           frevr1(ic1  ,ic0)-F(ic0,ic1-1))

            alphaj = alpha*((
     &           dabs(F(ic0,ic1  )-2.d0*F(ic0,ic1+1)+F(ic0,ic1+2))/(
     &           dabs(F(ic0,ic1  )-F(ic0,ic1+1))+
     &           dabs(F(ic0,ic1+1)-F(ic0,ic1+2))))**2.d0)

            cR = 0.25d0*minmod3(
     &           0.5d0*alphaj*minmod2(
     &           frevr1(ic1+1,ic0)-ffwrd1(ic1+1,ic0),
     &           frevr1(ic1+2,ic0)-ffwrd1(ic1+2,ic0)),
     &           F(ic0,ic1  )-ffwrd1(ic1+1,ic0),
     &           frevr1(ic1+2,ic0)-F(ic0,ic1+2))

            fface1(ic1+1,ic0) = ffwrd1(ic1+1,ic0) +
     &           (cL+cR)+sign(1.d0,+U(ic0,ic1  ,dim))*(cL-cR) +
     &           (cL+cR)+sign(1.d0,+U(ic0,ic1+1,dim))*(cL-cR)

         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute positive (right-going) x face-centered WENO fluxes
c     corresponding to the advection equation for a cell-centered
c     advection velocity and a cell-centered scalar-valued quantity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine weno_flux_plus02d(
     &     fplus0,fplus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/wenocoefs.i)dnl
c
c     Input.
c
      INTEGER fplus_gcw,F_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL F(CELL2d(ilower,iupper,F_gcw))
c
c     Input/Output.
c
      REAL fplus0(FACE2d0(ilower,iupper,fplus_gcw))
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER k,l,m
      REAL IS(0:R-1),IS_sub
      REAL alpha(0:R-1),alpha_sum
      REAL omega(0:R-1)
c
c     Compute the positive x face-centered fluxes on the "upper" side of
c     each cell (ic0,ic1).
c
      do ic1 = ilower1,iupper1
         do ic0 = ilower0-2,iupper0+1

            do k = 0,R-1
               IS(k) = 0.d0
               do m = 1,R-1
                  IS_sub = 0.d0
                  do l = 0,R-1
                     IS_sub = IS_sub + D3(k,l,m)*F(ic0+k+l-R+1,ic1)
                  enddo
                  IS(k) = IS(k) + IS_sub*IS_sub
               enddo
               alpha(k) = C3(k)/((EPSILON + IS(k))*(EPSILON + IS(k)))
            enddo

            alpha_sum = 0.d0
            do k = 0,R-1
               alpha_sum = alpha_sum + alpha(k)
            enddo
            do k = 0,R-1
               omega(k) = alpha(k)/alpha_sum
            enddo

            fplus0(ic0+1,ic1) = 0.d0
            do k = 0,R-1
               do l = 0,R-1
                  fplus0(ic0+1,ic1) = fplus0(ic0+1,ic1) +
     &                 omega(k)*A3(k,l)*F(ic0+k+l-R+1,ic1)
               enddo
            enddo

         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute negative (left-going) x face-centered WENO fluxes
c     corresponding to the advection equation for a cell-centered
c     advection velocity and a cell-centered scalar-valued quantity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine weno_flux_minus02d(
     &     fminus0,fminus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/wenocoefs.i)dnl
c
c     Input.
c
      INTEGER fminus_gcw,F_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL F(CELL2d(ilower,iupper,F_gcw))
c
c     Input/Output.
c
      REAL fminus0(FACE2d0(ilower,iupper,fminus_gcw))
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER k,l,m
      REAL IS(0:R-1),IS_sub
      REAL alpha(0:R-1),alpha_sum
      REAL omega(0:R-1)
c
c     Compute the negative x face-centered fluxes on the "upper" side of
c     each cell (ic0,ic1).
c
      do ic1 = ilower1,iupper1
         do ic0 = ilower0-2,iupper0+1

            do k = 0,R-1
               IS(k) = 0.d0
               do m = 1,R-1
                  IS_sub = 0.d0
                  do l = 0,R-1
                     IS_sub = IS_sub + D3(k,l,m)*F(ic0-k-l+R,ic1)
                  enddo
                  IS(k) = IS(k) + IS_sub*IS_sub
               enddo
               alpha(k) = C3(k)/((EPSILON + IS(k))*(EPSILON + IS(k)))
            enddo

            alpha_sum = 0.d0
            do k = 0,R-1
               alpha_sum = alpha_sum + alpha(k)
            enddo
            do k = 0,R-1
               omega(k) = alpha(k)/alpha_sum
            enddo

            fminus0(ic0+1,ic1) = 0.d0
            do k = 0,R-1
               do l = 0,R-1
                  fminus0(ic0+1,ic1) = fminus0(ic0+1,ic1) +
     &                 omega(k)*A3(k,l)*F(ic0-k-l+R,ic1)
               enddo
            enddo

         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute positive (right-going) y face-centered WENO fluxes
c     corresponding to the advection equation for a cell-centered
c     advection velocity and a cell-centered scalar-valued quantity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine weno_flux_plus12d(
     &     fplus1,fplus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/wenocoefs.i)dnl
c
c     Input.
c
      INTEGER fplus_gcw,F_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL F(CELL2d(ilower,iupper,F_gcw))
c
c     Input/Output.
c
      REAL fplus1(FACE2d1(ilower,iupper,fplus_gcw))
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER k,l,m
      REAL IS(0:R-1),IS_sub
      REAL alpha(0:R-1),alpha_sum
      REAL omega(0:R-1)
c
c     Compute the positive y face-centered fluxes on the "upper" side of
c     each cell (ic0,ic1).
c
      do ic0 = ilower0,iupper0
         do ic1 = ilower1-2,iupper1+1

            do k = 0,R-1
               IS(k) = 0.d0
               do m = 1,R-1
                  IS_sub = 0.d0
                  do l = 0,R-1
                     IS_sub = IS_sub + D3(k,l,m)*F(ic0,ic1+k+l-R+1)
                  enddo
                  IS(k) = IS(k) + IS_sub*IS_sub
               enddo
               alpha(k) = C3(k)/((EPSILON + IS(k))*(EPSILON + IS(k)))
            enddo

            alpha_sum = 0.d0
            do k = 0,R-1
               alpha_sum = alpha_sum + alpha(k)
            enddo
            do k = 0,R-1
               omega(k) = alpha(k)/alpha_sum
            enddo

            fplus1(ic1+1,ic0) = 0.d0
            do k = 0,R-1
               do l = 0,R-1
                  fplus1(ic1+1,ic0) = fplus1(ic1+1,ic0) +
     &                 omega(k)*A3(k,l)*F(ic0,ic1+k+l-R+1)
               enddo
            enddo

         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute negative (left-going) y face-centered WENO fluxes
c     corresponding to the advection equation for a cell-centered
c     advection velocity and a cell-centered scalar-valued quantity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine weno_flux_minus12d(
     &     fminus1,fminus_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
include(TOP_SRCDIR/src/advect/fortran/wenocoefs.i)dnl
c
c     Input.
c
      INTEGER fminus_gcw,F_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL F(CELL2d(ilower,iupper,F_gcw))
c
c     Input/Output.
c
      REAL fminus1(FACE2d1(ilower,iupper,fminus_gcw))
c
c     Local variables.
c
      INTEGER ic0,ic1
      INTEGER k,l,m
      REAL IS(0:R-1),IS_sub
      REAL alpha(0:R-1),alpha_sum
      REAL omega(0:R-1)
c
c     Compute the negative y face-centered fluxes on the "upper" side of
c     each cell (ic0,ic1).
c
      do ic0 = ilower0,iupper0
         do ic1 = ilower1-2,iupper1+1

            do k = 0,R-1
               IS(k) = 0.d0
               do m = 1,R-1
                  IS_sub = 0.d0
                  do l = 0,R-1
                     IS_sub = IS_sub + D3(k,l,m)*F(ic0,ic1-k-l+R)
                  enddo
                  IS(k) = IS(k) + IS_sub*IS_sub
               enddo
               alpha(k) = C3(k)/((EPSILON + IS(k))*(EPSILON + IS(k)))
            enddo

            alpha_sum = 0.d0
            do k = 0,R-1
               alpha_sum = alpha_sum + alpha(k)
            enddo
            do k = 0,R-1
               omega(k) = alpha(k)/alpha_sum
            enddo

            fminus1(ic1+1,ic0) = 0.d0
            do k = 0,R-1
               do l = 0,R-1
                  fminus1(ic1+1,ic0) = fminus1(ic1+1,ic0) +
     &                 omega(k)*A3(k,l)*F(ic0,ic1-k-l+R)
               enddo
            enddo

         enddo
      enddo
c
      return
      end
c
