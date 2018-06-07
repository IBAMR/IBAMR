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
c     Interpolate the components of a staggered velocity field onto the
c     faces of zones centered about the components of the velocity
c     field.
c
c     NOTES:
c
c     U0 and U1 are standard side-centered staggered grid velocities for
c     the patch [(ifirst0,ilast0),(ifirst1,ilast1)].
c
c     V00 and V01 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     V10 and V11 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_interp_comps2d(
     &     patch_ifirst0,patch_ilast0,
     &     patch_ifirst1,patch_ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_V0_gc0,n_V0_gc1,
     &     V00,V01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_V1_gc0,n_V1_gc1,
     &     V10,V11)
c
      implicit none
c
c     Input.
c
      INTEGER patch_ifirst0,patch_ilast0
      INTEGER patch_ifirst1,patch_ilast1

      INTEGER n_U_gc0,n_U_gc1

      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_V0_gc0,n_V0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_V1_gc0,n_V1_gc1

      REAL U0(
     &     SIDE2d0VECG(patch_ifirst,patch_ilast,n_U_gc)
     &     )
      REAL U1(
     &     SIDE2d1VECG(patch_ifirst,patch_ilast,n_U_gc)
     &     )
c
c     Input/Output.
c
      REAL V00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
      REAL V11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
c
c     Interpolate the components of the velocity at each zone face.
c
      gc0 = min(n_U_gc0-1,n_V0_gc0)
      gc1 = min(n_U_gc1  ,n_V0_gc1)

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            V00(i0,i1) = 0.5d0*(U0(i0-1,i1)+U0(i0,i1))
            V01(i1,i0) = 0.5d0*(U1(i0-1,i1)+U1(i0,i1))
         enddo
      enddo

      gc0 = min(n_U_gc0  ,n_V1_gc0)
      gc1 = min(n_U_gc1-1,n_V1_gc1)

      do    i0 = side1_ifirst0-gc0,side1_ilast0+gc0
         do i1 = side1_ifirst1-gc1,side1_ilast1+gc1
            V10(i0,i1) = 0.5d0*(U0(i0,i1-1)+U0(i0,i1))
            V11(i1,i0) = 0.5d0*(U1(i0,i1-1)+U1(i0,i1))
         enddo
      enddo
c
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Upwind the components of side-centered quantity onto the
c     faces of zones centered about the components of the velocity
c     field.
c
c     NOTES:
c
c     R0 and R1 are standard side-centered staggered grid quantities for
c     the patch [(ifirst0,ilast0),(ifirst1,ilast1)].
c
c     V00 and V01 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     V10 and V11 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
c     R00 and R01 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     R10 and R11 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vc_navier_stokes_upwind_quantity2d(
     &     patch_ifirst0,patch_ilast0,
     &     patch_ifirst1,patch_ilast1,
     &     n_R_gc0,n_R_gc1,
     &     R0,R1,
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_V0_gc0,n_V0_gc1,
     &     V00,V01,
     &     n_R0_gc0,n_R0_gc1,
     &     R00,R01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_V1_gc0,n_V1_gc1,
     &     V10,V11,
     &     n_R1_gc0,n_R1_gc1,
     &     R10,R11)
c
      implicit none
c
c     Input.
c
      INTEGER patch_ifirst0,patch_ilast0
      INTEGER patch_ifirst1,patch_ilast1

      INTEGER n_R_gc0,n_R_gc1

      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_V0_gc0,n_V0_gc1
      INTEGER n_R0_gc0,n_R0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_V1_gc0,n_V1_gc1
      INTEGER n_R1_gc0,n_R1_gc1

      REAL R0(
     &     SIDE2d0VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL R1(
     &     SIDE2d1VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL V00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
      REAL V11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
c
c     Input/Output.
c
      REAL R00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
      REAL R11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
c
c     Compute the upwinded quantity at each zone face.
c
      gc0 = n_R0_gc0
      gc1 = n_R0_gc1

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            if (V00(i0,i1) .ge. 0.d0) then
              R00(i0,i1) = R0(i0-1,i1)
            else
              R00(i0,i1) = R0(i0,i1)
            endif

            if (V01(i1,i0) .ge. 0.d0) then
              R01(i1,i0) = R0(i0,i1-1)
            else
              R01(i1,i0) = R0(i0,i1)
            endif

         enddo
      enddo

      gc0 = n_R1_gc0
      gc1 = n_R1_gc1

      do    i0 = side1_ifirst0-gc0,side1_ilast0+gc0
         do i1 = side1_ifirst1-gc1,side1_ilast1+gc1
            if (V10(i0,i1) .ge. 0.d0) then
              R10(i0,i1) = R1(i0-1,i1)
            else
              R10(i0,i1) = R1(i0,i1)
            endif

            if (V11(i1,i0) .ge. 0.d0) then
              R11(i1,i0) = R1(i0,i1-1)
            else
              R11(i1,i0) = R1(i0,i1)
            endif

         enddo
      enddo
c
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply cubic interpolation upwinding to place the components of
c     side-centered quantity onto the faces of zones centered about the
c     components of the velocity field.
c
c     Note that this scheme carries out nonlinear upwinding.
c
c     NOTES:
c
c     R0 and R1 are standard side-centered staggered grid quantities for
c     the patch [(ifirst0,ilast0),(ifirst1,ilast1)].
c
c     V00 and V01 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     V10 and V11 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
c     R00 and R01 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     R10 and R11 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine vc_navier_stokes_cui_quantity2d(
     &     patch_ifirst0,patch_ilast0,
     &     patch_ifirst1,patch_ilast1,
     &     n_R_gc0,n_R_gc1,
     &     R0,R1,
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_V0_gc0,n_V0_gc1,
     &     V00,V01,
     &     n_R0_gc0,n_R0_gc1,
     &     R00,R01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_V1_gc0,n_V1_gc1,
     &     V10,V11,
     &     n_R1_gc0,n_R1_gc1,
     &     R10,R11)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER patch_ifirst0,patch_ilast0
      INTEGER patch_ifirst1,patch_ilast1

      INTEGER n_R_gc0,n_R_gc1

      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_V0_gc0,n_V0_gc1
      INTEGER n_R0_gc0,n_R0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_V1_gc0,n_V1_gc1
      INTEGER n_R1_gc0,n_R1_gc1

      REAL R0(
     &     SIDE2d0VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL R1(
     &     SIDE2d1VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL V00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
      REAL V11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )

c
c     Input/Output.
c
      REAL R00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
      REAL R11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
      REAL RC,RU,RD
      REAL Rf_HR
c
c     Compute the cubic upwind interpolation of quantity at each zone face.
c
      gc0 = n_R0_gc0
      gc1 = n_R0_gc1

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            if (V00(i0,i1) .ge. 0.d0) then
              RC  = R0(i0-1,i1)
              RU  = R0(i0-2,i1)
              RD  = R0(i0,i1)
            else
              RC  = R0(i0,i1)
              RU  = R0(i0+1,i1)
              RD  = R0(i0-1,i1)
            endif

c           High-resolution scheme (HR)
            call interpolate_cui_hr_quantity2d(RU,RC,RD,Rf_HR)
            R00(i0,i1) = Rf_HR

            if (V01(i1,i0) .ge. 0.d0) then
              RC  = R0(i0,i1-1)
              RU  = R0(i0,i1-2)
              RD  = R0(i0,i1)
            else
              RC  = R0(i0,i1)
              RU  = R0(i0,i1+1)
              RD  = R0(i0,i1-1)
            endif

c           High-resolution scheme (HR)
            call interpolate_cui_hr_quantity2d(RU,RC,RD,Rf_HR)
            R01(i1,i0) = Rf_HR
         enddo
      enddo

      gc0 = n_R1_gc0
      gc1 = n_R1_gc1

      do    i0 = side1_ifirst0-gc0,side1_ilast0+gc0
         do i1 = side1_ifirst1-gc1,side1_ilast1+gc1
            if (V10(i0,i1) .ge. 0.d0) then
               RC  = R1(i0-1,i1)
               RU  = R1(i0-2,i1)
               RD  = R1(i0,i1)
            else
               RC  = R1(i0,i1)
               RU  = R1(i0+1,i1)
               RD  = R1(i0-1,i1)
            endif

c           High-resolution scheme (HR)
            call interpolate_cui_hr_quantity2d(RU,RC,RD,Rf_HR)
            R10(i0,i1) = Rf_HR

            if (V11(i1,i0) .ge. 0.d0) then
               RC  = R1(i0,i1-1)
               RU  = R1(i0,i1-2)
               RD  = R1(i0,i1)
            else
               RC  = R1(i0,i1)
               RU  = R1(i0,i1+1)
               RD  = R1(i0,i1-1)
            endif

c           High-resolution scheme (HR)
            call interpolate_cui_hr_quantity2d(RU,RC,RD,Rf_HR)
            R11(i1,i0) = Rf_HR
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply FBICS upwinding to place the components of
c     side-centered quantity onto the faces of zones centered about the
c     components of the velocity field.
c
c     Note that this scheme carries out nonlinear upwinding.
c
c     NOTES:
c
c     R0 and R1 are standard side-centered staggered grid quantities for
c     the patch [(ifirst0,ilast0),(ifirst1,ilast1)].
c
c     V00 and V01 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     V10 and V11 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
c     R00 and R01 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     R10 and R11 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine vc_navier_stokes_fbics_quantity2d(
     &     patch_ifirst0,patch_ilast0,
     &     patch_ifirst1,patch_ilast1,
     &     n_R_gc0,n_R_gc1,
     &     R0,R1,
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_V0_gc0,n_V0_gc1,
     &     V00,V01,
     &     n_R0_gc0,n_R0_gc1,
     &     R00,R01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_V1_gc0,n_V1_gc1,
     &     V10,V11,
     &     n_R1_gc0,n_R1_gc1,
     &     R10,R11)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER patch_ifirst0,patch_ilast0
      INTEGER patch_ifirst1,patch_ilast1

      INTEGER n_R_gc0,n_R_gc1

      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_V0_gc0,n_V0_gc1
      INTEGER n_R0_gc0,n_R0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_V1_gc0,n_V1_gc1
      INTEGER n_R1_gc0,n_R1_gc1

      REAL R0(
     &     SIDE2d0VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL R1(
     &     SIDE2d1VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL V00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
      REAL V11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )

c
c     Input/Output.
c
      REAL R00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
      REAL R11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
      REAL RC,RU,RD
      REAL Rf_HR
c
c     Compute the cubic upwind interpolation of quantity at each zone face.
c
      gc0 = n_R0_gc0
      gc1 = n_R0_gc1

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            if (V00(i0,i1) .ge. 0.d0) then
              RC  = R0(i0-1,i1)
              RU  = R0(i0-2,i1)
              RD  = R0(i0,i1)
            else
              RC  = R0(i0,i1)
              RU  = R0(i0+1,i1)
              RD  = R0(i0-1,i1)
            endif

c           High-resolution scheme (HR)
            call interpolate_fbics_hr_quantity2d(RU,RC,RD,Rf_HR)
            R00(i0,i1) = Rf_HR

            if (V01(i1,i0) .ge. 0.d0) then
               RC  = R0(i0,i1-1)
               RU  = R0(i0,i1-2)
               RD  = R0(i0,i1)
            else
               RC  = R0(i0,i1)
               RU  = R0(i0,i1+1)
               RD  = R0(i0,i1-1)
            endif

c           High-resolution scheme (HR)
            call interpolate_fbics_hr_quantity2d(RU,RC,RD,Rf_HR)
            R01(i1,i0) = Rf_HR
         enddo
      enddo

      gc0 = n_R1_gc0
      gc1 = n_R1_gc1

      do    i0 = side1_ifirst0-gc0,side1_ilast0+gc0
         do i1 = side1_ifirst1-gc1,side1_ilast1+gc1
            if (V10(i0,i1) .ge. 0.d0) then
               RC  = R1(i0-1,i1)
               RU  = R1(i0-2,i1)
               RD  = R1(i0,i1)
            else
               RC  = R1(i0,i1)
               RU  = R1(i0+1,i1)
               RD  = R1(i0-1,i1)
            endif

c           High-resolution scheme (HR)
            call interpolate_fbics_hr_quantity2d(RU,RC,RD,Rf_HR)
            R10(i0,i1) = Rf_HR

            if (V11(i1,i0) .ge. 0.d0) then
               RC  = R1(i0,i1-1)
               RU  = R1(i0,i1-2)
               RD  = R1(i0,i1)
            else
               RC  = R1(i0,i1)
               RU  = R1(i0,i1+1)
               RD  = R1(i0,i1-1)
            endif

c           High-resolution scheme (HR)
            call interpolate_fbics_hr_quantity2d(RU,RC,RD,Rf_HR)
            R11(i1,i0) = Rf_HR
         enddo
      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply modified gamma upwinding to place the components of
c     side-centered quantity onto the faces of zones centered about the
c     components of the velocity field.
c
c     Note that this scheme carries out nonlinear upwinding.
c
c     NOTES:
c
c     R0 and R1 are standard side-centered staggered grid quantities for
c     the patch [(ifirst0,ilast0),(ifirst1,ilast1)].
c
c     V00 and V01 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     V10 and V11 are face-centered staggered grid velocities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     velocities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
c     R00 and R01 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the x components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0+1),(ifirst1,ilast1)].
c
c     R10 and R11 are face-centered staggered grid quantities defined at
c     the faces of the control volumes centered about the y components
c     of the side-centered velocity, i.e., face-centered staggered grid
c     densities for the patch [(ifirst0,ilast0),(ifirst1,ilast1+1)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine vc_navier_stokes_mgamma_quantity2d(
     &     patch_ifirst0,patch_ilast0,
     &     patch_ifirst1,patch_ilast1,
     &     n_R_gc0,n_R_gc1,
     &     R0,R1,
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_V0_gc0,n_V0_gc1,
     &     V00,V01,
     &     n_R0_gc0,n_R0_gc1,
     &     R00,R01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_V1_gc0,n_V1_gc1,
     &     V10,V11,
     &     n_R1_gc0,n_R1_gc1,
     &     R10,R11)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
c     Input.
c
      INTEGER patch_ifirst0,patch_ilast0
      INTEGER patch_ifirst1,patch_ilast1

      INTEGER n_R_gc0,n_R_gc1

      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_V0_gc0,n_V0_gc1
      INTEGER n_R0_gc0,n_R0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_V1_gc0,n_V1_gc1
      INTEGER n_R1_gc0,n_R1_gc1

      REAL R0(
     &     SIDE2d0VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL R1(
     &     SIDE2d1VECG(patch_ifirst,patch_ilast,n_R_gc)
     &     )
      REAL V00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_V0_gc)
     &     )
      REAL V10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )
      REAL V11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_V1_gc)
     &     )

c
c     Input/Output.
c
      REAL R00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_R0_gc)
     &     )
      REAL R10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
      REAL R11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_R1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
      REAL RC,RU,RD
      REAL Rf_HR
c
c     Compute the cubic upwind interpolation of quantity at each zone face.
c
      gc0 = n_R0_gc0
      gc1 = n_R0_gc1

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            if (V00(i0,i1) .ge. 0.d0) then
              RC  = R0(i0-1,i1)
              RU  = R0(i0-2,i1)
              RD  = R0(i0,i1)
            else
              RC  = R0(i0,i1)
              RU  = R0(i0+1,i1)
              RD  = R0(i0-1,i1)
            endif

c           High-resolution scheme (HR)
            call interpolate_mgamma_hr_quantity2d(RU,RC,RD,Rf_HR)
            R00(i0,i1) = Rf_HR

            if (V01(i1,i0) .ge. 0.d0) then
               RC  = R0(i0,i1-1)
               RU  = R0(i0,i1-2)
               RD  = R0(i0,i1)
            else
               RC  = R0(i0,i1)
               RU  = R0(i0,i1+1)
               RD  = R0(i0,i1-1)
            endif

c           High-resolution scheme (HR)
            call interpolate_mgamma_hr_quantity2d(RU,RC,RD,Rf_HR)
            R01(i1,i0) = Rf_HR
         enddo
      enddo

      gc0 = n_R1_gc0
      gc1 = n_R1_gc1

      do    i0 = side1_ifirst0-gc0,side1_ilast0+gc0
         do i1 = side1_ifirst1-gc1,side1_ilast1+gc1
            if (V10(i0,i1) .ge. 0.d0) then
               RC  = R1(i0-1,i1)
               RU  = R1(i0-2,i1)
               RD  = R1(i0,i1)
            else
               RC  = R1(i0,i1)
               RU  = R1(i0+1,i1)
               RD  = R1(i0-1,i1)
            endif

c           High-resolution scheme (HR)
            call interpolate_mgamma_hr_quantity2d(RU,RC,RD,Rf_HR)
            R10(i0,i1) = Rf_HR

            if (V11(i1,i0) .ge. 0.d0) then
               RC  = R1(i0,i1-1)
               RU  = R1(i0,i1-2)
               RD  = R1(i0,i1)
            else
               RC  = R1(i0,i1)
               RU  = R1(i0,i1+1)
               RD  = R1(i0,i1-1)
            endif

c           High-resolution scheme (HR)
            call interpolate_mgamma_hr_quantity2d(RU,RC,RD,Rf_HR)
            R11(i1,i0) = Rf_HR
         enddo
      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interpolate face values based on the CUI high resolution scheme
c     qf will be set as some function of qU, qC, and qD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine interpolate_cui_hr_quantity2d(qU,qC,qD,qf)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      REAL qU,qC,qD
c
c     Input/Output.
c
      REAL qf

c
c     Local variables.
c
      REAL ac,af

c     High-resolution scheme (HR)
      ac = (qC - qU)/(qD - qU)

      if (qD - qU .eq. zero) then
c       ac will be NaN, but qf = qU, so af can be any number
        af = zero
      else if (zero .lt. ac .and. ac .le. 2.d0/13.d0) then
        af = three*ac
      else if (2.d0/13.d0 .lt. ac .and. ac .le. 4.d0/5.d0) then
        af = ac*5.d0/6.d0 + third
      else if (4.d0/5.d0 .lt. ac .and. ac .le. one) then
        af = one
      else
        af = ac
      endif
      qf = af*(qD - qU) + qU

      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interpolate face values based on the CUI compressive scheme
c     qf will be set as some function of qU, qC, and qD
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine interpolate_cui_bd_quantity2d(qU,qC,qD,qf)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      REAL qU,qC,qD
c
c     Input/Output.
c
      REAL qf

c
c     Local variables.
c
      REAL ac,af

c     High-resolution scheme (HR)
      ac = (qC - qU)/(qD - qU)

      if (qD - qU .eq. zero) then
c       ac will be NaN, but qf = qU, so af can be any number
        af = zero
      else if (zero .lt. ac .and. ac .le. third) then
        af = three*ac
      else if (third .lt. ac .and. ac .le. one) then
        af = one
      else
        af = ac
      endif
      qf = af*(qD - qU) + qU

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interpolate face values based on the FBICS high resolution scheme
c     qf will be set as some function of qU, qC, and qD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine interpolate_fbics_hr_quantity2d(qU,qC,qD,qf)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      REAL qU,qC,qD
c
c     Input/Output.
c
      REAL qf

c
c     Local variables.
c
      REAL ac,af

c     High-resolution scheme (HR)
      ac = (qC - qU)/(qD - qU)

      if (qD - qU .eq. zero) then
c       ac will be NaN, but qf = qU, so af can be any number
        af = zero
      else if (zero .lt. ac .and. ac .le. eighth) then
        af = three*ac
      else if (eighth .lt. ac .and. ac .le. threefourth) then
        af = ac + fourth
      else if (threefourth .lt. ac .and. ac .le. one) then
        af = one
      else
        af = ac
      endif
      qf = af*(qD - qU) + qU

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interpolate face values based on the M-Gamma high resolution scheme
c     qf will be set as some function of qU, qC, and qD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine interpolate_mgamma_hr_quantity2d(qU,qC,qD,qf)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl

c
c     Input.
c
      REAL qU,qC,qD
c
c     Input/Output.
c
      REAL qf

c
c     Local variables.
c
      REAL ac,af

c     High-resolution scheme (HR)
      ac = (qC - qU)/(qD - qU)

      if (qD - qU .eq. zero) then
c       ac will be NaN, but qf = qU, so af can be any number
        af = zero
      else if (zero .lt. ac .and. ac .le. fourth) then
        af = 5.d0/2.d0*ac
      else if (fourth .lt. ac .and. ac .le. one) then
        af = half*ac + half
      else
        af = ac
      endif
      qf = af*(qD - qU) + qU

      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Reset the face-centered advection velocity about the control
c     volumes for each component of the velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_reset_adv_velocity2d(
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_U_adv0_gc0,n_U_adv0_gc1,
     &     U_adv00,U_adv01,
     &     n_U_half0_gc0,n_U_half0_gc1,
     &     U_half00,U_half01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_U_adv1_gc0,n_U_adv1_gc1,
     &     U_adv10,U_adv11,
     &     n_U_half1_gc0,n_U_half1_gc1,
     &     U_half10,U_half11)
c
      implicit none
c
c     Input.
c
      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_U_adv0_gc0,n_U_adv0_gc1
      INTEGER n_U_half0_gc0,n_U_half0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_U_adv1_gc0,n_U_adv1_gc1
      INTEGER n_U_half1_gc0,n_U_half1_gc1

      REAL U_half00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_U_half0_gc)
     &     )
      REAL U_half01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_U_half0_gc)
     &     )
      REAL U_half10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_U_half1_gc)
     &     )
      REAL U_half11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_U_half1_gc)
     &     )
c
c     Input/Output.
c
      REAL U_adv00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_U_adv0_gc)
     &     )
      REAL U_adv01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_U_adv0_gc)
     &     )
      REAL U_adv10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_U_adv1_gc)
     &     )
      REAL U_adv11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_U_adv1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
c
c     Reset the advection velocity.
c
      gc0 = min(n_U_adv0_gc0,n_U_adv1_gc0,n_U_half0_gc0,n_U_half1_gc0)
      gc1 = min(n_U_adv0_gc1,n_U_adv1_gc1,n_U_half0_gc1,n_U_half1_gc1)

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            U_adv00(i0,i1) = U_half00(i0,i1)
            U_adv01(i1,i0) = U_half10(i0,i1)
         enddo
      enddo

      do    i1 = side1_ifirst1-gc1,side1_ilast1+gc1
         do i0 = side1_ifirst0-gc0,side1_ilast0+gc0
            U_adv10(i0,i1) = U_half01(i1,i0)
            U_adv11(i1,i0) = U_half11(i1,i0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the face-centered momentum about the control
c     volumes for each component of the velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vc_navier_stokes_compute_momentum2d(
     &     side0_ifirst0,side0_ilast0,
     &     side0_ifirst1,side0_ilast1,
     &     n_P_half0_gc0,n_P_half0_gc1,
     &     P_half00,P_half01,
     &     n_R_half0_gc0,n_R_half0_gc1,
     &     R_half00,R_half01,
     &     n_U_half0_gc0,n_U_half0_gc1,
     &     U_half00,U_half01,
     &     side1_ifirst0,side1_ilast0,
     &     side1_ifirst1,side1_ilast1,
     &     n_P_half1_gc0,n_P_half1_gc1,
     &     P_half10,P_half11,
     &     n_R_half1_gc0,n_R_half1_gc1,
     &     R_half10,R_half11,
     &     n_U_half1_gc0,n_U_half1_gc1,
     &     U_half10,U_half11)
c
      implicit none
c
c     Input.
c
      INTEGER side0_ifirst0,side0_ilast0
      INTEGER side0_ifirst1,side0_ilast1

      INTEGER n_P_half0_gc0,n_P_half0_gc1
      INTEGER n_R_half0_gc0,n_R_half0_gc1
      INTEGER n_U_half0_gc0,n_U_half0_gc1

      INTEGER side1_ifirst0,side1_ilast0
      INTEGER side1_ifirst1,side1_ilast1

      INTEGER n_P_half1_gc0,n_P_half1_gc1
      INTEGER n_R_half1_gc0,n_R_half1_gc1
      INTEGER n_U_half1_gc0,n_U_half1_gc1

      REAL U_half00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_U_half0_gc)
     &     )
      REAL U_half01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_U_half0_gc)
     &     )
      REAL U_half10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_U_half1_gc)
     &     )
      REAL U_half11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_U_half1_gc)
     &     )

      REAL R_half00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_R_half0_gc)
     &     )
      REAL R_half01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_R_half0_gc)
     &     )
      REAL R_half10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_R_half1_gc)
     &     )
      REAL R_half11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_R_half1_gc)
     &     )
c
c     Input/Output.
c
      REAL P_half00(
     &     FACE2d0VECG(side0_ifirst,side0_ilast,n_P_half0_gc)
     &     )
      REAL P_half01(
     &     FACE2d1VECG(side0_ifirst,side0_ilast,n_P_half0_gc)
     &     )
      REAL P_half10(
     &     FACE2d0VECG(side1_ifirst,side1_ilast,n_P_half1_gc)
     &     )
      REAL P_half11(
     &     FACE2d1VECG(side1_ifirst,side1_ilast,n_P_half1_gc)
     &     )
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gc0,gc1
c
c     Compute the momentum on the faces.
c
      gc0 = min(n_P_half0_gc0,n_P_half1_gc0)
      gc1 = min(n_P_half0_gc1,n_P_half1_gc1)

      do    i1 = side0_ifirst1-gc1,side0_ilast1+gc1
         do i0 = side0_ifirst0-gc0,side0_ilast0+gc0
            P_half00(i0,i1) = R_half00(i0,i1)*U_half00(i0,i1)
            P_half01(i1,i0) = R_half01(i1,i0)*U_half01(i1,i0)
         enddo
      enddo

      do    i1 = side1_ifirst1-gc1,side1_ilast1+gc1
         do i0 = side1_ifirst0-gc0,side1_ilast0+gc0
            P_half10(i0,i1) = R_half10(i0,i1)*U_half10(i0,i1)
            P_half11(i1,i0) = R_half11(i1,i0)*U_half11(i1,i0)
         enddo
      enddo
c
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes an updated density field using the three stage formula
c     R = a0*R0 + a1*R1 - a2*dt*div[rhalf*u]
c
c     a0,a1,a2 are coefficients for steps of a SSP Runge-Kutta update
c     R is a side-centered updated density field
c     R0,R1 are side-centered density fields from different RK stages
c     rhalf is the face-centered interpolation of R1
c     u is the face-centered advection velocity
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vc_update_density2d(
     &     dx,dt,a0,a1,a2,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nR0gc0,nR0gc1,
     &     R0,
     &     nR1gc0,nR1gc1,
     &     R1,
     &     nugc0,nugc1,
     &     u0,u1,
     &     nrhalfgc0,nrhalfgc1,
     &     rhalf0,rhalf1,
     &     nSgc0,nSgc1,
     &     S,
     &     nRgc0,nRgc1,
     &     R)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nR0gc0,nR0gc1
      INTEGER nR1gc0,nR1gc1
      INTEGER nugc0,nugc1
      INTEGER nrhalfgc0,nrhalfgc1
      INTEGER nRgc0,nRgc1
      INTEGER nSgc0,nSgc1

      REAL dx(0:NDIM-1),dt,a0,a1,a2

      REAL R0(CELL2dVECG(ifirst,ilast,nR0gc))
      REAL R1(CELL2dVECG(ifirst,ilast,nR1gc))
      REAL S(CELL2dVECG(ifirst,ilast,nSgc))
      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))
      REAL rhalf0(FACE2d0VECG(ifirst,ilast,nrhalfgc))
      REAL rhalf1(FACE2d1VECG(ifirst,ilast,nrhalfgc))
c
c     Input/Output.
c
      REAL R(CELL2dVECG(ifirst,ilast,nRgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL Px0,Px1
c
c     Compute R = a0*R0 + a1*R1 - a2*dt*(div[r_fc*u]+S).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            Px0 = (rhalf0(ic0+1,ic1)*u0(ic0+1,ic1) -
     &             rhalf0(ic0,ic1)*u0(ic0,ic1))/dx(0)
            R(ic0,ic1) = a0*R0(ic0,ic1) + a1*R1(ic0,ic1)
     &                   + a2 * dt * (-Px0 + S(ic0,ic1))
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            Px1 = (rhalf1(ic1+1,ic0)*u1(ic1+1,ic0) -
     &             rhalf1(ic1,ic0)*u1(ic1,ic0))/dx(1)
            R(ic0,ic1) = R(ic0,ic1) - a2* dt * Px1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Convert a side-centered vector field into a face-centered vector
c     field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_side_to_face2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     u_sc0,u_sc1,u_sc_gcw,
     &     u_fc0,u_fc1,u_fc_gcw)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER u_sc_gcw,u_fc_gcw

      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))
c
c     Input/Output.
c
      REAL u_fc0(FACE2d0(ifirst,ilast,u_fc_gcw))
      REAL u_fc1(FACE2d1(ifirst,ilast,u_fc_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gcw
c
c     Swap the velocity values.
c
      gcw = min(u_sc_gcw,u_fc_gcw)

      do i1 = ifirst1-gcw,ilast1+gcw
         do i0 = ifirst0-gcw,ilast0+gcw+1
            u_fc0(i0,i1) = u_sc0(i0,i1)
         enddo
      enddo

      do i1 = ifirst1-gcw,ilast1+gcw+1
         do i0 = ifirst0-gcw,ilast0+gcw
            u_fc1(i1,i0) = u_sc1(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Convert a face-centered vector field into a side-centered vector
c     field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_face_to_side2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     u_fc0,u_fc1,u_fc_gcw,
     &     u_sc0,u_sc1,u_sc_gcw)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER u_fc_gcw,u_sc_gcw

      REAL u_fc0(FACE2d0(ifirst,ilast,u_fc_gcw))
      REAL u_fc1(FACE2d1(ifirst,ilast,u_fc_gcw))
c
c     Input/Output.
c
      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER gcw
c
c     Swap the velocity values.
c
      gcw = min(u_sc_gcw,u_fc_gcw)

      do i1 = ifirst1-gcw,ilast1+gcw
         do i0 = ifirst0-gcw,ilast0+gcw+1
            u_sc0(i0,i1) = u_fc0(i0,i1)
         enddo
      enddo

      do i1 = ifirst1-gcw,ilast1+gcw+1
         do i0 = ifirst0-gcw,ilast0+gcw
            u_sc1(i0,i1) = u_fc1(i1,i0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy SAMRAI velocity and pressure patch data to PETSc Vec
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_to_patchlevel_vec_mac2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     first_local,last_local,
     &     p_cc,p_cc_gcw,
     &     u_sc0,u_sc1,u_sc_gcw,
     &     p_dof_cc,p_dof_cc_gcw,
     &     u_dof_sc0, u_dof_sc1,u_dof_sc_gcw,
     &     arr)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER first_local,last_local
      INTEGER p_cc_gcw,u_sc_gcw
      INTEGER p_dof_cc_gcw,u_dof_sc_gcw
c
c     Input/Output.
c
      REAL p_cc(CELL2d(ifirst,ilast,p_cc_gcw))
      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))

      INTEGER p_dof_cc(CELL2d(ifirst,ilast,p_dof_cc_gcw))
      INTEGER u_dof_sc0(SIDE2d0(ifirst,ilast,u_dof_sc_gcw))
      INTEGER u_dof_sc1(SIDE2d1(ifirst,ilast,u_dof_sc_gcw))

      REAL arr(first_local:last_local-1)
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER dof_idx
c
c     Copy velocity to array.
c

      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0+1
            dof_idx = u_dof_sc0(i0,i1)
            if ((dof_idx .ge. first_local) .and.
     &           (dof_idx .lt. last_local)) then
               arr(dof_idx) = u_sc0(i0,i1)
            end if
         enddo
      enddo

      do i1 = ifirst1,ilast1+1
         do i0 = ifirst0,ilast0
            dof_idx = u_dof_sc1(i0,i1)
            if ((dof_idx .ge. first_local) .and.
     &           (dof_idx .lt. last_local)) then
               arr(dof_idx) = u_sc1(i0,i1)
            end if
         enddo
      enddo
c

c
c     Copy pressure to array
c

      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            dof_idx = p_dof_cc(i0,i1)
            if ((dof_idx .ge. first_local) .and.
     &           (dof_idx .lt. last_local)) then
               arr(dof_idx) = p_cc(i0,i1)
            end if
         enddo
      enddo

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy PETSc Vec velocity and pressure to SAMRAI patch data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_from_patchlevel_vec_mac2d(
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     first_local,last_local,
     &     p_cc,p_cc_gcw,
     &     u_sc0,u_sc1,u_sc_gcw,
     &     p_dof_cc,p_dof_cc_gcw,
     &     u_dof_sc0, u_dof_sc1,u_dof_sc_gcw,
     &     arr)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0
      INTEGER ifirst1,ilast1

      INTEGER first_local,last_local
      INTEGER p_cc_gcw,u_sc_gcw
      INTEGER p_dof_cc_gcw,u_dof_sc_gcw
c
c     Input/Output.
c
      REAL p_cc(CELL2d(ifirst,ilast,p_cc_gcw))
      REAL u_sc0(SIDE2d0(ifirst,ilast,u_sc_gcw))
      REAL u_sc1(SIDE2d1(ifirst,ilast,u_sc_gcw))

      INTEGER p_dof_cc(CELL2d(ifirst,ilast,p_dof_cc_gcw))
      INTEGER u_dof_sc0(SIDE2d0(ifirst,ilast,u_dof_sc_gcw))
      INTEGER u_dof_sc1(SIDE2d1(ifirst,ilast,u_dof_sc_gcw))

      REAL arr(first_local:last_local-1)
c
c     Local variables.
c
      INTEGER i0,i1
      INTEGER dof_idx
c
c     Copy array velocity  to SAMRAI patch data.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0+1
            dof_idx = u_dof_sc0(i0,i1)
            if ((dof_idx .ge. first_local) .and.
     &           (dof_idx .lt. last_local)) then
               u_sc0(i0,i1) = arr(dof_idx)
            end if
         enddo
      enddo

      do i1 = ifirst1,ilast1+1
         do i0 = ifirst0,ilast0
            dof_idx = u_dof_sc1(i0,i1)
            if ((dof_idx .ge. first_local) .and.
     &           (dof_idx .lt. last_local)) then
               u_sc1(i0,i1) =  arr(dof_idx)
            end if
         enddo
      enddo
c

c
c     Copy array pressure to SAMRAI patch data.
c

      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            dof_idx = p_dof_cc(i0,i1)
            if ((dof_idx .ge. first_local) .and.
     &           (dof_idx .lt. last_local)) then
               p_cc(i0,i1) = arr(dof_idx)
            end if
         enddo
      enddo

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
