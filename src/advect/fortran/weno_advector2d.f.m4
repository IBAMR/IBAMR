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
c     NOTE: The following employs WENO-RF: WENO-Roe is used away from
c     sonic points and WENO-LLF is used at sonic points.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine weno_convective_fluxes2d(
     &     fface0,fface1,fface_gcw,
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
c     Parameters of the numerical method.
c
      INTEGER R
      REAL EPSILON
      parameter(R=3)
      parameter(EPSILON=1.0d-10)
c
c     Input.
c
      INTEGER fface_gcw,F_gcw,Q_gcw,U_gcw

      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1

      REAL Q(CELL2d(ilower,iupper,Q_gcw))
      REAL U(CELL2d(ilower,iupper,U_gcw),0:NDIM-1)
c
c     Input/Output.
c
      REAL fface0(FACE2d0(ilower,iupper,fface_gcw))
      REAL fface1(FACE2d1(ilower,iupper,fface_gcw))
      REAL F(CELL2d(ilower,iupper,F_gcw))
c
c     Local variables.
c
      INTEGER dim
      INTEGER ic0,ic1
      INTEGER k,l,m
      REAL a
      REAL IS,IS_sub
      REAL alpha(0:R),alpha_sum
      REAL omega(0:R)
      REAL qface(0:R)
c
c     Compute the x-face fluxes on the "upper" side of each cell
c     (ic0,ic1).
c
      dim = 0

      do ic1 = ilower1-R,iupper1+R
         do ic0 = ilower0-R,iupper0+R
            F(ic0,ic1) = U(ic0,ic1,dim)*Q(ic0,ic1)
         enddo
      enddo

      do ic1 = ilower1,iupper1
         do ic0 = ilower0-1,iupper0
            fface0(ic0+1,ic1) = 0.d0
         enddo
      enddo

      do ic1 = ilower1,iupper1
         do ic0 = ilower0-1,iupper0
            if (U(ic0,ic1,dim)*U(ic0+1,ic1,dim) .lt. 0.d0) then

               a = 0.d0
               do k = -R+1,R
                  a = max(a,abs(2.d0*U(ic0+k,ic1,dim))) ! allow for the case that Q = U
               enddo
c
c     Compute F+ = 0.5*(F(Q)+a*Q).
c
               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*0.5d0*(
     &                    F(ic0+k+l-R+1,ic1)+a*Q(ic0+k+l-R+1,ic1))
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*0.5d0*(
     &                       F(ic0+k+l-R+1,ic1)+Q(ic0+k+l-R+1,ic1))
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface0(ic0+1,ic1) =
     &                 fface0(ic0+1,ic1) + omega(k)*qface(k)
               enddo
c
c     Compute F- = 0.5*(F(Q)-a*Q).
c
               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*0.5d0*(
     &                    F(ic0-k-l+R,ic1)-a*Q(ic0-k-l+R,ic1))
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*0.5d0*(
     &                       F(ic0-k-l+R,ic1)-a*Q(ic0-k-l+R,ic1))
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface0(ic0+1,ic1) =
     &                 fface0(ic0+1,ic1) + omega(k)*qface(k)
               enddo

            else if (U(ic0,ic1,dim) .ge. 0.d0) then

               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*F(ic0+k+l-R+1,ic1)
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*F(ic0+k+l-R+1,ic1)
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
              enddo

               do k = 0,R-1
                  fface0(ic0+1,ic1) =
     &                 fface0(ic0+1,ic1) + omega(k)*qface(k)
               enddo

            else

               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*F(ic0-k-l+R,ic1)
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*F(ic0-k-l+R,ic1)
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface0(ic0+1,ic1) =
     &                 fface0(ic0+1,ic1) + omega(k)*qface(k)
               enddo

            endif
         enddo
      enddo
c
c     Compute the y-face fluxes on the "upper" side of each cell
c     (ic0,ic1).
c
      dim = 1

      do ic1 = ilower1-R,iupper1+R
         do ic0 = ilower0-R,iupper0+R
            F(ic0,ic1) = U(ic0,ic1,dim)*Q(ic0,ic1)
         enddo
      enddo

      do ic0 = ilower0,iupper0
         do ic1 = ilower1-1,iupper1
            fface1(ic1+1,ic0) = 0.d0
         enddo
      enddo

      do ic0 = ilower0,iupper0
         do ic1 = ilower1-1,iupper1
            if (U(ic0,ic1,dim)*U(ic0,ic1+1,dim) .lt. 0.d0) then

               a = 0.d0
               do k = -R+1,R
                  a = max(a,abs(2.d0*U(ic0,ic1+k,dim))) ! allow for the case that Q = V
               enddo
c
c     Compute F+ = 0.5*(F(Q)+a*Q).
c
               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*0.5d0*(
     &                    F(ic0,ic1+k+l-R+1)+a*Q(ic0,ic1+k+l-R+1))
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*0.5d0*(
     &                       F(ic0,ic1+k+l-R+1)+a*Q(ic0,ic1+k+l-R+1))
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface1(ic1+1,ic0) =
     &                 fface1(ic1+1,ic0) + omega(k)*qface(k)
               enddo
c
c     Compute F- = 0.5*(F(Q)-a*Q).
c
                 do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*0.5d0*(
     &                    F(ic0,ic1-k-l+R)-Q(ic0,ic1-k-l+R))
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*0.5d0*(
     &                       F(ic0,ic1-k-l+R)-Q(ic0,ic1-k-l+R))
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface1(ic1+1,ic0) =
     &                 fface1(ic1+1,ic0) + omega(k)*qface(k)
               enddo

            else if (U(ic0,ic1,dim) .ge. 0.d0) then

               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*F(ic0,ic1+k+l-R+1)
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*F(ic0,ic1+k+l-R+1)
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface1(ic1+1,ic0) =
     &                 fface1(ic1+1,ic0) + omega(k)*qface(k)
               enddo

            else

               do k = 0,R-1
                  qface(k) = 0.d0
                  do l = 0,R-1
                     qface(k) = qface(k) + A3(k,l)*F(ic0,ic1-k-l+R)
                  enddo
               enddo

               do k = 0,R-1
                  IS = 0.d0
                  do m = 1,R-1
                     IS_sub = 0.d0
                     do l = 0,R-1
                        IS_sub = IS_sub + D3(k,l,m)*F(ic0,ic1-k-l+R)
                     enddo
                     IS = IS + IS_sub*IS_sub
                  enddo
                  alpha(k) = C3_JS(k)/(EPSILON + IS)
               enddo

               alpha_sum = sum(alpha)
               do k = 0,R-1
                  omega(k) = alpha(k)/alpha_sum
               enddo

               do k = 0,R-1
                  fface1(ic1+1,ic0) =
     &                 fface1(ic1+1,ic0) + omega(k)*qface(k)
               enddo

            endif
         enddo
      enddo
c
      return
      end
