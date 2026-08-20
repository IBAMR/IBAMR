c ---------------------------------------------------------------------
c dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,2)dnl
define(real,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Reference Map Technique
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c------------------------------ div_t_2d ---------------------------------
      subroutine div_tmain2d(
     &     f11, f12, f21, f22, f_gcw, GGS, nu_s,
     &     Txx_s, Txy_s, Tyx_s, Tyy_s, ts_gcw,
     &     phi, phi_gcw, H_phi, H_phi_gcw,
     &     xi0, xi1, xi_gcw,
     &     ilower0, iupper0, ilower1, iupper1,
     &     dx, div_t_gcw,
     &     div_t1, div_t2)
 
      implicit none

c-----------------------------------------------------------------------
c  arguments
c-----------------------------------------------------------------------
      integer ilower0, iupper0
      integer ilower1, iupper1
      integer f_gcw, ts_gcw, div_t_gcw
      integer phi_gcw, H_phi_gcw, xi_gcw
      real    dx(0:1)
      real    GGS
      real    nu_s

      real  f11(CELL2d(ilower,iupper,f_gcw))
      real  f12(CELL2d(ilower,iupper,f_gcw))
      real  f21(CELL2d(ilower,iupper,f_gcw))
      real  f22(CELL2d(ilower,iupper,f_gcw))

      real  Txx_s(CELL2d(ilower,iupper,ts_gcw))
      real  Txy_s(CELL2d(ilower,iupper,ts_gcw))
      real  Tyx_s(CELL2d(ilower,iupper,ts_gcw))
      real  Tyy_s(CELL2d(ilower,iupper,ts_gcw))

      real  phi(CELL2d(ilower,iupper,phi_gcw))
      real  H_phi(CELL2d(ilower,iupper,H_phi_gcw))
      real  xi0(CELL2d(ilower,iupper,xi_gcw))
      real  xi1(CELL2d(ilower,iupper,xi_gcw))

      real  div_t1(SIDE2d0(ilower,iupper,div_t_gcw))
      real  div_t2(SIDE2d1(ilower,iupper,div_t_gcw))

c-----------------------------------------------------------------------
c  locals
c-----------------------------------------------------------------------
      integer i0,i1
      integer imax1_0, imax1_1, imax2_0, imax2_1
      real fac0,fac1
      double precision a, max1, max2
      double precision thr
      logical valid

      double precision sxx_l, sxx_r
      double precision syx_b, syx_t
      double precision syy_b, syy_t
      double precision sxy_l, sxy_r
c-----------------------------------------------------------------------
c  setup
c-----------------------------------------------------------------------
      fac0 = 1.d0/dx(0)
      fac1 = 1.d0/dx(1)

c     simple threshold for reporting
      thr = 1.d-2

c-----------------------------------------------------------------------
c  compute solid stress tensor
c-----------------------------------------------------------------------
      call solid_stress_tensor_2d(Txx_s,Txy_s,Tyx_s,Tyy_s,
     &     ilower0,iupper0,ilower1,iupper1,ts_gcw,
     &     phi,phi_gcw,xi0,xi1,xi_gcw,H_phi,H_phi_gcw,dx,
     &     f11,f12,f21,f22,f_gcw,GGS, nu_s)

c-----------------------------------------------------------------------
c  zero outputs first
c-----------------------------------------------------------------------
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0+1
            div_t1(i0,i1) = 0.d0
         enddo
      enddo

      do i1 = ilower1, iupper1+1
         do i0 = ilower0, iupper0
            div_t2(i0,i1) = 0.d0
         enddo
      enddo
 
c-----------------------------------------------------------------------
c  div_t1 at x-faces
c-----------------------------------------------------------------------
      do i1 = ilower1+1, iupper1-1
         do i0 = ilower0+1, iupper0

c            do i1 = ilower1, iupper1
c         do i0 = ilower0, iupper0+1

            sxx_r = H_phi(i0  ,i1) * Txx_s(i0  ,i1)
            sxx_l = H_phi(i0-1,i1) * Txx_s(i0-1,i1)

            syx_t = 0.25d0 * (
     &          H_phi(i0  ,i1+1) * Tyx_s(i0  ,i1+1)
     &        + H_phi(i0-1,i1+1) * Tyx_s(i0-1,i1+1) )

            syx_b = 0.25d0 * (
     &          H_phi(i0  ,i1-1) * Tyx_s(i0  ,i1-1)
     &        + H_phi(i0-1,i1-1) * Tyx_s(i0-1,i1-1) )

            div_t1(i0,i1) = 1.0 *((sxx_r - sxx_l) * fac0
     &                     + (syx_t - syx_b) * fac1)

c      if (i1 .ge. iupper1-3) then
c         div_t1(i0,i1) = 0.d0
c         endif

         enddo
      enddo

c-----------------------------------------------------------------------
c  div_t2 at y-faces
c-----------------------------------------------------------------------
      do i1 = ilower1+1, iupper1
         do i0 = ilower0+1, iupper0-1

c          do i1 = ilower1, iupper1+1
c         do i0 = ilower0, iupper0

            syy_t = H_phi(i0,i1  ) * Tyy_s(i0,i1  )
            syy_b = H_phi(i0,i1-1) * Tyy_s(i0,i1-1)

            sxy_r = 0.25d0 * (
     &          H_phi(i0+1,i1  ) * Txy_s(i0+1,i1  )
     &        + H_phi(i0+1,i1-1) * Txy_s(i0+1,i1-1) )

            sxy_l = 0.25d0 * (
     &          H_phi(i0-1,i1  ) * Txy_s(i0-1,i1  )
     &        + H_phi(i0-1,i1-1) * Txy_s(i0-1,i1-1) )

            div_t2(i0,i1) = 1.0 *( (syy_t - syy_b) * fac1
     &                     + (sxy_r - sxy_l) * fac0)


c      if (i1 .ge. iupper1-3) then
c         div_t2(i0,i1) = 0.d0
c         endif

         enddo
      enddo

c-----------------------------------------------------------------------
c  compact debug summary
c-----------------------------------------------------------------------
      max1 = 0.d0
      imax1_0 = ilower0
      imax1_1 = ilower1
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0+1
            a = dabs(dble(div_t1(i0,i1)))
            if (a .gt. max1) then
               max1 = a
               imax1_0 = i0
               imax1_1 = i1
            endif
         enddo
      enddo

      max2 = 0.d0
      imax2_0 = ilower0
      imax2_1 = ilower1
      do i1 = ilower1, iupper1+1
         do i0 = ilower0, iupper0
            a = dabs(dble(div_t2(i0,i1)))
            if (a .gt. max2) then
               max2 = a
               imax2_0 = i0
               imax2_1 = i1
            endif
         enddo
      enddo

c      print *, 'DBG divT: patch=', ilower0, iupper0, ilower1, iupper1
c      print *, 'DBG divT: max|div_t1|=', max1, ' at ', imax1_0, imax1_1
c      print *, 'DBG divT: max|div_t2|=', max2, ' at ', imax2_0, imax2_1

c      if (max1 .gt. thr) then
c         print *, 'DBG divT1 spike: phi/H = ',
c     &        phi(imax1_0,imax1_1), H_phi(imax1_0,imax1_1)
c      endif

c      if (max2 .gt. thr) then
c         print *, 'DBG divT2 spike: phi/H = ',
c     &        phi(imax2_0,imax2_1), H_phi(imax2_0,imax2_1)
c      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-------------------- solid_stress_tensor_2d --------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solid_stress_tensor_2d(Txx_s, Txy_s, Tyx_s, Tyy_s,
     &     ilower0, iupper0, ilower1, iupper1, ts_gcw, 
     &     phi, phi_gcw, xi0, xi1, xi_gcw, H_phi, H_phi_gcw, dx,
     &     f11, f12, f21, f22, f_gcw, GGS, nu_s)

      implicit none

      integer ilower0, iupper0
      integer ilower1, iupper1
      integer i0, i1
      integer ts_gcw, phi_gcw, xi_gcw, H_phi_gcw, f_gcw

      real dx(0:NDIM-1)
      real GGS
      real nu_s

c     Inputs
      real phi(CELL2d(ilower,iupper,phi_gcw))
      real xi0(CELL2d(ilower,iupper,xi_gcw))
      real xi1(CELL2d(ilower,iupper,xi_gcw))

c     Outputs
      real Txx_s(CELL2d(ilower,iupper,ts_gcw))
      real Txy_s(CELL2d(ilower,iupper,ts_gcw))
      real Tyx_s(CELL2d(ilower,iupper,ts_gcw))
      real Tyy_s(CELL2d(ilower,iupper,ts_gcw))

c     Local arrays
      real H_phi(CELL2d(ilower,iupper,H_phi_gcw))
      real f11(CELL2d(ilower,iupper,f_gcw))
      real f12(CELL2d(ilower,iupper,f_gcw))
      real f21(CELL2d(ilower,iupper,f_gcw))
      real f22(CELL2d(ilower,iupper,f_gcw))

      real nu, lam
      real E11, E22, E12, trE
      double precision eps, margin
      double precision fac0, fac1
      double precision g1, g2, g3, g4

      double precision txx, txy, tyx, tyy
      integer n_T_nan, n_T_huge, n_T_finite
      double precision T_absmax
      logical Tok

c     --- parameters ---
      nu  = nu_s
      lam = 2.d0*GGS*nu/(1.d0 - 2.d0*nu)

      eps    = 1.5d0 * dmax1(dble(dx(0)), dble(dx(1)))
      margin = 0.5d0 * eps

      fac0 = 1.d0/dble(dx(0))
      fac1 = 1.d0/dble(dx(1))

      n_T_nan    = 0
      n_T_huge   = 0
      n_T_finite = 0
      T_absmax   = 0.0d0

c ---- ZERO stress on the whole patch interior ----
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
            Txx_s(i0,i1) = 0.d0
            Txy_s(i0,i1) = 0.d0
            Tyx_s(i0,i1) = 0.d0
            Tyy_s(i0,i1) = 0.d0
         enddo
      enddo

c     Build H_phi
      call smoothed_heaviside_2d(phi, H_phi, dx,
     &     ilower0,iupper0, ilower1, iupper1, phi_gcw, H_phi_gcw)

c     Compute F (inverse grad xi)
      call deformation_gradient_tensor_F_2d(xi0, xi1,
     &     phi, phi_gcw,
     &     f11, f12, f21, f22, xi_gcw, f_gcw,
     &     dx, ilower0,iupper0, ilower1,iupper1)

c     ---- linear elastic stress using U = F - I (only where F is valid) ----
      do i1 = ilower1+1, iupper1-1
         do i0 = ilower0+1, iupper0-1


c            if (i1 .ge. iupper1-2) cycle

c            if (dble(phi(i0,i1))   .gt. -margin) cycle
c            if (dble(phi(i0+1,i1)) .gt. -margin) cycle
c            if (dble(phi(i0-1,i1)) .gt. -margin) cycle
c            if (dble(phi(i0,i1+1)) .gt. -margin) cycle
c            if (dble(phi(i0,i1-1)) .gt. -margin) cycle

c            E11 = f11(i0,i1) - 1.d0
c            E22 = f22(i0,i1) - 1.d0
c            E12 = 0.5d0 * (f12(i0,i1) + f21(i0,i1))
c            trE = E11 + E22

c            Txx_s(i0,i1) = 2.d0*GGS*E11 + lam*trE
c            Tyy_s(i0,i1) = 2.d0*GGS*E22 + lam*trE
c            Txy_s(i0,i1) = 2.d0*GGS*E12
c            Tyx_s(i0,i1) = Txy_s(i0,i1)


c      Txx_s(i0,i1) =
c     & 2.d0*GGS*
c     & ( f11(i0,i1)*f11(i0,i1)
c     & + f12(i0,i1)*f12(i0,i1) -1.d0 )

c      Tyy_s(i0,i1) =
c     & 2.d0*GGS*
c     & ( f21(i0,i1)*f21(i0,i1)
c     & + f22(i0,i1)*f22(i0,i1) - 1.d0 )

c      Txy_s(i0,i1) =
c     & 2.d0*GGS*
c     & ( f11(i0,i1)*f21(i0,i1)
c     & + f12(i0,i1)*f22(i0,i1) )

c      Tyx_s(i0,i1) = Txy_s(i0,i1)


c     grad(xi) directly, as in Appendix C
      g1 = 0.5d0*(dble(xi0(i0+1,i1)) - dble(xi0(i0-1,i1)))*fac0
      g2 = 0.5d0*(dble(xi0(i0,i1+1)) - dble(xi0(i0,i1-1)))*fac1
      g3 = 0.5d0*(dble(xi1(i0+1,i1)) - dble(xi1(i0-1,i1)))*fac0
      g4 = 0.5d0*(dble(xi1(i0,i1+1)) - dble(xi1(i0,i1-1)))*fac1


c     Paper Appendix B/C form, using alpha_x, alpha_y, beta_x, beta_y
c     g1 = alpha_x, g2 = alpha_y, g3 = beta_x, g4 = beta_y

      Txx_s(i0,i1) =
     & 2.d0*GGS*
     & ( g2*g2 + g4*g4 - 1.d0 )

      Txy_s(i0,i1) =
     & -2.d0*GGS*
     & ( g1*g2 + g3*g4 )

      Tyx_s(i0,i1) = Txy_s(i0,i1)

      Tyy_s(i0,i1) =
     & 2.d0*GGS*
     & ( g1*g1 + g3*g3 - 1.d0 )



c           --- minimal stress stats ---
            txx = Txx_s(i0,i1)
            txy = Txy_s(i0,i1)
            tyx = Tyx_s(i0,i1)
            tyy = Tyy_s(i0,i1)

            Tok = .true.
            if (txx .ne. txx) Tok = .false.
            if (txy .ne. txy) Tok = .false.
            if (tyx .ne. tyx) Tok = .false.
            if (tyy .ne. tyy) Tok = .false.

            if (dabs(txx) .ge. 1.0d90) Tok = .false.
            if (dabs(txy) .ge. 1.0d90) Tok = .false.
            if (dabs(tyx) .ge. 1.0d90) Tok = .false.
            if (dabs(tyy) .ge. 1.0d90) Tok = .false.

            if (Tok) then
               n_T_finite = n_T_finite + 1
               T_absmax = max(T_absmax, dabs(txx))
               T_absmax = max(T_absmax, dabs(txy))
               T_absmax = max(T_absmax, dabs(tyx))
               T_absmax = max(T_absmax, dabs(tyy))
            else
               if ( (txx .ne. txx) .or. (txy .ne. txy) .or.
     &              (tyx .ne. tyx) .or. (tyy .ne. tyy) ) then
                  n_T_nan = n_T_nan + 1
               else
                  n_T_huge = n_T_huge + 1
               endif
            endif

         enddo
      enddo

c ---- print once per call ----
c      write(*,*) 'DBG T: n_T_nan   =', n_T_nan
c      write(*,*) 'DBG T: n_T_huge  =', n_T_huge
c      write(*,*) 'DBG T: n_T_finite=', n_T_finite
c      write(*,*) 'DBG T: T_absmax  =', T_absmax

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-------------------------- smoothed_heaviside_2d -------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            subroutine smoothed_heaviside_2d(phi, H_phi, dx, 
     &  ilower0,iupper0, ilower1, iupper1, phi_gcw , H_phi_gcw)

      implicit none
      integer ilower0, iupper0, ilower1, iupper1
      integer phi_gcw, H_phi_gcw
      real phi(CELL2d(ilower,iupper,phi_gcw))
      real H_phi(CELL2d(ilower,iupper,H_phi_gcw))
      integer :: i0, i1
      real :: epsilon, ph
      real, parameter :: pi = 3.14159265358979323846
      real dx(0:NDIM-1)
      real x
      
      epsilon = 1.5d0 * max(dx(0), dx(1))

      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0

            ph = phi(i0,i1)

            if (ph <= -epsilon) then
               H_phi(i0,i1) = 1.d0
            else if (ph >=  epsilon) then
               H_phi(i0,i1) = 0.d0
            else
              H_phi(i0,i1) = 0.5d0 * ( 1.0d0 - ph/epsilon -
     &            1/pi * sin(pi* ph/epsilon) )
      endif


         end do
      end do


      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-------------------- deformation_gradient_tensor_F_2d --------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deformation_gradient_tensor_F_2d(xi0, xi1,
     &     phi, phi_gcw,
     &     f11, f12, f21, f22, xi_gcw, f_gcw,
     &     dx, ilower0,iupper0, ilower1, iupper1)

      implicit none

      integer phi_gcw
      integer xi_gcw, f_gcw
      integer ilower0, iupper0, ilower1, iupper1
      integer i0, i1

      real phi(CELL2d(ilower,iupper,phi_gcw))
      real xi0(CELL2d(ilower,iupper,xi_gcw))
      real xi1(CELL2d(ilower,iupper,xi_gcw))

      real f11(CELL2d(ilower,iupper,f_gcw))
      real f12(CELL2d(ilower,iupper,f_gcw))
      real f21(CELL2d(ilower,iupper,f_gcw))
      real f22(CELL2d(ilower,iupper,f_gcw))

      real dx(0:NDIM-1)

c     DP scalars for det math
      double precision fac0, fac1
      double precision g0, g1, g2, g3, g4, det
      double precision eps, margin

      integer nfail, nint
      integer n_det_nan, n_det_finite, n_det_huge
      double precision det_min, det_max
      logical det_ok 

c-----------------------------------------------------------------------
c  setup
c-----------------------------------------------------------------------
      fac0 = 1.0d0/dble(dx(0))
      fac1 = 1.0d0/dble(dx(1))

      eps    = 1.5d0 * dmax1(dble(dx(0)), dble(dx(1)))
      margin = 0.5d0 * eps

      nfail = 0
      nint  = 0

      n_det_nan    = 0
      n_det_finite = 0
      n_det_huge   = 0
      det_min      = 1.0d99
      det_max      = -1.0d99

c-----------------------------------------------------------------------
c  initialize F on patch interior
c-----------------------------------------------------------------------
      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
            f11(i0,i1) = 1.d0
            f12(i0,i1) = 0.d0
            f21(i0,i1) = 0.d0
            f22(i0,i1) = 1.d0
         enddo
      enddo

c-----------------------------------------------------------------------
c  interior loop
c-----------------------------------------------------------------------
      do i1 = ilower1+1, iupper1-1
         do i0 = ilower0+1, iupper0-1

c     require stencil neighbors also solid
c            if (dble(phi(i0  ,i1  )) .gt. -margin) goto 1000
c            if (dble(phi(i0+1,i1  )) .gt. -margin) goto 1000
c            if (dble(phi(i0-1,i1  )) .gt. -margin) goto 1000
c            if (dble(phi(i0  ,i1+1)) .gt. -margin) goto 1000
c            if (dble(phi(i0  ,i1-1)) .gt. -margin) goto 1000

c     grad(xi)
            g1 = 0.5d0*(dble(xi0(i0+1,i1)) - dble(xi0(i0-1,i1)))*fac0
            g2 = 0.5d0*(dble(xi0(i0,i1+1)) - dble(xi0(i0,i1-1)))*fac1
            g3 = 0.5d0*(dble(xi1(i0+1,i1)) - dble(xi1(i0-1,i1)))*fac0
            g4 = 0.5d0*(dble(xi1(i0,i1+1)) - dble(xi1(i0,i1-1)))*fac1

c     det = det(grad xi)
            det = g1*g4 - g2*g3

c     det stats
            if (det .ne. det) then
               n_det_nan = n_det_nan + 1
            else if (dabs(det) .ge. 1.0d90) then
               n_det_huge = n_det_huge + 1
            else
               n_det_finite = n_det_finite + 1
               if (det .lt. det_min) det_min = det
               if (det .gt. det_max) det_max = det
            endif

c     usable det?
            det_ok = .true.
            if (det .ne. det) det_ok = .false.
            if (dabs(det) .ge. 1.0d90) det_ok = .false.
            if (dabs(det) .lt. 1.0d-14) det_ok = .false.

            if (.not. det_ok) then
               nfail = nfail + 1
               goto 1000
            endif

c     F = inv(grad xi)
            g0 = 1.d0/det
            f11(i0,i1) = dble(g4*g0)
            f12(i0,i1) = dble(-g2*g0)
            f21(i0,i1) = dble(-g3*g0)
            f22(i0,i1) = dble(g1*g0)



c            if (i0 .eq. (ilower0+iupper0)/2 .and.
c     &          i1 .eq. (ilower1+iupper1)/2) then
c               write(*,*) 'DBG XI CENTER i0 i1=', i0, i1
c               write(*,*) 'DBG gradxi g1 g2 g3 g4=', g1, g2, g3, g4
c               write(*,*) 'DBG F f11 f12 f21 f22=',
c     &              f11(i0,i1), f12(i0,i1), f21(i0,i1), f22(i0,i1)
c               write(*,*) 'DBG det_gradxi=', det
c            endif



            nint = nint + 1
            goto 2000

 1000       continue
c           leave identity already set

 2000       continue
         enddo
      enddo

c-----------------------------------------------------------------------
c  print once per call
c-----------------------------------------------------------------------
c      write(*,*) 'DBG DET: n_det_nan   =', n_det_nan
c      write(*,*) 'DBG DET: n_det_huge  =', n_det_huge
c      write(*,*) 'DBG DET: n_det_finite=', n_det_finite
c      if (n_det_finite .gt. 0) then
c         write(*,*) 'DBG det_min=', det_min, ' det_max=', det_max
c      endif
c      write(*,*) 'DBG nfail=', nfail, ' nint=', nint


     
            return
      end
