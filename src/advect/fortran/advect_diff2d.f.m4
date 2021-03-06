c ---------------------------------------------------------------------
c
c Copyright (c) 2006 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the advective flux corresponding to a face centered value
c     and a face centered advective velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_flux2d(
     &     dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     nfluxgc0,nfluxgc1,
     &     u0,u1,
     &     qhalf0,qhalf1,
     &     flux0,flux1)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nugc0,nugc1
      INTEGER nqhalfgc0,nqhalfgc1
      INTEGER nfluxgc0,nfluxgc1

      REAL dt

      REAL u0(FACE2d0VECG(ifirst,ilast,nugc))
      REAL u1(FACE2d1VECG(ifirst,ilast,nugc))

      REAL qhalf0(FACE2d0VECG(ifirst,ilast,nqhalfgc))
      REAL qhalf1(FACE2d1VECG(ifirst,ilast,nqhalfgc))
c
c     Input/Output.
c
      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Compute the time integral of the advective flux.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0-1,ilast0
            flux0(ic0+1,ic1) = dt*u0(ic0+1,ic1)*qhalf0(ic0+1,ic1)
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1-1,ilast1
            flux1(ic1+1,ic0) = dt*u1(ic1+1,ic0)*qhalf1(ic1+1,ic0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_consdiff2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qval)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1)

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
c
c     Update a quantity using flux differencing.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            qval(ic0,ic1) = qval(ic0,ic1)
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing but include the proper
c     source term to account for a non-discretely divergence free
c     advection velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_consdiffwithdivsource2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqfluxgc0,nqfluxgc1,
     &     nufluxgc0,nufluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qflux0,qflux1,
     &     uflux0,uflux1,
     &     qval)
c
      implicit none
      REAL fourth
      parameter (fourth=0.25d0)
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nfluxgc0,nfluxgc1
      INTEGER nqfluxgc0,nqfluxgc1
      INTEGER nufluxgc0,nufluxgc1
      INTEGER nqvalgc0,nqvalgc1

      REAL dx(0:NDIM-1),dt

      REAL flux0(FACE2d0VECG(ifirst,ilast,nfluxgc))
      REAL flux1(FACE2d1VECG(ifirst,ilast,nfluxgc))

      REAL qflux0(FACE2d0VECG(ifirst,ilast,nqfluxgc))
      REAL qflux1(FACE2d1VECG(ifirst,ilast,nqfluxgc))

      REAL uflux0(FACE2d0VECG(ifirst,ilast,nufluxgc))
      REAL uflux1(FACE2d1VECG(ifirst,ilast,nufluxgc))
c
c     Input/Output.
c
      REAL qval(CELL2dVECG(ifirst,ilast,nqvalgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL divsource
c
c     Update a quantity using flux differencing.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            divsource = (fourth/dt)*
     &           ( qflux0(ic0+1,ic1) + qflux0(ic0,ic1)
     &           + qflux1(ic1+1,ic0) + qflux1(ic1,ic0) )*
     &           ( (uflux0(ic0+1,ic1)-uflux0(ic0,ic1))/dx(0)
     &           + (uflux1(ic1+1,ic0)-uflux1(ic1,ic0))/dx(1) )

            qval(ic0,ic1) = qval(ic0,ic1) + divsource
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the advective derivative N = [u_ADV*grad(q)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_derivative2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nuadvgc0,nuadvgc1,
     &     nqgc0,nqgc1,
     &     uadv0,uadv1,
     &     q0,q1,
     &     nNgc0,nNgc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nuadvgc0,nuadvgc1
      INTEGER nqgc0,nqgc1
      INTEGER nNgc0,nNgc1

      REAL dx(0:NDIM-1)

      REAL uadv0(FACE2d0VECG(ifirst,ilast,nuadvgc))
      REAL uadv1(FACE2d1VECG(ifirst,ilast,nuadvgc))

      REAL q0(FACE2d0VECG(ifirst,ilast,nqgc))
      REAL q1(FACE2d1VECG(ifirst,ilast,nqgc))
c
c     Input/Output.
c
      REAL N(CELL2dVECG(ifirst,ilast,nNgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL U,V
      REAL Qx0,Qx1
c
c     Compute (U,V)*grad(q).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            U = 0.5d0*(uadv0(ic0+1,ic1)+uadv0(ic0,ic1))
            Qx0 = (q0(ic0+1,ic1)-q0(ic0,ic1))/dx(0)
            N(ic0,ic1) = U*Qx0
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            V = 0.5d0*(uadv1(ic1+1,ic0)+uadv1(ic1,ic0))
            Qx1 = (q1(ic1+1,ic0)-q1(ic1,ic0))/dx(1)
            N(ic0,ic1) = N(ic0,ic1) + V*Qx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the convective derivative N = div[q*u_ADV].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine convect_derivative2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nuadvgc0,nuadvgc1,
     &     nqgc0,nqgc1,
     &     uadv0,uadv1,
     &     q0,q1,
     &     nNgc0,nNgc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nuadvgc0,nuadvgc1
      INTEGER nqgc0,nqgc1
      INTEGER nNgc0,nNgc1

      REAL dx(0:NDIM-1)

      REAL uadv0(FACE2d0VECG(ifirst,ilast,nuadvgc))
      REAL uadv1(FACE2d1VECG(ifirst,ilast,nuadvgc))

      REAL q0(FACE2d0VECG(ifirst,ilast,nqgc))
      REAL q1(FACE2d1VECG(ifirst,ilast,nqgc))
c
c     Input/Output.
c
      REAL N(CELL2dVECG(ifirst,ilast,nNgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL QUx0,QVx1
c
c     Compute div[q*(U,V)].
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            QUx0 = (uadv0(ic0+1,ic1)*q0(ic0+1,ic1)-
     &           uadv0(ic0,ic1)*q0(ic0,ic1))/dx(0)
            N(ic0,ic1) = QUx0
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            QVx1 = (uadv1(ic1+1,ic0)*q1(ic1+1,ic0)-
     &           uadv1(ic1,ic0)*q1(ic1,ic0))/dx(1)
            N(ic0,ic1) = N(ic0,ic1) + QVx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the skew-symmetric derivative N = 0.5([u_ADV*grad(q)] +
c     div[q*u_ADV]).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine skew_sym_derivative2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nuadvgc0,nuadvgc1,
     &     nqgc0,nqgc1,
     &     uadv0,uadv1,
     &     q0,q1,
     &     nNgc0,nNgc1,
     &     N)
c
      implicit none
c
c     Input.
c
      INTEGER ifirst0,ilast0,ifirst1,ilast1

      INTEGER nuadvgc0,nuadvgc1
      INTEGER nqgc0,nqgc1
      INTEGER nNgc0,nNgc1

      REAL dx(0:NDIM-1)

      REAL uadv0(FACE2d0VECG(ifirst,ilast,nuadvgc))
      REAL uadv1(FACE2d1VECG(ifirst,ilast,nuadvgc))

      REAL q0(FACE2d0VECG(ifirst,ilast,nqgc))
      REAL q1(FACE2d1VECG(ifirst,ilast,nqgc))
c
c     Input/Output.
c
      REAL N(CELL2dVECG(ifirst,ilast,nNgc))
c
c     Local variables.
c
      INTEGER ic0,ic1
      REAL U,V
      REAL Qx0,Qx1
      REAL QUx0,QVx1
c
c     Compute 0.5*((U,V)*grad(q) + div[q*(U,V)]).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            U = 0.5d0*(uadv0(ic0+1,ic1)+uadv0(ic0,ic1))
            Qx0 = (q0(ic0+1,ic1)-q0(ic0,ic1))/dx(0)
            QUx0 = (uadv0(ic0+1,ic1)*q0(ic0+1,ic1)-
     &           uadv0(ic0,ic1)*q0(ic0,ic1))/dx(0)
            N(ic0,ic1) = 0.5d0*(U*Qx0+QUx0)
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            V = 0.5d0*(uadv1(ic1+1,ic0)+uadv1(ic1,ic0))
            Qx1 = (q1(ic1+1,ic0)-q1(ic1,ic0))/dx(1)
            QVx1 = (uadv1(ic1+1,ic0)*q1(ic1+1,ic0)-
     &           uadv1(ic1,ic0)*q1(ic1,ic0))/dx(1)
            N(ic0,ic1) = N(ic0,ic1) + 0.5d0*(V*Qx1+QVx1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
