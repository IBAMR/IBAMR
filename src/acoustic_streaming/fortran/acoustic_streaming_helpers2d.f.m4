c ---------------------------------------------------------------------
c
c Copyright (c) 2008 - 2023 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

define(NDIM,2)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy staggered velocity field onto another staggered velocity field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_staggered_data_2d(U0, U1,
     &     U_gcw, V0, V1, V_gcw, ilower0, iupper0,
     &     ilower1, iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw, V_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1

      REAL U0(SIDE2d0(ilower,iupper,U_gcw))
      REAL U1(SIDE2d1(ilower,iupper,U_gcw))
c
c     Input/Output.
c
      REAL V0(SIDE2d0(ilower,iupper,V_gcw))
      REAL V1(SIDE2d1(ilower,iupper,V_gcw))
c
c     Local variables.
c
      INTEGER i0,i1

c
c     Copy the components of U to V.
c

      do i1 = ilower1,iupper1
        do i0 = ilower0, iupper0+1
          V0(i0,i1) = U0(i0,i1)
        enddo  
      enddo

      do i1 = ilower1,iupper1+1
        do i0 = ilower0, iupper0
          V1(i0,i1) = U1(i0,i1)
        enddo  
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Copy cell-centered field onto another cell-centered field.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine copy_cell_data_2d(U,
     &     U_gcw, V, V_gcw, ilower0, iupper0,
     &     ilower1, iupper1)
c
      implicit none
c
c     Input.
c
      INTEGER U_gcw, V_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1

      REAL U(CELL2d(ilower,iupper,U_gcw))

c
c     Input/Output.
c
      REAL V(CELL2d(ilower,iupper,V_gcw))
c
c     Local variables.
c
      INTEGER i0,i1

c
c     Copy the components of U to V.
c

      do i1 = ilower1,iupper1
        do i0 = ilower0, iupper0
          V(i0,i1) = U(i0,i1)
        enddo  
      enddo

c
      return
      end
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the forcing term for the RHS of second order momentum equation due   
c     to the first order system f1 = - < div . (rho (U1 x U1))>
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      


      subroutine acoustic_momentum_coupling_2d(U0r,U1r, 
     &     Ur_gcw, U0i, U1i, Ui_gcw, rho0, rho1, 
     &     rho_gcw, f0, f1, f_gcw, ilower0, iupper0, 
     &     ilower1, iupper1, dx)

c
      implicit none
c
c     Input.
c
      INTEGER Ur_gcw, Ui_gcw, rho_gcw, f_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1

      REAL U0r(SIDE2d0(ilower,iupper,Ur_gcw))
      REAL U1r(SIDE2d1(ilower,iupper,Ur_gcw))
      REAL U0i(SIDE2d0(ilower,iupper,Ui_gcw))
      REAL U1i(SIDE2d1(ilower,iupper,Ui_gcw))
      REAL rho0(SIDE2d0(ilower,iupper,rho_gcw))
      REAL rho1(SIDE2d1(ilower,iupper,rho_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c
      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

c     Local variables
      
      REAL fac0, fac1, g1, g2
      INTEGER i0, i1

      fac0 = 1.0d0/dx(0)
      fac1 = 1.0d0/dx(1)
      
c
c     Compute the coupling source term.
c

      do i1 = ilower1,iupper1
        do i0 = ilower0, iupper0+1


          g1 = rho0(i0+1,i1)*(U0r(i0+1,i1)**2 + U0i(i0+1,i1)**2) - 
     &         rho0(i0-1,i1)*(U0r(i0-1,i1)**2 + U0i(i0-1,i1)**2)
          g1 = 0.5d0*fac0*g1 

          g2 = 0.5d0*(rho0(i0,i1+1) + rho0(i0,i1)) *  
     &        (0.5d0*(U0r(i0,i1+1) + U0r(i0,i1)) *
     &         0.5d0*(U1r(i0,i1+1) + U1r(i0-1,i1+1)) + 
     &         0.5d0*(U0i(i0,i1+1) + U0i(i0,i1)) *
     &         0.5d0*(U1i(i0,i1+1) + U1i(i0-1,i1+1))) 
     &       - 
     &         0.5d0*(rho0(i0,i1) + rho0(i0,i1-1)) * 
     &         (0.5d0*(U0r(i0,i1) + U0r(i0,i1-1))  *
     &          0.5d0*(U1r(i0,i1) + U1r(i0-1,i1))  +  
     &          0.5d0*(U0i(i0,i1) + U0i(i0,i1-1))  *
     &          0.5d0*(U1i(i0,i1) + U1i(i0-1,i1)))
          g2 = fac1*g2;  

          f0(i0,i1) = -0.5d0*(g1+g2)

        enddo  
      enddo

      do i1 = ilower1,iupper1+1
        do i0 = ilower0, iupper0

          g1 = rho1(i0,i1+1)*(U1r(i0,i1+1)**2 + U1i(i0,i1+1)**2) -
     &         rho1(i0,i1-1)*(U1r(i0,i1-1)**2 + U1i(i0,i1-1)**2) 
          g1 = 0.5d0*fac1*g1

          g2 = 0.5d0*(rho1(i0+1,i1) + rho1(i0,i1))*
     &        (0.5d0*(U1r(i0+1,i1) + U1r(i0,i1))  *
     &         0.5d0*(U0r(i0+1,i1) + U0r(i0+1,i1-1)) + 
     &         0.5d0*(U1i(i0+1,i1) + U1i(i0,i1)) *
     &         0.5d0*(U0i(i0+1,i1) + U0i(i0+1,i1-1)))
     &       -
     &         0.5d0*(rho1(i0,i1) + rho1(i0-1,i1))*
     &        (0.5d0*(U1r(i0,i1) + U1r(i0-1,i1))  *
     &         0.5d0*(U0r(i0,i1) + U0r(i0,i1-1))  + 
     &         0.5d0*(U1i(i0,i1) + U1i(i0-1,i1))  * 
     &         0.5d0*(U0i(i0,i1) + U0i(i0,i1-1)))   

          g2 = fac0*g2

          f1(i0,i1) = -0.5d0*(g1+g2)       


        enddo  
      enddo
c     
      return
      end
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the forcing term for the RHS of second order continuity equation due   
c     to the first order system f2 =  < div . (rho1 U1)>
c     Note that we are removing the -ve sign from f2 because the second order system 
c     considers -div (rho0 U2) as the continuity equation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    


      subroutine acoustic_mass_coupling_2d(U0r, U1r,
     &     Ur_gcw, U0i, U1i, Ui_gcw, pr,
     &     pr_gcw, pi,  pi_gcw, c0,
     &     m, m_gcw, ilower0, iupper0, ilower1, iupper1, dx)

c
      implicit none
c
c     Input.
c
      INTEGER Ur_gcw, Ui_gcw, pr_gcw, pi_gcw, m_gcw
      INTEGER ilower0, iupper0
      INTEGER ilower1, iupper1


      REAL U0r(SIDE2d0(ilower,iupper,Ur_gcw))
      REAL U1r(SIDE2d1(ilower,iupper,Ur_gcw))
      REAL U0i(SIDE2d0(ilower,iupper,Ui_gcw))
      REAL U1i(SIDE2d1(ilower,iupper,Ui_gcw))
      REAL pr(CELL2d(ilower,iupper,pr_gcw))
      REAL pi(CELL2d(ilower,iupper,pi_gcw))
      REAL c0
      REAL dx(0:NDIM-1)

c
c     Input/Output.
c
      REAL m(CELL2d(ilower,iupper,m_gcw))

c     Local variables

      REAL fac0, fac1, fac2, g1, g2
      INTEGER i0,i1

      fac0 = 1.d0/dx(0)
      fac1 = 1.d0/dx(1) 
      fac2 = 0.5d0/c0**2

c
c     Compute the coupling source term.
c

      do i1 = ilower1,iupper1
        do i0 = ilower0, iupper0

          g1  = (pr(i0+1,i1) + pr(i0,i1))*U0r(i0+1,i1) - 
     &          (pr(i0-1,i1) + pr(i0,i1))*U0r(i0,i1)   +
     &          (pi(i0+1,i1) + pi(i0,i1))*U0i(i0+1,i1) - 
     &          (pi(i0-1,i1) + pi(i0,i1))*U0i(i0,i1)     
          g1 = fac0*fac2*g1

          g2 = (pr(i0,i1+1) + pr(i0,i1))*U1r(i0,i1+1) - 
     &         (pr(i0,i1) + pr(i0,i1-1))*U1r(i0,i1)   + 
     &         (pi(i0,i1+1) + pi(i0,i1))*U1i(i0,i1+1) - 
     &         (pi(i0,i1) + pi(i0,i1-1))*U1i(i0,i1)
          g2 = fac1*fac2*g2

          m(i0,i1) = 0.5d0*(g1+g2) 

        enddo  
      enddo
      
c     
      return
      end

