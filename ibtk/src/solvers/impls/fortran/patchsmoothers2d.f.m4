c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2022 by the IBAMR developers
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
c     Wavefront traversal of a box.
c
c     WAVEFRONT_2D(update,lower0,upper0,lower1,upper1) applies the
c     cell-update macro named by update to every cell of the box, in the
c     form update(j0,j1).  Four consecutive rows are updated per step,
c     each row lagging the previous one by one column, so that every
c     cell reads the same neighbor values as in a lexicographic sweep
c     (rows in order of increasing i1, each row in order of increasing
c     i0).  Boxes narrower than four columns and the rows left over
c     after the last complete block of four rows are swept row by row.
c     The four updates of a step do not depend on one another.  They are
c     written out rather than looped over so that the compiler keeps
c     them as four separate sequences of operations.
c
c     The caller declares the integers i0, i1, j0, j1, r and i1end.
c
define(WAVEFRONT_2D,`if (($3)-($2)+1 .ge. 4) then
         i1end = ($4) + 4*((($5)-($4)+1)/4) - 1
      else
         i1end = ($4) - 1
      endif

      do i1 = $4,i1end,4
c        Ramp up.
         do i0 = $2,($2)+2
            do r = 0,i0-($2)
               j0 = i0 - r
               j1 = i1 + r
               $1(j0,j1)
            enddo
         enddo
c        All four rows active.
         do i0 = ($2)+3,$3
            j0 = i0
            j1 = i1
            $1(j0,j1)
            j0 = i0 - 1
            j1 = i1 + 1
            $1(j0,j1)
            j0 = i0 - 2
            j1 = i1 + 2
            $1(j0,j1)
            j0 = i0 - 3
            j1 = i1 + 3
            $1(j0,j1)
         enddo
c        Ramp down.
         do i0 = ($3)+1,($3)+3
            do r = i0-($3),3
               j0 = i0 - r
               j1 = i1 + r
               $1(j0,j1)
            enddo
         enddo
      enddo

      do i1 = i1end+1,$5
         do i0 = $2,$3
            $1(i0,i1)
         enddo
      enddo')dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = D div grad U +
c     C U. Both D and C coefficients are constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(CONST_DC_UPDATE,`U($1,$2) = fac*(
     &           fac0*(U($1-1,$2)+U($1+1,$2)) +
     &           fac1*(U($1,$2-1)+U($1,$2+1)) -
     &           F($1,$2))')dnl
      subroutine smooth_gs_const_dc_2d(
     &     U,U_gcw,
     &     D,C,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw

      REAL D,C

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL    fac0,fac1,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*C)

      WAVEFRONT_2D(`CONST_DC_UPDATE',ilower0,iupper0,ilower1,iupper1)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = D
c     div grad U + C U. Both D and C coefficients
c     are constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_const_dc_2d(
     &     U,U_gcw,
     &     D,C,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw
      INTEGER red_or_black

      REAL D,C

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*C)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( mod(i0+i1,2) .eq. red_or_black ) then
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = alpha div grad U +
c     beta U with masking of certain degrees of freedom.
c
c     NOTE: The solution U is unmodified at masked degrees of freedom.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(GS_MASK_UPDATE,`if (mask($1,$2) .eq. 0) then
               U($1,$2) = fac*(
     &              fac0*(U($1-1,$2)+U($1+1,$2)) +
     &              fac1*(U($1,$2-1)+U($1,$2+1)) -
     &              F($1,$2))
            endif')dnl
      subroutine gssmoothmask2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     mask,mask_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,mask_gcw

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      INTEGER mask(ilower0-mask_gcw:iupper0+mask_gcw,
     &             ilower1-mask_gcw:iupper1+mask_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL    fac0,fac1,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      WAVEFRONT_2D(`GS_MASK_UPDATE',ilower0,iupper0,ilower1,iupper1)
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = alpha
c     div grad U + beta U with masking of certain degrees of freedom.
c
c     NOTE: The solution U is unmodified at masked degrees of freedom.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rbgssmoothmask2d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     mask,mask_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,mask_gcw
      INTEGER red_or_black

      REAL alpha,beta

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      INTEGER mask(ilower0-mask_gcw:iupper0+mask_gcw,
     &             ilower1-mask_gcw:iupper1+mask_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac = 0.5d0/(fac0+fac1-0.5d0*beta)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( (mod(i0+i1,2) .eq. red_or_black) .and.
     &           (mask(i0,i1) .eq. 0) ) then
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c              Variable coefficient patch smoothers

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = D div grad U +
c     C U. C is a cell-centered variable and
c     D is constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(CONST_D_VAR_C_UPDATE,`fac = 0.5d0/(fac0+fac1-0.5d0*C($1,$2))
            U($1,$2) = fac*(
     &           fac0*(U($1-1,$2)+U($1+1,$2)) +
     &           fac1*(U($1,$2-1)+U($1,$2+1)) -
     &           F($1,$2))')dnl
      subroutine smooth_gs_const_d_var_c_2d(
     &     U,U_gcw,
     &     D,
     &     C,C_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,C_gcw

      REAL D
      REAL C(CELL2d(ilower,iupper,C_gcw))

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL    fac0,fac1,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))


      WAVEFRONT_2D(`CONST_D_VAR_C_UPDATE',ilower0,iupper0,ilower1,
         iupper1)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = D
c     div grad U + C U. C is cell-centered variable and
c     D is constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_const_d_var_c_2d(
     &     U,U_gcw,
     &     D,
     &     C, C_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,C_gcw
      INTEGER red_or_black

      REAL D
      REAL C(CELL2d(ilower,iupper,C_gcw))

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    fac0,fac1,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( mod(i0+i1,2) .eq. red_or_black ) then
              fac = 0.5d0/(fac0+fac1-0.5d0*C(i0,i1))
               U(i0,i1) = fac*(
     &              fac0*(U(i0-1,i1)+U(i0+1,i1)) +
     &              fac1*(U(i0,i1-1)+U(i0,i1+1)) -
     &              F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = div D grad U +
c     C U.
c
c     The smoother is written for cell-centered U and side-centered
c     D = (D0,D1) with constant C coefficient.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(VAR_D_CONST_C_UPDATE,`facu0 = D0($1+1,$2)/(hx*hx)
            facl0 = D0($1,$2)/(hx*hx)
            facu1 = D1($1,$2+1)/(hy*hy)
            facl1 = D1($1,$2)/(hy*hy)
            fac   = 1.d0/(facu0+facl0+facu1+facl1-C)
            U($1,$2) = fac*(
     &           facu0*U($1+1,$2) +
     &           facl0*U($1-1,$2) +
     &           facu1*U($1,$2+1) +
     &           facl1*U($1,$2-1) -
     &           F($1,$2))')dnl
      subroutine smooth_gs_var_d_const_c_2d(
     &     U,U_gcw,
     &     D0,D1,D_gcw,
     &     C,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,D_gcw

      REAL C

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL D0(SIDE2d0(ilower,iupper,D_gcw))
      REAL D1(SIDE2d1(ilower,iupper,D_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL    hx,hy
      REAL    facu0,facl0
      REAL    facu1,facl1
      REAL    fac
c
c     Perform a single Gauss-Seidel sweep.
c
      hx = dx(0)
      hy = dx(1)

      WAVEFRONT_2D(`VAR_D_CONST_C_UPDATE',ilower0,iupper0,ilower1,
         iupper1)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F =
c     div D grad U + C U.
c
c     The smoother is written for cell-centered U and side-centered
c     D = (D0,D1) with constant C coefficient.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_var_d_const_c_2d(
     &     U,U_gcw,
     &     D0,D1,D_gcw,
     &     C,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,D_gcw
      INTEGER red_or_black
      REAL C

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL D0(SIDE2d0(ilower,iupper,D_gcw))
      REAL D1(SIDE2d1(ilower,iupper,D_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    hx,hy
      REAL    facu0,facl0
      REAL    facu1,facl1
      REAL    fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      hx = dx(0)
      hy = dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( mod(i0+i1,2) .eq. red_or_black ) then
               facu0 = D0(i0+1,i1)/(hx*hx)
               facl0 = D0(i0,i1)/(hx*hx)
               facu1 = D1(i0,i1+1)/(hy*hy)
               facl1 = D1(i0,i1)/(hy*hy)
               fac   = 1.d0/(facu0+facl0+facu1+facl1-C)
               U(i0,i1) = fac*(
     &             facu0*U(i0+1,i1) +
     &             facl0*U(i0-1,i1) +
     &             facu1*U(i0,i1+1) +
     &             facl1*U(i0,i1-1) -
     &             F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = div D grad U +
c     C U.
c
c     The smoother is written for cell-centered U, side-centered
c     D = (D0,D1) and cell-centered C coefficient.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(VAR_DC_UPDATE,`facu0 = D0($1+1,$2)/(hx*hx)
            facl0 = D0($1,$2)/(hx*hx)
            facu1 = D1($1,$2+1)/(hy*hy)
            facl1 = D1($1,$2)/(hy*hy)
            fac   = 1.d0/(facu0+facl0+facu1+facl1-C($1,$2))
            U($1,$2) = fac*(
     &           facu0*U($1+1,$2) +
     &           facl0*U($1-1,$2) +
     &           facu1*U($1,$2+1) +
     &           facl1*U($1,$2-1) -
     &           F($1,$2))')dnl
      subroutine smooth_gs_var_dc_2d(
     &     U,U_gcw,
     &     D0,D1,D_gcw,
     &     C,C_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,D_gcw,C_gcw

      REAL C(CELL2d(ilower,iupper,C_gcw))

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL D0(SIDE2d0(ilower,iupper,D_gcw))
      REAL D1(SIDE2d1(ilower,iupper,D_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL    hx,hy
      REAL    facu0,facl0
      REAL    facu1,facl1
      REAL    fac
c
c     Perform a single Gauss-Seidel sweep.
c
      hx = dx(0)
      hy = dx(1)

      WAVEFRONT_2D(`VAR_DC_UPDATE',ilower0,iupper0,ilower1,iupper1)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F =
c     div D grad U + C U.
c
c     The smoother is written for cell-centered U, side-centered
c     D = (D0,D1) and cell-centered C coefficient.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_var_dc_2d(
     &     U,U_gcw,
     &     D0,D1,D_gcw,
     &     C,C_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER U_gcw,F_gcw,D_gcw,C_gcw
      INTEGER red_or_black
      REAL C(CELL2d(ilower,iupper,C_gcw))

      REAL F(ilower0-F_gcw:iupper0+F_gcw,
     &       ilower1-F_gcw:iupper1+F_gcw)

      REAL D0(SIDE2d0(ilower,iupper,D_gcw))
      REAL D1(SIDE2d1(ilower,iupper,D_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL U(ilower0-U_gcw:iupper0+U_gcw,
     &       ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      INTEGER i0,i1
      REAL    hx,hy
      REAL    facu0,facl0
      REAL    facu1,facl1
      REAL    fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      hx = dx(0)
      hy = dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            if ( mod(i0+i1,2) .eq. red_or_black ) then
               facu0 = D0(i0+1,i1)/(hx*hx)
               facl0 = D0(i0,i1)/(hx*hx)
               facu1 = D1(i0,i1+1)/(hy*hy)
               facl1 = D1(i0,i1)/(hy*hy)
               fac   = 1.d0/(facu0+facl0+facu1+facl1-C(i0,i1))
               U(i0,i1) = fac*(
     &             facu0*U(i0+1,i1) +
     &             facl0*U(i0-1,i1) +
     &             facu1*U(i0,i1+1) +
     &             facl1*U(i0,i1-1) -
     &             F(i0,i1))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single Gauss-Seidel sweep for
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1).
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
c
c     Each component of u is updated in turn; the other components are
c     not modified while one component is updated.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(VC_UPDATE_U0,`c = beta
         if (var_c .eq. 1) then
            c = c0($1,$2)*beta
         endif

         if (use_harmonic_interp .eq. 1) then
            mu_upper = h_avg4(mu($1,$2),mu($1+1,$2),
     &                       mu($1,$2+1),mu($1+1,$2+1))
            mu_lower = h_avg4(mu($1,$2),mu($1-1,$2),
     &                       mu($1,$2+1),mu($1-1,$2+1))
         else
            mu_upper = a_avg4(mu($1,$2),mu($1+1,$2),
     &                       mu($1,$2+1),mu($1+1,$2+1))
            mu_lower = a_avg4(mu($1,$2),mu($1-1,$2),
     &                       mu($1,$2+1),mu($1-1,$2+1))
         endif

         dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &       fac1**2.d0*(mu($1,$2+1) + mu($1,$2))) - c

         nmr = -f0($1,$2) + alpha*(fac*(
     &     mu_upper*u0($1+1,$2) + mu_lower*u0($1-1,$2))+
     &     fac1**2.d0*(mu($1,$2+1)*u0($1,$2+1) + mu($1,$2)*u0($1,$2-1))+
     &     fac0*fac1*(mu($1,$2+1)*(u1($1,$2+1)-u1($1-1,$2+1))-
     &        mu($1,$2)*(u1($1,$2)-u1($1-1,$2))))

         u0($1,$2) = nmr/dnr')dnl
define(VC_UPDATE_U1,`c = beta
         if (var_c .eq. 1) then
            c = c1($1,$2)*beta
         endif

         if (use_harmonic_interp .eq. 1) then
            mu_upper = h_avg4(mu($1,$2),mu($1+1,$2),
     &                       mu($1,$2+1),mu($1+1,$2+1))
            mu_lower = h_avg4(mu($1,$2),mu($1+1,$2),
     &                       mu($1,$2-1),mu($1+1,$2-1))
         else
            mu_upper = a_avg4(mu($1,$2),mu($1+1,$2),
     &                       mu($1,$2+1),mu($1+1,$2+1))
            mu_lower = a_avg4(mu($1,$2),mu($1+1,$2),
     &                       mu($1,$2-1),mu($1+1,$2-1))
         endif

         dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &      fac0**2.d0*(mu($1+1,$2) + mu($1,$2))) - c

         nmr = -f1($1,$2) + alpha*(fac*(
     &    mu_upper*u1($1,$2+1) + mu_lower*u1($1,$2-1))+
     &    fac0**2.d0*(mu($1+1,$2)*u1($1+1,$2) + mu($1,$2)*u1($1-1,$2))+
     &    fac0*fac1*(mu($1+1,$2)*(u0($1+1,$2)-u0($1+1,$2-1))-
     &       mu($1,$2)*(u0($1,$2)-u0($1,$2-1))))

         u1($1,$2) = nmr/dnr')dnl
      subroutine vcgssmooth2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha, beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      REAL a_avg4, h_avg4
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw
      INTEGER var_c,use_harmonic_interp

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))

      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      WAVEFRONT_2D(`VC_UPDATE_U0',ilower0,iupper0+1,ilower1,iupper1)

      fac = 2.d0*fac1**2.d0
      WAVEFRONT_2D(`VC_UPDATE_U1',ilower0,iupper0,ilower1,iupper1+1)

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single "red" or "black" Gauss-Seidel sweep for
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1).
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcrbgssmooth2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp,
     &     red_or_black)
c
      implicit none
c
c     Functions.
c
      REAL a_avg4, h_avg4
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw
      INTEGER var_c,use_harmonic_interp
      INTEGER red_or_black

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))

      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single "red" or "black"  Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            if (mod(i0+i1,2) .eq. red_or_black) then

            c = beta
            if (var_c .eq. 1) then
               c = c0(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg4(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg4(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            endif

            dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &         fac1**2.d0*(mu(i0,i1+1) + mu(i0,i1))) - c

            nmr = -f0(i0,i1) + alpha*(fac*(
     &         mu_upper*u0(i0+1,i1) + mu_lower*u0(i0-1,i1))+
     &         fac1**2.d0*(mu(i0,i1+1)*u0(i0,i1+1)+
     &            mu(i0,i1)*u0(i0,i1-1))+
     &         fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-
     &            u1(i0-1,i1+1))-
     &            mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1))))

            u0(i0,i1) = nmr/dnr

          endif
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            if (mod(i0+i1,2) .eq. red_or_black) then

            c = beta
            if (var_c .eq. 1) then
               c = c1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            else
                mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            endif

            dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &         fac0**2.d0*(mu(i0+1,i1) + mu(i0,i1))) - c

            nmr = -f1(i0,i1) + alpha*(fac*(
     &         mu_upper*u1(i0,i1+1) + mu_lower*u1(i0,i1-1))+
     &         fac0**2.d0*(mu(i0+1,i1)*u1(i0+1,i1)+
     &            mu(i0,i1)*u1(i0-1,i1))+
     &         fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)
     &            -u0(i0+1,i1-1))-
     &         mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1))))

            u1(i0,i1) = nmr/dnr

           endif
         enddo
      enddo

c
      return
      end
c

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single Gauss-Seidel sweep for
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1),
c  with masking of certain degrees of freedom.
c
c     NOTE: The solution (u0,u1) is unmodified at masked degrees of freedom.
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
c
c     Each component of u is updated in turn; the other components are
c     not modified while one component is updated.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
define(VC_MASK_UPDATE_U0,`if (mask0($1,$2) .eq. 0) then

            c = beta
            if (var_c .eq. 1) then
               c = c0($1,$2)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu($1,$2),mu($1+1,$2),
     &                           mu($1,$2+1),mu($1+1,$2+1))
                mu_lower = h_avg4(mu($1,$2),mu($1-1,$2),
     &                           mu($1,$2+1),mu($1-1,$2+1))
            else
                mu_upper = a_avg4(mu($1,$2),mu($1+1,$2),
     &                           mu($1,$2+1),mu($1+1,$2+1))
                mu_lower = a_avg4(mu($1,$2),mu($1-1,$2),
     &                           mu($1,$2+1),mu($1-1,$2+1))
            endif

            dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &         fac1**2.d0*(mu($1,$2+1) + mu($1,$2))) - c

            nmr = -f0($1,$2) + alpha*(fac*(
     &         mu_upper*u0($1+1,$2) + mu_lower*u0($1-1,$2))+
     &         fac1**2.d0*(mu($1,$2+1)*u0($1,$2+1)+
     &            mu($1,$2)*u0($1,$2-1))+
     &         fac0*fac1*(mu($1,$2+1)*(u1($1,$2+1)-
     &            u1($1-1,$2+1))-
     &            mu($1,$2)*(u1($1,$2)-u1($1-1,$2))))

            u0($1,$2) = nmr/dnr

            endif')dnl
define(VC_MASK_UPDATE_U1,`if (mask1($1,$2) .eq. 0) then

            c = beta
            if (var_c .eq. 1) then
               c = c1($1,$2)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu($1,$2),mu($1+1,$2),
     &                           mu($1,$2+1),mu($1+1,$2+1))
                mu_lower = h_avg4(mu($1,$2),mu($1+1,$2),
     &                           mu($1,$2-1),mu($1+1,$2-1))
            else
                mu_upper = a_avg4(mu($1,$2),mu($1+1,$2),
     &                           mu($1,$2+1),mu($1+1,$2+1))
                mu_lower = a_avg4(mu($1,$2),mu($1+1,$2),
     &                           mu($1,$2-1),mu($1+1,$2-1))
            endif

            dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &         fac0**2.d0*(mu($1+1,$2) + mu($1,$2))) - c

            nmr = -f1($1,$2) + alpha*(fac*(
     &         mu_upper*u1($1,$2+1) + mu_lower*u1($1,$2-1))+
     &         fac0**2.d0*(mu($1+1,$2)*u1($1+1,$2)+
     &            mu($1,$2)*u1($1-1,$2))+
     &         fac0*fac1*(mu($1+1,$2)*(u0($1+1,$2)-
     &            u0($1+1,$2-1))-
     &         mu($1,$2)*(u0($1,$2)-u0($1,$2-1))))

            u1($1,$2) = nmr/dnr

            endif')dnl
      subroutine vcgssmoothmask2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     mask0,mask1,mask_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha, beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      REAL a_avg4, h_avg4
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw,mask_gcw
      INTEGER var_c,use_harmonic_interp

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))

      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      INTEGER mask0(SIDE2d0(ilower,iupper,mask_gcw))
      INTEGER mask1(SIDE2d1(ilower,iupper,mask_gcw))

      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,j0,j1,r,i1end
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      WAVEFRONT_2D(`VC_MASK_UPDATE_U0',ilower0,iupper0+1,ilower1,
         iupper1)

      fac = 2.d0*fac1**2.d0
      WAVEFRONT_2D(`VC_MASK_UPDATE_U1',ilower0,iupper0,ilower1,
         iupper1+1)

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single "red" or "black" Gauss-Seidel sweep for
c     (f0,f1) = alpha div mu (grad (u0,u1) + grad (u0, u1)^T) + beta c (u0,u1),
c  with masking of certain degrees of freedom.
c
c     NOTE: The solution (u0,u1) is unmodified at masked degrees of freedom.
c
c  The smoother is written for side-centered vector fields (u0, u1) and (f0, f1)
c  with node-centered coefficient mu and side-centered coefficient (c0,c1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcrbgssmoothmask2d(
     &     u0,u1,u_gcw,
     &     f0,f1,f_gcw,
     &     mask0,mask1,mask_gcw,
     &     c0,c1,c_gcw,
     &     mu,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_c,
     &     use_harmonic_interp,
     &     red_or_black)
c
      implicit none
c
c     Functions.
c
      REAL a_avg4, h_avg4
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER u_gcw,f_gcw,c_gcw,mu_gcw,mask_gcw
      INTEGER var_c,use_harmonic_interp
      INTEGER red_or_black

      REAL alpha,beta

      REAL mu(NODE2d(ilower,iupper,mu_gcw))

      REAL f0(SIDE2d0(ilower,iupper,f_gcw))
      REAL f1(SIDE2d1(ilower,iupper,f_gcw))

      INTEGER mask0(SIDE2d0(ilower,iupper,mask_gcw))
      INTEGER mask1(SIDE2d1(ilower,iupper,mask_gcw))

      REAL c0(SIDE2d0(ilower,iupper,c_gcw))
      REAL c1(SIDE2d1(ilower,iupper,c_gcw))

      REAL dx(0:NDIM-1)

c
c     Input/Output.
c

      REAL u0(SIDE2d0(ilower,iupper,u_gcw))
      REAL u1(SIDE2d1(ilower,iupper,u_gcw))
c
c     Local variables.
c
      INTEGER i0,i1
      REAL fac0,fac1,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single "red" or "black"  Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      fac = 2.d0*fac0**2.d0
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            if ( (mod(i0+i1,2) .eq. red_or_black) .and.
     &           (mask0(i0,i1) .eq. 0) ) then

            c = beta
            if (var_c .eq. 1) then
               c = c0(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg4(mu(i0,i1),mu(i0-1,i1),
     &                           mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg4(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            endif

            dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &         fac1**2.d0*(mu(i0,i1+1) + mu(i0,i1))) - c

            nmr = -f0(i0,i1) + alpha*(fac*(
     &         mu_upper*u0(i0+1,i1) + mu_lower*u0(i0-1,i1))+
     &         fac1**2.d0*(mu(i0,i1+1)*u0(i0,i1+1)+
     &            mu(i0,i1)*u0(i0,i1-1))+
     &         fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-
     &            u1(i0-1,i1+1))-
     &            mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1))))

            u0(i0,i1) = nmr/dnr

          endif
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            if ( (mod(i0+i1,2) .eq. red_or_black) .and.
     &           (mask1(i0,i1) .eq. 0) ) then

            c = beta
            if (var_c .eq. 1) then
               c = c1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            else
                mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                           mu(i0,i1-1),mu(i0+1,i1-1))
            endif

            dnr = alpha*(fac*(mu_upper+ mu_lower)+
     &         fac0**2.d0*(mu(i0+1,i1) + mu(i0,i1))) - c

            nmr = -f1(i0,i1) + alpha*(fac*(
     &         mu_upper*u1(i0,i1+1) + mu_lower*u1(i0,i1-1))+
     &         fac0**2.d0*(mu(i0+1,i1)*u1(i0+1,i1)+
     &            mu(i0,i1)*u1(i0-1,i1))+
     &         fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)
     &            -u0(i0+1,i1-1))-
     &         mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1))))

            u1(i0,i1) = nmr/dnr

           endif
         enddo
      enddo

c
      return
      end
c

c
