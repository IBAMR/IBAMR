c ---------------------------------------------------------------------
c
c Copyright (c) 2026 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

dnl This file only exists to be included into the cell-centered linear
dnl refine and coarsen kernels (cart_cell_refine{2,3}d.f.m4 and
dnl cart_cell_coarsen{2,3}d.f.m4); it is never compiled on its own.

dnl coarsen_index(i,i_c,i_f,ratio) sets i_c to the index of the coarse cell that
dnl contains the fine cell i and i_f to the lowest fine index of that coarse cell.
define(coarsen_index,`dnl
if ($1.lt.0) then
            $2=($1+1)/$4-1
         else
            $2=$1/$4
         endif
         $3=$2*$4
')dnl'

dnl linear_stencil(i_c,xi,lower_cell,upper_cell,clower,cupper,touches_lower,
dnl                touches_upper,axis)
dnl sets the two coarse indices idx<axis>(0:1) and the weights w<axis>(0:1)
dnl of linear interpolation to a fine point at offset xi from the center of
dnl the coarse cell i_c. At a regular physical boundary, the stencil is one
dnl sided: lower_cell and upper_cell are the first and last cells of the
dnl region that is not extended past the boundary.
define(linear_stencil,`dnl
if ($5 .eq. $6 .and.
     &          $7 .ne. 0 .and.
     &          $8 .ne. 0) then
            idx$9(0) = $1
            idx$9(1) = $1
            w$9(0) = 1.d0
            w$9(1) = 0.d0
         else if ($7 .ne. 0 .and.
     &            $1 .eq. $3 .and. $2 .le. 0.d0) then
            idx$9(0) = $1
            idx$9(1) = $1+1
            w$9(0) = 1.d0
            w$9(1) = 0.d0
         else if ($8 .ne. 0 .and.
     &            $1 .eq. $4 .and. $2 .ge. 0.d0) then
            idx$9(0) = $1-1
            idx$9(1) = $1
            w$9(0) = 0.d0
            w$9(1) = 1.d0
         else if ($2 .le. 0.d0) then
            idx$9(0) = $1-1
            idx$9(1) = $1
            w$9(0) = -$2
            w$9(1) = 1.d0+$2
         else
            idx$9(0) = $1
            idx$9(1) = $1+1
            w$9(0) = 1.d0-$2
            w$9(1) = $2
         endif
')dnl'
