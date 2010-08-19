c
c     Copyright (c) 2002-2010 Boyce Griffith
c
c     Permission is hereby granted, free of charge, to any person
c     obtaining a copy of this software and associated documentation
c     files (the "Software"), to deal in the Software without
c     restriction, including without limitation the rights to use, copy,
c     modify, merge, publish, distribute, sublicense, and/or sell copies
c     of the Software, and to permit persons to whom the Software is
c     furnished to do so, subject to the following conditions:
c
c     The above copyright notice and this permission notice shall be
c     included in all copies or substantial portions of the Software.
c
c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
c     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
c     DEALINGS IN THE SOFTWARE.
c
dnl Process this file with m4 to produce FORTRAN source code
define(NDIM,1)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Detect sharp spatial gradients.
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine advect_detectgrad1d(
     &  ifirst0,ilast0,
     &  vghost0,tagghost0,ttagghost0,
     &  dx,
     &  gradtol,
     &  dotag,donttag,
     &  var,
     &  tags,temptags)
c
      implicit none
include(TOP_SRCDIR/src/fortran/const.i)dnl
c
      INTEGER
     &  ifirst0,ilast0,
     &  dotag,donttag,
     &  vghost0,
     &  tagghost0,
     &  ttagghost0
      REAL
     &  dx(0:NDIM-1),
     &  gradtol
      REAL
     &  var(CELL1dVECG(ifirst,ilast,vghost))
      INTEGER
     &  tags(CELL1dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL1dVECG(ifirst,ilast,ttagghost))
c
      REAL tol
      REAL facejump, loctol
      REAL presm1,presp1
      logical tagcell
      INTEGER ic0
c
      tol = gradtol
      
      do ic0=ifirst0,ilast0
         
         if (tags(ic0) .ne. 0) then
            loctol = 0.125*tol
         else 
            loctol = tol
         endif
         
         tagcell = .false.
         
         presm1 = var(ic0-1)
         presp1 = var(ic0+1)
         facejump = abs(var(ic0)-presm1)
         facejump = max(facejump,abs(var(ic0)-presp1))
         tagcell = ((facejump).gt.(loctol*dx(0)))
         
         if ( tagcell ) then
            temptags(ic0) = dotag
         endif
      enddo
c     
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
