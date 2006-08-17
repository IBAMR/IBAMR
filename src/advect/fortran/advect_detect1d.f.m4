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
