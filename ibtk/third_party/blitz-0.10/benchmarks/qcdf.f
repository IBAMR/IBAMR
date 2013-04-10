C  Initial Fortran 77 version of the Lattice QCD benchmark.

       subroutine qcdf(M, res, src, V, iters)
       integer V, iters, i, site, spin, col
       complex*16 M(3,3,V), res(3,2,V), src(3,2,V)
  
       DO i=1,iters
         DO site=1,V
           DO spin=1,2
             DO col=1,3
                res(col,spin,site) = M(col,1,site) * src(1,spin,site) 
     .            + M(col,2,site) * src(2,spin,site) 
     .            + M(col,3,site) * src(3,spin,site)
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       return
       end

C  Hand-tuned version
C  Changes:
C    o Ordering of array altered to improve layout of data in memory
C    o col and spin loops unwound; it was found that unwinding the
C      col loop inside the spin loop was marginally faster (by 1.1%)
C    o Unwinding both loops was faster than unwinding just one.

       subroutine qcdf2(M, res, src, V, iters)
       integer V, iters, i, site
       complex*16 M(3,3,V), res(3,2,V), src(3,2,V)

       DO i=1,iters
         DO site=1,V

C             col=1, spin=1

              res(1,1,site) = M(1,1,site) * src(1,1,site)
     .            + M(1,2,site) * src(2,1,site)
     .            + M(1,3,site) * src(3,1,site)

C             col=2, spin=1

              res(2,1,site) = M(2,1,site) * src(1,1,site)
     .            + M(2,2,site) * src(2,1,site)
     .            + M(2,3,site) * src(3,1,site)

C             col=3, spin=1

              res(3,1,site) = M(3,1,site) * src(1,1,site)
     .            + M(3,2,site) * src(2,1,site)
     .            + M(3,3,site) * src(3,1,site)

C             col=1, spin=2

              res(1,2,site) = M(1,1,site) * src(1,2,site)
     .            + M(1,2,site) * src(2,2,site)
     .            + M(1,3,site) * src(3,2,site)

C             col=2, spin=2

              res(2,2,site) = M(2,1,site) * src(1,2,site)
     .            + M(2,2,site) * src(2,2,site)
     .            + M(2,3,site) * src(3,2,site)

C             col=3, spin=2

              res(3,2,site) = M(3,1,site) * src(1,2,site)
     .            + M(3,2,site) * src(2,2,site)
     .            + M(3,3,site) * src(3,2,site)

         ENDDO
       ENDDO

       return
       end

