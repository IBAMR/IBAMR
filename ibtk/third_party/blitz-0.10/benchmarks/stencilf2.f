! Tiled fortran stencil

       SUBROUTINE stencilftiled(A, B, N, iters)
       INTEGER N, iters
       DOUBLE PRECISION A(N,N,N), B(N,N,N)
       INTEGER z
       DOUBLE PRECISION c
       c = 1 / 7.

       DO z=1,iters
         CALL tiledStencil(A, B, N, c)
         CALL tiledStencil(B, A, N, c)
       ENDDO 
 
       RETURN
       END

       SUBROUTINE tiledStencil(A, B, N, c)
       INTEGER N
       DOUBLE PRECISION A(N,N,N), B(N,N,N)
       DOUBLE PRECISION c
       INTEGER i,j,k,bi,bj,bk,ni,nj,nk
       INTEGER blockSize

       blockSize = 16

       DO bi=2,N-1,blockSize
         DO bj=2,N-1,blockSize
           DO bk=2,N-1,blockSize
             ni = min(bi+blockSize-1,N-1)
             nj = min(bj+blockSize-1,N-1)
             nk = min(bk+blockSize-1,N-1)

             DO k=bk,nk
               DO j=bj,nj
                 DO i=bi,ni
                   A(i,j,k) = c * (B(i,j,k) + B(i+1,j,k) + B(i-1,j,k) 
     .            + B(i,j+1,k) + B(i,j-1,k) + B(i,j,k+1) + B(i,j,k-1))
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO

       RETURN
       END

