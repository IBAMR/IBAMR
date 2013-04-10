!
! Tuned Fortran 77 version
! Optimizations:
! - Rather than four separate arrays, one 4D array is allocated.  This
!   allows P1, P2, P3, and C to be interlaced in memory, improving data
!   locality
! - The stencil is tiled to improve cache usage
! - Instead of copying, the indices into the 4D array are shuffled.

!      INTEGER N, iters
!      REAL check

!      N = 112
!      iters = 210
!      CALL acoustic3d_f77Tuned(N,iters,check)
!      PRINT *, check
!      END

      SUBROUTINE acoustic3d_f77Tuned(N, niters, check)
      INTEGER N, niters, iter
      REAL check
      REAL P1(N,N,N), P2(N,N,N), P3(N,N,N), C(N,N,N)

      CALL acoust3d_f77Tuned_setup(P1,P2,P3,C,N)

      DO iter = 1, niters, 3
        CALL acoustic3d_f77Tuned_stencil(P1,P2,P3,C,N)
        CALL acoustic3d_f77Tuned_stencil(P2,P3,P1,C,N)
        CALL acoustic3d_f77Tuned_stencil(P3,P1,P2,C,N)
      END DO

      check = P1(N/2,N/2,N/2)

      RETURN
      END



      SUBROUTINE acoustic3d_f77Tuned_stencil(P1,P2,P3,C,N)
      INTEGER N
      REAL P1(N,N,N), P2(N,N,N), P3(N,N,N), C(N,N,N)
      INTEGER i,j,k
      INTEGER bi,bj,bk,ni,nj,nk,blockSize

      blockSize = 8

      DO bk=2,N-1,blockSize
        nk = min(bk+blockSize-1,N-1)
        DO bj=2,N-1,blockSize
          nj = min(bj+blockSize-1,N-1)
          DO bi=2,N-1,blockSize
            ni = min(bi+blockSize-1,N-1)
            DO k=bk,nk
              DO j=bj,nj
                DO i=bi,ni
                  P3(i,j,k) = (2-6*C(i,j,k))*P2(i,j,k) 
     .              + C(i,j,k) * (P2(i,j-1,k) + P2(i,j+1,k) 
     .              + P2(i-1,j,k) + P2(i+1,j,k) + P2(i,j,k-1)
     .              + P2(i,j,k+1)) - P1(i,j,k)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

      RETURN
      END




      SUBROUTINE acoust3d_f77Tuned_setup(P1, P2, P3, c, N)
      INTEGER N
      REAL P1(N,N,N), P2(N,N,N), P3(N,N,N), c(N,N,N)
      INTEGER i, j, k
      INTEGER cavityLeft, cavityRight, cavityFront, cavityBack,
     .        cavityTop, cavityBottom, cavityTop2, cavityBottom2
      REAL ci, cj, ck, s2

!     Set up velocity field

      DO k=1,N
        DO j=1,N
          DO i=1,N/2
            c(i,j,k) = 0.05
          END DO
          DO i=N/2+1,N
            c(i,j,k) = 0.3
          END DO
        END DO
      END DO

!     Cavities

      cavityLeft = 3 * N / 7.0 
      cavityRight = 4 * N / 7.0 
      cavityFront = 3 * N / 7.0
      cavityBack = 4 * N / 7.0
      cavityTop = 5 * N / 7.0
      cavityBottom = 6 * N / 7.0
      cavityTop2 = 1 * N / 7.0
      cavityBottom2 = 2 * N / 7.0

      DO k=cavityFront, cavityBack
        DO j=cavityLeft, cavityRight
          DO i=cavityTop, cavityBottom
            c(i,j,k) = 0.02
          END DO
          DO i=cavityTop2, cavityBottom2
            c(i,j,k) = 0.001
          END DO
        END DO
      END DO

!     Initial pressure distribution

      ci = N/2.0
      cj = N/2.0
      ck = N/2.0
      s2 = 64.0 * 9.0 / ((N/2.0) ** 2)
  
      DO k=1,N
        DO j=1,N
          DO i=1,N
            P1(i,j,k) = 0.0
            P2(i,j,k) = exp(- ((i-ci)**2+(j-cj)**2+(k-ck)**2) * s2)
            P3(i,j,k) = 0.0
          END DO
        END DO
      END DO

      CALL acoust3d_f77Tuned_setup_check(P2, N)
      CALL acoust3d_f77Tuned_setup_check(c, N)

      RETURN
      END



      SUBROUTINE acoust3d_f77Tuned_setup_check(P, N)
      INTEGER N
      REAL P(N,N,N)
      INTEGER i, j, k
      REAL check

      check = 0.0
      DO k=1,N
        DO j=1,N
          DO i=1,N
            check = check + P(i,j,k) * (i+N*j+N*N*k)
          END DO
        END DO
      END DO
     
      PRINT *, 'Array check: ', check 

      RETURN
      END

