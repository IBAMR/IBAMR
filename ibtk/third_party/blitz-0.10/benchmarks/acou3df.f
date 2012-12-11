!      INTEGER N, iters
!      REAL check

!      N = 112
!      iters = 210
!      CALL acoustic3d_f77(N,iters,check)
!      PRINT *, check
!      END

      SUBROUTINE acoustic3d_f77(N, niters, check)
      INTEGER N, niters, iter
      REAL check
      REAL P1(N,N,N), P2(N,N,N), P3(N,N,N), C(N,N,N)
      INTEGER i, j, k

      CALL acoust3d_f77_setup(P1, P2, P3, C, N)

      DO iter=1, niters
        DO k=2,N-1
          DO j=2,N-1
            DO i=2,N-1
              P3(i,j,k) = (2-6*C(i,j,k))*P2(i,j,k) + C(i,j,k)
     .          * (P2(i,j-1,k) + P2(i,j+1,k) + P2(i-1,j,k) 
     .           + P2(i+1,j,k) + P2(i,j,k-1) + P2(i,j,k+1)) 
     .           - P1(i,j,k)
            END DO
          END DO
        END DO

        DO k=1,N
          DO j=1,N
            DO i=1,N
              P1(i,j,k) = P2(i,j,k)
              P2(i,j,k) = P3(i,j,k)
            END DO
          END DO
        END DO
      END DO

      check = P1(N/2,N/2,N/2)

      RETURN
      END



      SUBROUTINE acoust3d_f77_setup(P1, P2, P3, c, N)
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

      CALL acoust3d_f77_setup_check(P2, N)
      CALL acoust3d_f77_setup_check(c, N)

      RETURN
      END



      SUBROUTINE acoust3d_f77_setup_check(A, N)
      INTEGER N
      REAL A(N,N,N)
      INTEGER i, j, k
      REAL check

      check = 0.0
      DO k=1,N
        DO j=1,N
          DO i=1,N
            check = check + A(i,j,k) * (i+N*j+N*N*k)
          END DO
        END DO
      END DO
     
      PRINT *, 'Array check: ', check 

      RETURN
      END

