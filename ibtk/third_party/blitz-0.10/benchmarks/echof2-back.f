! Tuned Fortran 77 version
! Optimizations:
!   - blocked stencil algorithm to improve cache use
!   - arrays interlaced by making one big 3-dimensional array
!   - copying of arrays avoided by cycling indices into the
!     3-d array

      SUBROUTINE echo_f77Tuned(N, niters, check)
      INTEGER N, niters, iter
      REAL check
      REAL A(N,N,4)
!     P1 = A(N,N,1), P2 = A(N,N,2), P3 = A(N,N,3), C = A(N,N,4)
      INTEGER P1, P2, P3, C
      INTEGER i, j
      INTEGER bi,bj,ni,nj,blockSize

      P1 = 1
      P2 = 2
      P3 = 3
      C = 4

      CALL echo_f77Tuned_setInitialConditions(A,C,P1,P2,P3,N)

      blockSize = 128

      DO iter=1, niters
        DO bj=2,N-1,blockSize
          nj = min(bj+blockSize-1,N-1)
          DO bi=2,N-1,blockSize
            ni = min(bi+blockSize-1,N-1)
            DO j=bj,nj
              DO i=bi,ni
                A(i,j,P3) = (2-4*A(i,j,C))*A(i,j,P2) + A(i,j,C)
     .           *(A(i,j-1,P2) + A(i,j+1,P2) + A(i-1,j,P2) 
     .           + A(i+1,j,P2)) - A(i,j,P1)
              END DO
            END DO
          END DO
        END DO
        P1 = P2
        P2 = P3
      END DO

      check = A(N/2,N/2,P1)

      RETURN
      END




      SUBROUTINE echo_f77Tuned_setInitialConditions(A, C, P1, P2, P3, N)
      INTEGER N
      REAL A(N,N,4)
      INTEGER C, P1, P2, P3
      INTEGER i, j, blockLeft, blockRight, blockTop, blockBottom
      INTEGER channelLeft, channelRight, channel1Height, channel2Height
      INTEGER cr, cc
      REAL s2

!     Default velocity in the air

      DO j=1,N
        DO i=1,N
          A(i,j,C) = 0.2;
        END DO
      END DO
    
!     Solid block with which the pulse collids
 
      blockLeft = 1
      blockRight = 2 * N / 5.0
      blockTop = N / 3.0
      blockBottom = 2 * N / 3.0

      DO j=blockLeft,blockRight
        DO i=blockTop,blockBottom
          A(i,j,C) = 0.5
        END DO
      END DO

!     Channel directing the pulse leftwards

      channelLeft = 4 * N / 5.0
      channelRight = N
      channel1Height = 3 * N / 8.0
      channel2Height = 5 * N / 8.0

      DO j = channelLeft,channelRight
        A(channel1Height,j,C) = 0.0;
        A(channel2Height,j,C) = 0.0;
      END DO

!     Initial pressure distribution: a gaussian pulse inside the channel

      cr = N / 2
      cc = 7 * N / 8.0
      s2 = 64.0 * 9.0 / ((N/2.0) ** 2)

      DO j=1,N
        DO i=1,N
          A(i,j,P1) = 0.0
          A(i,j,P2) = exp(-((i-cr)**2 + (j-cc)**2) * s2);
          A(i,j,P3) = 0.0
        END DO
      END DO

      CALL checkArray2(A,P2,N)
      CALL checkArray2(A,C,N)

      RETURN
      END





      SUBROUTINE checkArray2(A, P, N)
      INTEGER N, P
      REAL A(N,N,4)

      INTEGER i,j
      REAL check
      check = 0.0
      DO j=1,N
        DO i=1,N
          check = check + (i*n+j)*A(i,j,P)
        END DO
      END DO

      PRINT *, 'Array check: ', check
      RETURN
      END
