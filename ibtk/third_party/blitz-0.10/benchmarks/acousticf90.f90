!INTEGER N, niters
!REAL check
!N = 128
!niters = 128*3
!CALL echo_f90(N, niters, check)
!PRINT *, check
!END

SUBROUTINE echo_f90(N, niters, check)
  INTEGER, INTENT( IN ) :: N, niters
  REAL, INTENT( OUT ) :: check

  REAL, DIMENSION (N,N) :: P1, P2, P3, c
  INTEGER iter

  CALL echo_f90_setupInitialConditions(c, P1, P2, P3, N)
  CALL checkArray_f90(P2, N)
  CALL checkArray_f90(c, N)

  DO iter=1, niters
    P3(2:N-1,2:N-1) = (2-4*c(2:N-1,2:N-1)) * P2(2:N-1,2:N-1) &
        + c(2:N-1,2:N-1)*(P2(1:N-2,2:N-1) + P2(3:N,2:N-1)    &
        + P2(2:N-1,1:N-2) + P2(2:N-1,3:N)) - P1(2:N-1,2:N-1)
    P1 = P2
    P2 = P3
  END DO

  check = P1(N/2,7*N/8)

  RETURN
END





SUBROUTINE echo_f90_setupInitialConditions(c, P1, P2, P3, N)
  INTEGER, INTENT( IN ) :: N
  REAL, DIMENSION (N,N) :: P1(N,N), P2(N,N), P3(N,N), c(N,N)

  INTEGER blockLeft, blockRight, blockTop, blockBottom
  INTEGER channelLeft, channelRight, channel1Height, channel2Height
  REAL cr, cc
  INTEGER i, j
  REAL s2

  ! Set the velocity field
  c = 0.2

  ! Solid block with which the pulse collides
  blockLeft = 1
  blockRight = 2 * N / 5.0
  blockTop = N / 3.0
  blockBottom = 2 * N / 3.0
  c(blockTop:blockBottom, blockLeft:blockRight) = 0.5

  ! Channel directing the pulse leftwards
  channelLeft = 4 * N / 5.0
  channelRight = N
  channel1Height = 3 * N / 8.0
  channel2Height = 5 * N / 8.0

  c(channel1Height,channelLeft:channelRight) = 0.0;
  c(channel2Height,channelLeft:channelRight) = 0.0;

  ! Initial pressure distribution: a gaussian pulse inside the channel
  cr = N / 2.0
  cc = 7.0 * N / 8.0
  s2 = 64.0 * 9.0 / ((N / 2.0) ** 2)

  DO j=1,N
    DO i=1,N
      P2(i,j) = exp(-((i-cr)**2 + (j-cc)**2) * s2)
    END DO
  END DO

  P1 = 0.0
  P3 = 0.0
END





      SUBROUTINE checkArray_f90(A, N)
      INTEGER N
      REAL, DIMENSION(N,N) :: A

      INTEGER i,j
      REAL check
      check = 0.0
      DO j=1,N
        DO i=1,N
          check = check + (i*n+j)*A(i,j)
        END DO
      END DO

      PRINT *, 'Array check: ', check
      RETURN
      END

