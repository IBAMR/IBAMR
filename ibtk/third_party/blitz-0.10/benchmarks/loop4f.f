      SUBROUTINE loop4_F77(niters, N, x, a, b)
      INTEGER i, iter, niters, N
      REAL*8 x(N), a(N), b(N)

      DO iter=1,niters
        DO i=1,N
          x(i)=a(i)+b(i)
        END DO
      END DO

      RETURN
      END
