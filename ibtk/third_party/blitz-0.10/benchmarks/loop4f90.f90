      SUBROUTINE loop4_F90(niters, N, x, a, b)
      INTEGER i, iter, niters, N
      REAL*8 x(N), a(N), b(N)

      DO iter=1,niters
        x=a+b
      END DO

      RETURN
      END
