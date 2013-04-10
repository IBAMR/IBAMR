
      SUBROUTINE loop15_F90(N, x, a, b, c)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N)

      x = c + a*b
      RETURN
      END


      SUBROUTINE loop15_F90Overhead(N, x, a, b, c)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), c(N)

      RETURN
      END
