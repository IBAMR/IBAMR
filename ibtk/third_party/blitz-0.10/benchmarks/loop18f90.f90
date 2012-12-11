
      SUBROUTINE loop18_F90(N, x, a, b, u, v)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), u, v

      x = (u+a)*(v+b)
      RETURN
      END


      SUBROUTINE loop18_F90Overhead(N, x, a, b, u, v)
      INTEGER i, N
      REAL*8 x(N), a(N), b(N), u, v

      RETURN
      END
