
      SUBROUTINE loop100_F90(N, a, b, c, d, x, u, v, w, z)
      INTEGER i, N
      REAL*8 a(N), b(N), c(N), d(N), x(N), u, v, w, z

      x=(1.0-c*c)/((4*w)*sin(1.0+c*c-2*v*c))*a*b*u*exp(-z*d)
      RETURN
      END


      SUBROUTINE loop100_F90Overhead(N, a, b, c, d, x, u, v, w, z)
      INTEGER i, N
      REAL*8 a(N), b(N), c(N), d(N), x(N), u, v, w, z

      RETURN
      END
