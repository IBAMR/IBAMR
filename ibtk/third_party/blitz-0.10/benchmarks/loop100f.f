
      SUBROUTINE loop100_F77(N, a, b, c, d, x, u, v, w, z)
      INTEGER i, N
      REAL*8 a(N), b(N), c(N), d(N), x(N), u, v, w, z

      DO i=1,N
          x(i)=(1.0-c(i)*c(i))/((4*w)*sin(1.0+c(i)*c(i)-2*v*c(i)))*a(i
     !)*b(i)*u*exp(-z*d(i));
      END DO
      RETURN
      END


      SUBROUTINE loop100_F77Overhead(N, a, b, c, d, x, u, v, w, z)
      INTEGER i, N
      REAL*8 a(N), b(N), c(N), d(N), x(N), u, v, w, z
      RETURN
      END
