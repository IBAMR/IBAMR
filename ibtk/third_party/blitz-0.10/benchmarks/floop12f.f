
      SUBROUTINE floop12_F77(N, y, x, a, b, c, d, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), d(N), u

      DO i=1,N
          y(i) = u+a(i); x(i) = a(i)+b(i)+c(i)+d(i);
      END DO
      RETURN
      END


      SUBROUTINE floop12_F77Overhead(N, y, x, a, b, c, d, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), d(N), u
      RETURN
      END
