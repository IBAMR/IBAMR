
      SUBROUTINE floop10_F77(N, x, a, b, c, u)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), u

      DO i=1,N
          x(i) = u+a(i)+b(i)+c(i);
      END DO
      RETURN
      END


      SUBROUTINE floop10_F77Overhead(N, x, a, b, c, u)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), u
      RETURN
      END
