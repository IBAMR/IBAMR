      SUBROUTINE f90daxpy(a, x, y, n, iters)
      INTEGER n, iters
      DOUBLE PRECISION, DIMENSION (n) :: x, y
      DOUBLE PRECISION a, b

      b = - a

      DO i=1,iters
        y = y + a * x
        y = y + b * x
      END DO
      END SUBROUTINE
