SUBROUTINE stencilf90(A, B, n, iters)
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: n, iters
    DOUBLE PRECISION, DIMENSION (n,n,n) :: A, B
    DOUBLE PRECISION :: c
    INTEGER :: count

    c = 1 / 7.

    DO count=1,iters
        A(2:N-1,2:N-1,2:N-1) = c * (B(2:N-1,2:N-1,2:N-1) + B(3:N,2:N-1,2:N-1) &
           + B(1:N-2,2:N-1,2:N-1) + B(2:N-1,3:N,2:N-1) + B(2:N-1,1:N-2,2:N-1) &
           + B(2:N-1,2:N-1,3:N) + B(2:N-1,2:N-1,1:N-2))
     
        B(2:N-1,2:N-1,2:N-1) = c * (A(2:N-1,2:N-1,2:N-1) + A(3:N,2:N-1,2:N-1) &
           + A(1:N-2,2:N-1,2:N-1) + A(2:N-1,3:N,2:N-1) + A(2:N-1,1:N-2,2:N-1) &
           + A(2:N-1,2:N-1,3:N) + A(2:N-1,2:N-1,1:N-2))
    END DO
END SUBROUTINE
