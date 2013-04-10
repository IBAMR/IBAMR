
      PROGRAM CFDF

      PARAMETER         ( N = 50, niters = 10 )

      REAL VX(N,N,N), VY(N,N,N), VZ(N,N,N)
      REAL VX2(N,N,N), VY2(N,N,N), VZ2(N,N,N)
      REAL P(N,N,N)
      REAL FX(N,N,N), FY(N,N,N), FZ(N,N,N)
      REAL AX(N,N,N), AY(N,N,N), AZ(N,N,N)

      INTEGER iter

C     Initialize arrays

      CALL initialize(N, VX, VY, VZ, VX2, VY2, VZ2, P, FX, FY, FZ,
     .  AX, AY, AZ)
   
C     Apply the stencil a few times
 
      DO iter=1,niters
        CALL cfdStencil(N, VX, VY, VZ, VX2, VY2, VZ2, P, 
     .    FX, FY, FZ, AX, AY, AZ)
      END DO 
       
      STOP
      END

      SUBROUTINE cfdStencil(N, VX, VY, VZ, VX2, VY2, VZ2, P,
     .  FX, FY, FZ, AX, AY, AZ)

      INTEGER N
      REAL VX(N,N,N), VY(N,N,N), VZ(N,N,N)
      REAL VX2(N,N,N), VY2(N,N,N), VZ2(N,N,N)
      REAL P(N,N,N)
      REAL FX(N,N,N), FY(N,N,N), FZ(N,N,N)
      REAL AX(N,N,N), AY(N,N,N), AZ(N,N,N)

      PARAMETER ( delta_t = 0.001, recip_rho = 1.0e-3,
     .                    eta = 1.0e-6, c1 = 0.1, c2 = 0.1 )

      DO i=3,N-2
       DO j=3,N-2
        DO k=3,N-2
         VX2(i,j,k) = VX(i,j,k)+delta_t*(recip_rho*(eta*
     .    c1 * (-90*VX(i,j,k)-VX(i-2,j,k)+16*VX(i-1,j,k)
     .    +16*VX(i+1,j,k)-VX(i+2,j,k)-VX(i,j-2,k)+16*VX(i,j-1,k)
     .    +16*VX(i,j+1,k)-VX(i,j+2,k)-VX(i,j,k-2)+16*VX(i,j,k-1)
     .    +16*VX(i,j,k+1)-VX(i,j,k+2))+c2*(P(i-2,j,k)
     .    -8*P(i-1,j,k)+8*P(i+1,j,k)+P(i+2,j,k))+FX(i,j,k))
     .    -AX(i,j,k))
         VY2(i,j,k) = VY(i,j,k)+delta_t*(recip_rho*(eta*
     .    c1 * (-90*VY(i,j,k)-VY(i-2,j,k)+16*VY(i-1,j,k)
     .    +16*VY(i+1,j,k)-VY(i+2,j,k)-VY(i,j-2,k)+16*VY(i,j-1,k)
     .    +16*VY(i,j+1,k)-VY(i,j+2,k)-VY(i,j,k-2)+16*VY(i,j,k-1)
     .    +16*VY(i,j,k+1)-VY(i,j,k+2))+c2*(P(i,j-2,k)
     .    -8*P(i,j-1,k)+8*P(i,j+1,k)+P(i,j+2,k))+FY(i,j,k))
     .    -AY(i,j,k))
         VZ2(i,j,k) = VZ(i,j,k)+delta_t*(recip_rho*(eta*
     .    c1 * (-90*VZ(i,j,k)-VZ(i-2,j,k)+16*VZ(i-1,j,k)
     .    +16*VZ(i+1,j,k)-VZ(i+2,j,k)-VZ(i,j-2,k)+16*VZ(i,j-1,k)
     .    +16*VZ(i,j+1,k)-VZ(i,j+2,k)-VZ(i,j,k-2)+16*VZ(i,j,k-1)
     .    +16*VZ(i,j,k+1)-VZ(i,j,k+2))+c2*(P(i,j,k-2)
     .    -8*P(i,j,k-1)+8*P(i,j,k+1)+P(i,j,k+2))+FZ(i,j,k))
     .    -AZ(i,j,k))
        END DO
       END DO
      END DO

      RETURN
      END



      SUBROUTINE initialize(N, VX, VY, VZ, VX2, VY2, VZ2, P,
     .  FX, FY, FZ, AX, AY, AZ)

      INTEGER N
      REAL VX(N,N,N), VY(N,N,N), VZ(N,N,N)
      REAL VX2(N,N,N), VY2(N,N,N), VZ2(N,N,N)
      REAL P(N,N,N)
      REAL FX(N,N,N), FY(N,N,N), FZ(N,N,N)
      REAL AX(N,N,N), AY(N,N,N), AZ(N,N,N)

      DO i=3,N-2
       DO j=3,N-2
        DO k=3,N-2
         VX(i,j,k) = 0
         VY(i,j,k) = 0
         VZ(i,j,k) = 0
         P(i,j,k) = 0
         FX(i,j,k) = 0
         FY(i,j,k) = 0
         FZ(i,j,k) = 0
         AX(i,j,k) = 0
         AY(i,j,k) = 0
         AZ(i,j,k) = 0
        END DO
       END DO
      END DO

      END

