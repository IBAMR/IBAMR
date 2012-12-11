      REAL*8 AA(513,513), DD(513,513)
      REAL*8 X (513,513), Y (513,513)
      REAL*8 RX(513,513), RY(513,513)
      DO J = 2,N-1
         DO I = 2,N-1
            XX = X(I+1,J)-X(I-1,J)
            YX = Y(I+1,J)-Y(I-1,J)
            XY = X(I,J+1)-X(I,J-1)
            YY = Y(I,J+1)-Y(I,J-1)
            A = 0.25 * (XY*XY+YY*YY)
            B = 0.25* (XX*XX+YX*YX)
            C = 0.125 * (XX*XY+YX*YY)
            AA(I,J) = -B
            DD(I,J) = B+B+A*REL
            PXX = X(I+1,J)-2.*X(I,J)+X(I-1,J)
            QXX = Y(I+1,J)-2.*Y(I,J)+Y(I-1,J)
            PYY = X(I,J+1)-2.*X(I,J)+X(I,J-1)
            QYY = Y(I,J+1)-2.*Y(I,J)+Y(I,J-1)
            PXY = X(I+1,J+1)-X(I+1,J-1)-X(I-1,J+1)+X(I-1,J-1)
            QXY = Y(I+1,J+1)-Y(I+1,J-1)-Y(I-1,J+1)+Y(I-1,J-1)
            RX(I,J) = A*PXX+B*PYY-C*PXY
            RY(I,J) = A*QXX+B*QYY-C*QXY
         END DO
      END DO
