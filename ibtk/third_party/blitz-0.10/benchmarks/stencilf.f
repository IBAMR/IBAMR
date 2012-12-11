       subroutine stencilf(A, B, N, iters)
       integer N, iters
       double precision A(N,N,N), B(N,N,N)
       integer i,j,k,z
       double precision c
       c = 1 / 7.

       do z=1,iters
       
         do k=2,N-1
           do j=2,N-1
             do i=2,N-1
               A(i,j,k) = c * (B(i,j,k) + B(i+1,j,k) + B(i-1,j,k) 
     .            + B(i,j+1,k) + B(i,j-1,k) + B(i,j,k+1) + B(i,j,k-1))
             enddo
           enddo
         enddo

         do k=2,N-1
           do j=2,N-1
             do i=2,N-1
               B(i,j,k) = c * (A(i,j,k) + A(i+1,j,k) + A(i-1,j,k) 
     .            + A(i,j+1,k) + A(i,j-1,k) + A(i,j,k+1) + A(i,j,k-1))
             enddo
           enddo
         enddo

       enddo

       return
       end

