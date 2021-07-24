
module debugging
         contains
         !Every subroutine has an error counter 
        
         !initialize a debug variable. If it is true, then the debu gqill be performed
         subroutine init_debug(debug,printer)
            implicit none
            character :: yes, var 
            logical :: debug,printer !our debugging variable and printer
            printer=.false. !Default no warmings messages
            yes="y"
            print *, "Do you want to start the debugging mode?[y]/[n]"
            read *, var
            if (var.eq.yes) then
               debug=.true.
               print *, "=================="
               print *, "Debugging mode: ON"
               print *, "=================="
               print *, "Do you want  to print checkpoints?[y]/[n]"
               read *, var
               if(var.eq.yes) printer=.true.
            else
               debug=.false.
               print *, "==================="
               print *, "Debugging mode: OFF"
               print *, "==================="
            endif
            
          end subroutine init_debug
          

      !check if the dimensiions are positive
         subroutine check_dimensions(nn,mm,debug,printer,count)
            implicit none
            integer :: nn,mm !dimensions of the matrix
            logical :: debug,printer
            integer :: count !error counter
            if(debug.eqv..true.) then
               if (nn.ge.1 .and. mm.ge.1) then
                  if(printer.eqv..true.) then
                  print *,  "=================================="
                  print *,  "Check dimensions: OK"
                  print *,  "=================================="
                  endif
               else
                  if(printer.eqv..true.) then
                  print *, "====================================="
                  print *, "Check dimensions of: ERROR"
                  print *, "====================================="
                  endif
                  nn=1 !set default  values for dimensions if it encounters some error
                  mm=1 !try to interpret dimensions anyway
                  count=count+1
                  stop
               endif
            endif
          end subroutine check_dimensions
          
      
     
    
      
      !Check the allocation status
         subroutine check_allocation(stat,debug,printer)
            implicit none
            integer :: stat     !it takes and integr value. 0 is the allocation succeeded
            logical :: debug,printer
            if(debug.eqv..true.) then
               if(stat.ne.0) then
                  if(printer.eqv..true.) then
                  print *, "Allocation:  FAILED"
                  print *, "Program aborted"
               endif
               stop !if allocation fails, the program must be stopped
               else
                  if(printer.eqv..true.) then
                     print *, "Allocation: SUCCEEDED"
                  endif
               endif
            endif
          end subroutine check_allocation
          

!check eigenvalues  computation status
      subroutine check_eigen(info,debug,printer,count)
      implicit none
      integer :: info
      logical :: debug,printer
      integer :: count !error counter
      if(debug.eqv..true.) then
         if(info.eq.0) then
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: SUCCEEDED"
            endif
         else
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: FAILED"
               print *, "INFO=",info
            endif
            count=count+1
         endif
      endif
      end subroutine check_eigen
    
         !This subroutine checks if the computation time exceeds a certain threshold.
         subroutine check_time(start,threshold,debug,printer)
            implicit none
            logical :: debug,printer
            double precision :: start, finish,threshold
            if(debug.eqv..true.) then
               call cpu_time(finish)
               if(finish-start .gt. threshold) then
                  if(printer.eqv..true.) then
                  print *, "Computation time exceeded"
                  print *,"Ran in: ",finish-start
                  print *, "Program aborted"
                  endif
                  stop
               endif
            endif
          end subroutine check_time

         subroutine check_diff(nn,mm,mat1,mat2,threshold,debug,printer,count)
            implicit none
            integer :: nn,mm
            double complex :: mat1(nn,mm),mat2(nn,mm)
            logical :: debug,printer
            double precision :: threshold,diff
            integer :: ii,jj,temp,count
            temp=0
            if(debug.eqv..true.) then
               do jj=1,mm
                  do ii=1,nn
                     diff=sqrt(dble((mat1(ii,jj)-mat2(ii,jj))*conjg(mat1(ii,jj)-mat2(ii,jj))))
                     if(diff.gt.threshold) then
                        temp=temp+1
                     endif    
                  enddo
               enddo
               if(temp.eq.0) then
                  if(printer.eqv..true.) then
                     print *, "Check equality: OK"
                  endif         
               else
                  if(printer.eqv..true.) then
                     print *, "Check equality: ERROR"
                  endif
                  count=count+1
               endif 
            endif 
          end subroutine check_diff


          subroutine check_hermiticity(A, threshold, logfile)
            implicit none
            double complex :: A(:,:)
            double precision :: threshold
            integer :: NN, unit
            character*20 :: logfile

            NN = size(A,1)

            inquire(file=logfile, number=unit)
            if(unit .eq. -1) then
               print *, "Error: logfile should be open"
               stop
            endif
            

            if (real(sum((A-transpose(conjg(A)))**2)).ge.NN**2*threshold**2) then
               write(unit,*)  "Matrix is noy hermitian: stop computation."
               stop
            else
               write(unit,*) "Matrix is hermitian."
            endif
            

          end subroutine check_hermiticity

          subroutine check_inverse(A, Inverse, threshold, logfile)
            implicit none
            double complex :: A(:,:), Inverse(:,:), prod(size(A,1), size(A,2))
            double precision :: threshold
            integer :: NN, ii, unit
            character*20 :: logfile


            inquire(file=logfile, number=unit)
            if(unit .eq. -1) then
               print *, "Error: logfile should be open"
               stop
            endif

            
            NN = size(A,1)

            prod = matmul(A, Inverse)

            do ii=1,NN
               prod(ii,ii) = prod(ii,ii) - complex(1d0,0d0)
            enddo

            if (dble(sum(prod**2)).ge.NN**2*threshold**2) then
               write(unit,*) "AA^(-1) product failed"
               stop
            else
               write(unit,*) "AA^(-1) product succeded. Continue..."
            endif
            

            prod = matmul(Inverse,A)

            do ii=1,NN
               prod(ii,ii) = prod(ii,ii) - complex(1d0,0d0)
            enddo

            if (dble(sum(prod**2)).ge.NN**2*threshold**2) then
               write(unit,*) "A^(-1)A product failed"
               stop
            else
               write(unit,*) "A^(-1)A product succeded. Continue..."
            endif
            
           

          end subroutine check_inverse
          
          
            
            
      
         
      end module 


      module matrices
      contains

        subroutine inv_herm(A)
          implicit none
          double complex, intent(inout) :: A(:,:)
          double complex ::  AP(size(A,1)*(size(A,1)+1)/2)
          integer :: NN, ii,jj, info
      
          ! Define thesize of the matrix
          NN = size(A,1)

          ! Compute AP UPLO = 'U'
          do jj = 1,NN
             do ii=1,jj
                AP(ii+(jj-1)*jj/2) = A(ii,jj)
             enddo
          enddo
          
          !Compute LU factorization
          call zpptrf('U', NN, AP, info)

          !Check singularity
          if (info.ne.0)  stop "LU decomposition failed: matrix numerically singular"
    

          !Compute inverse
          call zpptri('U',NN,AP,info)

          ! Compute Inverse Inverse
          do jj = 1,NN
             do ii=1,jj
                A(ii,jj) = AP(ii+(jj-1)*jj/2)
                A(jj,ii) = conjg(A(ii,jj))
             enddo
          enddo

          !Check inversion
          if (info.ne.0) stop "Inversion failed: matrix numerically singular"

        end subroutine inv_herm


        subroutine print_matrix(A)
          implicit none
          integer :: NN, MM, ii, jj
          double complex, intent(in) :: A(:,:)

          NN = size(A,1)
          MM = size(A,2)

          do ii=1,NN
             print *, (A(ii,jj), jj=1,MM)
          enddo
        end subroutine print_matrix
        
          

          
      end module matrices
      
          
          
          
