   program main
      use debugging
      use matrices
      
      implicit none
      integer :: NN, dmn, N_max,N_min,ii,jj
      integer, allocatable :: seed(:)
      double precision ::  start, finish, f, s
      double complex, allocatable :: matrix(:,:), O(:), Odag(:)



      ! Seeding
      call random_seed(size=NN)
      allocate(seed(NN))
      seed = 10
      call random_seed(put=seed)

      ! Parameters
      print *, "Insert N_max:"
      read *, N_max
      print *, "Insert N_min:"
      read *, N_min

      if(N_min>N_max) then
         print *, "Error: Invalid limits."
         stop
      end if
      

     

      ! Open diagonalization file
      open(10,file="data/diagonalization.txt", status="unknown", access="append")

      
      

      call cpu_time(s)

      do dmn=N_min,N_max
         allocate(matrix(dmn,dmn))
         allocate(O(dmn))
         allocate(Odag(dmn))
         O = complex(1d0,2d0)
         Odag = conjg(O)
         do ii=1,dmn
            do jj=1,dmn
               matrix(ii,jj) = Odag(ii)*O(jj)
               if(ii.eq.jj) then
                  matrix(ii,ii) = matrix(ii,ii) +0.0001
               endif
               
            enddo
         enddo
      
         call cpu_time(start)
         call inv_herm(matrix)
         call cpu_time(finish)
         write(10,*) dmn, finish-start
         deallocate(matrix, O,Odag)

      enddo

      call cpu_time(f)

      close(10)

      print *, "Total computation time (s):", f-s

      
      
    end program main
    
    
    
    

      
