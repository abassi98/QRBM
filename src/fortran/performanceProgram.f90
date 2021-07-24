   program main
      use debugging
      use QRBM
      use symmQRBM
      use matrices
      
      implicit none
      type(RBM) :: Net
      type(sRBM) :: sNet
      integer :: Nv,Nh,N_sweeps,N_samples, steps, model, alpha, NN, N_max,N_min 
      integer, allocatable :: seed(:)
      double precision ::  h, start, finish, sstart, sfinish, f, s
      double precision, dimension(:), allocatable :: E
      logical :: pbc, debug


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
      

      steps = 1000
      N_samples = 100
      N_sweeps = 10
      model = 1 !Ising
      h = 0.0
      alpha = 1
      pbc = .true.
      debug = .false.

      allocate(E(steps))

      ! Open performance file
      open(10,file="data/performance.txt", status="unknown", access="append")

      
      

      call cpu_time(s)

      do Nv=N_min,N_max
         Nh = int(Nv*alpha)
         
         !Init Net
         call init(Net, Nv, Nh, spin= dble(1.0))
         !SR
         call cpu_time(start)
         call  SR(Net, model, h, pbc, steps, N_sweeps, N_samples, E, debug)
         call cpu_time(finish)
         call delete(Net)

         !Init symmetric  Net
         call s_init(sNet, Nv, alpha, spin=dble(1.0))
         !SR
         call cpu_time(sstart)
         call s_SR(sNet, model, h, pbc, steps, N_sweeps, N_samples, E, debug)
         call cpu_time(sfinish)
         call s_delete(sNet)

         write(10,*) Nv,  finish-start, sfinish-sstart

      enddo

      call cpu_time(f)

      close(10)

      print *, "Total computation time (s):", f-s

      
      
    end program main
    
    
    
    

      
