   program main
      use debugging
      use QRBM
      use matrices
      
      implicit none
      type(RBM) :: Net
      integer :: Nv,Nh,N_sweeps,N_samples, steps, model, alpha, N_av, NN,ii,jj
      integer, allocatable :: seed(:)
      double precision ::  h, start, finish, E_exp, std
      double precision, dimension(:), allocatable :: E
      logical :: pbc, debug


      ! Seeding
      call random_seed(size=NN)
      allocate(seed(NN))
      seed = 10
      call random_seed(put=seed)


      !Read parameters
      print *, "Insert number of spins:"
      read *, Nv

      print *, "Insert alpha:"
      read *, alpha

      print *, "Insert number of steps:"
      read *, steps

      print *, "Insert model (0 -> Heisenberg, 1-> Ising):"
      read *, model
      if(model.eq.1) then
         print *, "Insert h_field:"
         read *, h
      else
         h = 0d0
      endif

      N_sweeps = 10
      N_samples = 100
      pbc = .true.


      ! Write parameters in a file for plots
      open(10, file="temp/parameters.txt", status="unknown")
      write(10,*) Nv,alpha, steps, model, h, N_sweeps, N_samples
      close(10)
      
      
      ! Read debug to print checks on logfile
      print *, "Do you want to debug?(True/False):"
      read *, debug

      Nh=int(Nv*alpha)
      
      ! Initialize the network
      call init(Net, Nv, Nh, spin= dble(1.0))

      ! Allocate energy vector
      allocate(E(steps))

      

      call cpu_time(start)
      ! Stochastic reconfiguration
      call SR(Net, model, h, pbc, steps, N_sweeps, N_samples, E, debug)

      call cpu_time(finish)


      N_av = 100
      ! Compute the average  and std of last N_av steps
      E_exp = sum(E(steps-N_av+1:steps))/N_av
      std = sqrt(sum(E(steps-N_av+1:steps)**2)/N_av-E_exp**2)
      
         

      print *, "Estimated ground state energy and std: ", E_exp, std
      
      
      
     

      open(12, file="data/ground_state.txt", status='unknown')
      write(12,*) E
      close(12)

      open(12,file="data/a.txt", status="unknown")
      write(12,*) Net%a
      close(12)

      open(12,file="data/b.txt", status="unknown")
      write(12,*) Net%b
      close(12)

      open(12,file="data/W.txt", status="unknown")
      do ii=1,Nv
         write(12,*) (Net%W(ii,jj), jj=1,Nh)
      enddo
      
      close(12)

      call delete(Net)


      
      
      print *, "RBM computation time (s): ", finish-start
      

      
      
    end program main
    
    
    
    

      
