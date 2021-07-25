


module QRBM
  use debugging
  use matrices
  

  implicit none

  private :: logfile, compute_theta, wf, compute_VarOp
  public :: init,delete, SR

  character*20 :: logfile="outlog.txt"


  type :: RBM
     integer :: Nv, Nh, dmn
     double precision ::  spin
     double complex, allocatable :: a(:),b(:)
     double complex, allocatable :: W(:,:)
  end type RBM


contains

  subroutine init(self, Nv, Nh, spin)
    implicit none
    type(RBM):: self
    integer, intent(in) :: Nv, Nh
    integer ::  status, ii, jj
    double precision, intent(in) :: spin
    double precision ::  re, im
    logical :: ca = .false., cb = .false., cW = .false.

    
    self%Nv = Nv
    self%Nh = Nh
    self%dmn = Nv + Nh +Nv*Nh
    self%spin = spin

    allocate(self%a(Nv), stat=status)
    if (status .eq. 0) ca = .true.
       
    allocate(self%b(Nh), stat = status)
    if (status .eq. 0) cb = .true.

    allocate(self%W(Nv, Nh), stat = status)
    if (status .eq. 0) cW = .true.

    if (ca .and. cb .and. cW .eqv. .true.) then
       print *, "RBM succesfully initialized"
    else
       print *, "RBM allocation: FAILED"
       stop
    endif

    do ii=1,Nv
       call random_number(re)
       call random_number(im)
       self%a(ii) = complex(re*0.2-0.1,im*0.2-0.1)
    enddo

    do jj=1,Nh
       call random_number(re)
       call random_number(im)
       self%b(jj) = complex(re*0.2-0.1,im*0.2-0.1)
    enddo

    do jj=1,Nh
       do ii=1,Nv
          call random_number(re)
          call random_number(im)
          self%W(ii,jj) = complex(re*0.2-0.1,im*0.2-0.1)
       enddo
    enddo
    
    
  end subroutine init


  subroutine delete(self)
    implicit none
    type(RBM) :: self
    deallocate(self%a, self%b, self%W)

    print *, "RBM succesfully deleted"

  end subroutine delete

  subroutine compute_theta(self, S, theta)
    implicit none
    type(RBM), intent(in) :: self
    double complex, intent(in) :: S(self%Nv)
    double complex, intent(out) :: theta(self%Nh)

    theta = self%b+matmul(S, self%W)

  end subroutine compute_theta

  double complex function wf(self,S)
    implicit none
    type(RBM), intent(in) :: self
    integer :: jj
    double complex, dimension(self%Nv), intent(in) :: S
    double complex, dimension(self%Nh) :: theta, act
    

    ! Compute theta
    call compute_theta(self, S, theta)

    !print *, "theta", double precision(sum((conjg(theta)*theta)**2))

    ! Compute activation
    act = cosh(theta)

    !print *, "act", double precision(sum((conjg(act)*act)**2))

    ! Add exponential part
    wf = exp(dot_product(self%a, S))

   
    

    ! Add product of cosh
    do jj=1,self%Nh
       wf = wf*act(jj) ! Cut the 2 if numbers are to high
    enddo

    

  end function wf

  double complex function Heisenberg_loc(self,S, pbc)
    implicit none
    type(RBM), intent(in) :: self
    integer :: ii, jj
    double complex:: frac
    double complex, dimension(self%Nv), intent(in) :: S
    double complex, dimension(self%Nh) :: theta, theta_new, act, act_new
    logical, intent(in) :: pbc

    ! Initialize local energy
    Heisenberg_loc = complex(0d0,0d0)

    ! Commute theta for S
    call compute_theta(self, S, theta)

    ! Bulk term
    do  ii=1,self%Nv-1

       !Compute theta of new configuration
       theta_new = theta - 2.*self%W(ii,:)*S(ii)-2.*self%W(ii+1,:)*S(ii+1)

       ! Compute activations
       act = cosh(theta)
       act_new = cosh(theta_new)

       ! Compute the fraction between exp part of wave functions
       frac = exp(-2.*S(ii)*self%a(ii)-2.*S(ii+1)*self%a(ii+1))

       ! Multiply by cosh 
       do jj=1,self%Nh
          frac = frac*act_new(jj)/act(jj)
       enddo

       ! Update local energy
       Heisenberg_loc = Heisenberg_loc + frac*(1-S(ii)*S(ii+1)) + S(ii)*S(ii+1)

    enddo

    ! Periodic boundary conditions

    if (pbc .eqv. .true.) then
       ! Copmute theta of new configuration
       theta_new = theta - 2.*self%W(1,:)*S(1)-2.*self%W(self%Nv,:)*S(self%Nv)

       ! Compute activations
       act = cosh(theta)
       act_new = cosh(theta_new)

       ! Compute the fraction between exp part of wave functions
       frac = exp(-2.*S(1)*self%a(1)-2.*S(self%Nv)*self%a(self%Nv))

       ! Multiply by cosh 
       do jj=1,self%Nh
          frac = frac*act_new(jj)/act(jj)
       enddo

       ! Update local energy
       Heisenberg_loc = Heisenberg_loc + frac*(1-S(1)*S(self%Nv)) + S(1)*S(self%Nv)
       
    endif

  end function Heisenberg_loc

  double complex function Ising_loc(self, S, h, pbc)
    implicit none
    type(RBM), intent(in) :: self
    integer :: ii, jj
    double precision, intent(in) :: h
    double complex:: frac
    double complex, dimension(self%Nv), intent(in) :: S
    double complex, dimension(self%Nh) :: theta, theta_new, act, act_new
    logical, intent(in) :: pbc

    ! Initialize local energy
    Ising_loc = complex(0d0,0d0)

    ! Compute theta for S
    call compute_theta(self, S, theta)

    ! Linear term
    do  ii=1,self%Nv

       !Compute theta of new configuration
       theta_new = theta - 2.*self%W(ii,:)*S(ii)

       ! Compute activations
       act = cosh(theta)
       act_new = cosh(theta_new)

       ! Compute the fraction between exp part of wave functions
       frac = exp(-2.*S(ii)*self%a(ii))

       ! Multiply by cosh 
       do jj=1,self%Nh
          frac = frac*act_new(jj)/act(jj)
       enddo

       ! Update local energy
       Ising_loc = Ising_loc - frac

    enddo

    ! Add field
    Ising_loc = Ising_loc*h

    ! Interaction term
    do ii=1,self%Nv-1
       Ising_loc = Ising_loc -S(ii)*S(ii+1)
    end do
   

    ! Periodic boundary conditions

    if (pbc .eqv. .true.) then
       Ising_loc = Ising_loc - S(1)*S(self%Nv)
    endif

  end function Ising_loc
  
 

    

  double precision function sample_wf(self, S, N_sweeps)
    implicit none
    type(RBM), intent(in) :: self
     !f2py3 intent(in,out) :: self
    integer, intent(in) :: N_sweeps
    integer :: ii,jj, site
    double precision :: uu, prob
    double complex:: frac
    double complex, dimension(self%Nv):: S 
    double complex, dimension(self%Nh) :: theta, theta_new, act, act_new

    do ii=1,self%Nv*N_sweeps

       ! Compute theta for old configuration
       call compute_theta(self, S, theta)

       ! Generate a site in (1,.., Nv)
       call random_number(uu)
       site = int(floor(self%Nv*uu)+1.)
       

       ! Update  theta
       theta_new = theta - 2.*self%W(site,:)*S(site)

       ! Compute activations
       act = cosh(theta)
       act_new = cosh(theta_new)

       ! Compute the fraction between exp part of wave functions
       frac = exp(-2.*S(site)*self%a(site))

       ! Multiply by cosh 
       do jj=1,self%Nh
          frac = frac*act_new(jj)/act(jj)
       enddo

       ! Compute the probability
       prob = min(1., dble(conjg(frac)*frac))

    

       ! Metropolis algorithm
       call random_number(uu)
       if (uu.le.prob) then
          S(site) = -1.*S(site)
       endif

    enddo
    
       
    ! Compute the square wave function of sampled configuration
    sample_wf = dble(conjg(wf(self,S))*wf(self,S))

  end function sample_wf

  
  subroutine compute_VarOp(self, S, O)
    implicit none
    type(RBM), intent(in) :: self
     !f2py3 intent(in,out) :: self
    integer :: ii,jj
    double complex, intent(in) :: S(self%Nv)
    double complex :: Oh(self%Nh)
    double complex, intent(out) :: O(self%dmn)

    ! Copmute variational operator
    O(1:self%Nv) = S
    Oh = tanh(self%b+matmul(S,self%W))
    O(self%Nv+1:self%Nv+self%Nh) = Oh

    do jj=1,self%Nh
       do ii=1,self%Nv
          O(self%Nv+self%Nh+(jj-1)*self%Nv+ii) = S(ii)*Oh(jj)
       enddo
    enddo

  end subroutine compute_VarOp

  subroutine SR(self, model, h, pbc, steps, N_sweeps, N_samples, E, debug)
    implicit none
    type(RBM) :: self
     !f2py3 intent(in,out) :: self
    integer :: step, sample, ii,jj
    integer, intent(in) :: model, steps, N_sweeps, N_samples
    double precision :: norm, square_wf, lr, uu, threshold, start, finish
    double precision, intent(in) :: h
    double precision, intent(out) :: E(steps)
    double complex:: Eloc_av, Eloc, S(self%Nv)
    double complex, dimension(self%dmn) :: O, Odag, O_av, Odag_av, ElocOdag, ElocOdag_av, Force
    double complex, dimension(self%dmn, self%dmn) :: OdagO,OdagO_av, t_dot, Fisher, Inverse
    logical, intent(in) :: pbc, debug

    threshold = 0.00001

    open(9, file=logfile, status="unknown")

    if(debug) then
       write(9,*) "Logfile for RBM ansatz"
       write(9,*) "    "
       if(model.eq.0) then
          write(9,*) "Model: Heisenberg"
          write(9,*) "Parameters:"
          write(9,*) "Number of visible spins:   ", self%Nv
          write(9,*) "Number of hidden spins:   ", self%Nh
          write(9,*) "Number of variational parameters:   ", self%dmn
          write(9,*) "Periodic boundary conditions:   ", pbc
          write(9,*) "Number  of SR iterations:   ", steps
          write(9,*) "Number of MC sweep through the lattice:   ", N_sweeps
          write(9,*) "Number of samples:   ", N_samples
       elseif(model.eq.1) then
          write(9,*) "Model: Ising"
          write(9,*) "Parameters:"
          write(9,*) "Number of spins:   ", self%Nv
          write(9,*) "Number of hidden spins:   ", self%Nh
          write(9,*) "Number of variational parameters:   ", self%dmn
          write(9,*) "H_field:   ", h
          write(9,*) "Periodic boundary conditions:   ", pbc
          write(9,*) "Number  of SR iterations:   ", steps
          write(9,*) "Number of MC sweep through the lattice:   ", N_sweeps
          write(9,*) "Number of samples:   ", N_samples
       else
          write(9,*) "Unrecognize model"
          write(9,*) "Computation stopped."
          stop
       endif

       write(9,*) "*************************************************************************************************"

       write(9,*) "Starting Sochastic Reconfiguration..."
    endif

    do step=1,steps

       if(debug) then
          write(9,*) "Epoch: ", step
       end if

       ! Initialize tensors
       OdagO_av = complex(0d0,0d0)
       Odag_av = complex(0d0,0d0)
       O_av = complex(0d0,0d0)
       ElocOdag_av = complex(0d0,0d0)
       Eloc_av = complex(0d0,0d0)
       norm = 0.0

       !Generate the standard sequence
       do ii=1,self%Nv
          call random_number(uu)
          uu = (2.*floor(2*uu)-1.)*self%spin
          S(ii) = uu
       enddo

       if(debug) then
          write(9,*) "Square a:   ", dble(sum(conjg(self%a)*self%a))
          write(9,*) "Square b:   ", dble(sum(conjg(self%b)*self%b))
          write(9,*) "Square W:   ", dble(sum(conjg(self%W)*self%W))
       endif
     
       ! Sample sequence and probability according to MCMC
       square_wf = sample_wf(self, S, 100)

       if(debug) then
          write(9,*)   "-------------------------------------------------------------------------------------------------"
          write(9,*) "Starting sampling..."
       endif

       call cpu_time(start)
       do sample=1,N_samples

         
          
          ! Sample sequence and probability according to MCMC
          square_wf = sample_wf(self, S, N_sweeps)

         

          ! Compute O and Odag
          call compute_VarOp(self, S, O)
          Odag = conjg(O)

          !print *, "VarOp", O

          ! Compute OdagO
          do jj=1,self%dmn
             do ii=1,self%dmn
                OdagO(ii,jj) = Odag(ii)*O(jj)
             enddo
          enddo

          ! Compute Eloc
          if (model .eq. 0) then
             Eloc = Heisenberg_loc(self,S, pbc)
          elseif (model .eq. 1) then
             Eloc = Ising_loc(self, S, h, pbc)
          else
             print *, "Undefined model"
             stop
          endif
          
             

          ! Copmute ElocOdag
          ElocOdag = Eloc*Odag

          ! Update averages
          norm = norm+square_wf
          ElocOdag_av = ElocOdag_av+ElocOdag*square_wf
          Eloc_av =  Eloc_av+Eloc*square_wf

          OdagO_av = OdagO_av + OdagO*square_wf
          O_av = O_av +O*square_wf
          Odag_av = Odag_av +Odag*square_wf

       enddo

       call cpu_time(finish)

        if(debug) then
          write(9,*) "Sampling recap:"
          write(9,*) "Sampling time (s):   ", finish-start
          write(9,*) "Last sampled square_wf:   ", square_wf
          write(9,*) "Normalization constant:   ", norm
          write(9,*) "End sampling."
          write(9,*)  "-------------------------------------------------------------------------------------------------"
       endif

      



       !Compute t_dot
       do jj=1,self%dmn
          do ii=1,self%dmn
             t_dot(ii,jj) = Odag_av(ii)*O_av(jj)
          enddo
       enddo

       ! Compute Fisher matrix
       Fisher = OdagO_av/norm-t_dot/norm**2


       if(debug) then
          write(9,*) "Sum squared Fisher matrix:   ", sum(conjg(Fisher)*Fisher)
          !Check hermiticity
          call check_hermiticity(Fisher, threshold, logfile)
       endif

       
       

       ! Compute learning rate
       lr =max(0.01*(0.9**step), 0.001)


       ! Regularize Fisher matrix
       do ii=1,self%dmn
          Fisher(ii,ii) = Fisher(ii,ii)+ 0.0001
       enddo


       !call print_matrix(Fisher)

       if(debug) then
          write(9,*) "Sum squared regularized Fisher matrix:   ", sum(conJG(Fisher)*Fisher)
       endif

       !print *, "Fisher", Fisher
       

       call cpu_time(start)
       Inverse = Fisher
       ! Invert Fisher matrix
       call inv_herm(Inverse)

       call cpu_time(finish)

       if(debug) then
          write(9,*) "Inversion time (s):   ", finish-start
          ! Check inverse
          call check_inverse(Fisher, Inverse, threshold, logfile)
       endif

  

       ! Compute the force
       Force = ElocOdag_av/norm-Eloc_av*Odag_av/norm**2

       if(debug) then
          write(9,*) "Sum squared Force:   ", sum(conjg(Force)*Force)
          write(9,*) "Sum squared Inverse:   ", sum(conjg(Inverse)*Inverse)
       end if
      
       
      

       ! Compute new weights (save in Force)
       Force = lr*matmul(Inverse, Force)

       if(debug) then
          write(9,*) "Sum squared learning force:   ", sum(conjg(Force)*Force)
          !write(9,*) "Learning force:   ", Force
       endif

      
      
       ! Update weights and biases
       self%a = self%a -Force(1:self%Nv)
       self%b = self%b - Force(self%Nv+1: self%Nv+self%Nh)

     

       do jj=1,self%Nh
          do ii=1,self%Nv
             self%W(ii,jj) = self%W(ii,jj) - Force(self%Nv+self%Nh+(jj-1)*self%Nv+ii)
          enddo
       enddo

       if(debug) then
          write(9,*) "Expected energy:   ", Eloc_av/norm
       endif
       

       E(step) = dble(Eloc_av/norm)

       if(debug) then
          write(9,*) "End of epoch:   ", step, "."
          write(9,*)  "*************************************************************************************************"
       endif
      

      

    enddo

  end subroutine SR
  
    

end module QRBM


  
     
  

   

